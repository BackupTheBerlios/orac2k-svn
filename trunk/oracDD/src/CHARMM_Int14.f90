!!$/---------------------------------------------------------------------\
!!$   Copyright  © 2006-2007 Massimo Marchi <Massimo.Marchi at cea.fr>   |
!!$                                                                      |
!!$    This software is a computer program named oracDD whose            |
!!$    purpose is to simulate and model complex molecular systems.       |
!!$    The code is written in fortran 95 compliant with Technical        |
!!$    Report TR 15581, and uses MPI-1 routines for parallel             |
!!$    coding.                                                           |
!!$                                                                      |
!!$    This software is governed by the CeCILL license under             |
!!$    French law and abiding by the rules of distribution of            |
!!$    free software.  You can  use, modify and/ or redistribute         |
!!$    the software under the terms of the CeCILL icense as              |
!!$    circulated by CEA, CNRS and INRIA at the following URL            |
!!$    "http://www.cecill.info".                                         |
!!$                                                                      |
!!$    As a counterpart to the access to the source code and rights      |
!!$    to copy, modify and redistribute granted by the license,          |
!!$    users are provided only with a limited warranty and the           |
!!$    software's author, the holder of the economic rights, and         |
!!$    the successive licensors have only limited liability.             |
!!$                                                                      |
!!$    The fact that you are presently reading this means that you       |
!!$    have had knowledge of the CeCILL license and that you accept      |
!!$    its terms.                                                        |
!!$                                                                      |
!!$    You should have received a copy of the CeCILL license along       |
!!$    with this program; if not, you can collect copies on the URL's    |
!!$    "http://www.cecill.info/licences/Licence_CeCILL_V2-en.html"       |
!!$    "http://www.cecill.info/licences/Licence_CeCILL_V2-fr.html"       |
!!$                                                                      |
!!$----------------------------------------------------------------------/
SUBROUTINE Int14(fpx,fpy,fpz)
  REAL(8) :: fpx(:),fpy(:),fpz(:)

  REAL(8), POINTER :: xp0(:),yp0(:),zp0(:)
  INTEGER, SAVE :: MyCalls=0
  REAL(8), ALLOCATABLE, SAVE :: ecc6(:),ecc12(:)
  INTEGER, ALLOCATABLE, SAVE :: Id_ij(:,:)
  REAL(8), PARAMETER :: a1=0.2548296D0,a2=-0.28449674D0,a3&
       &=1.4214137D0,a4=-1.453152D0,a5=1.0614054D0,qp=0.3275911D0
  REAL(8), SAVE :: alphal,twrtpi
  INTEGER, POINTER :: Slv(:),Id(:)
  REAL(8), POINTER :: Charge(:)
  INTEGER :: i,i1,i2,j,lij,type,nint14
  REAL(8) :: xpi,ypi,zpi,xa,ya,za,rsq,rsp,rsqi,qforce,xpj,ypj,zpj&
       &,ssvir,r6,r12,chrgei,chrgej,qt,expcst,erfcst 
  REAL(8) :: rspqi,alphar,furpar,elj,ucon,ucou,ucoul(2),uconf(2)


  IF(Calls == 0) THEN
     Calls=Calls+1
     Conv_Fact=1000.0D0*4.184/(unite*avogad)
  END IF
  IF(MyCalls == 0) THEN
     MyCalls=MyCalls+1
     CALL Init
  END IF
  xp0=>Atoms(:) % x
  yp0=>Atoms(:) % y
  zp0=>Atoms(:) % z
  Id =>Atoms(:) % Id_Type
  Slv=>Atoms(:) % Id_Slv
  Charge=>Atoms(:) % chg

  nint14=SIZE(Indx_Int14,2)
  ucoul=0.0_8
  uconf=0.0_8

  DO i=1,nint14
     i1=Indx_Int14(1,i)
     i2=Indx_Int14(2,i)
     type=Slv(i1)
     lij=Id_ij(Id(i1),Id(i2))
     chrgei=charge(i1)
     chrgej=charge(i2)

     xpi=xp0(i1)
     ypi=yp0(i1)
     zpi=zp0(i1)
     xpj=xp0(i2)
     ypj=yp0(i2)
     zpj=zp0(i2)
     
     xa=xpi-xpj
     ya=ypi-ypj
     za=zpi-zpj
     rsq=xa**2+ya**2+za**2
     rsp=DSQRT(rsq)
     rsqi=1.0d0/rsq
     rspqi=rsqi/rsp
     r6=rsqi*rsqi*rsqi
     r12=r6*r6
     ssvir=12.0d0*ecc12(lij)*r12-6.0d0*ecc6(lij)*r6
     qforce=ssvir*rsqi
     ucon=(ecc12(lij)*r12-ecc6(lij)*r6)
     uconf(type)=uconf(type)+ucon
     alphar=alphal*rsp
     qt=1.0d0/(1.0d0+qp*alphar)
     expcst=dexp(-alphar*alphar)
     erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*expcst
     furpar=chrgei*chrgej
     ucou=furpar*erfcst/rsp
     ucoul(type)=ucoul(type)+ucou
     
!!$     phi(i1)=phi(i1)+chrgej*erfcst/rsp
!!$     phi(i2)=phi(i2)+chrgei*erfcst/rsp

     qforce=qforce+furpar*(erfcst+twrtpi*alphar*expcst)*rspqi
     fpx(i1)=fpx(i1)+qforce*xa
     fpy(i1)=fpy(i1)+qforce*ya
     fpz(i1)=fpz(i1)+qforce*za
     fpx(i2)=fpx(i2)-qforce*xa
     fpy(i2)=fpy(i2)-qforce*ya
     fpz(i2)=fpz(i2)-qforce*za
  END DO
  Uint14coul_slt=ucoul(1)
  Uint14coul_slv=ucoul(2)
  Uint14conf_slt=uconf(1)
  Uint14conf_slv=uconf(2)
CONTAINS
  SUBROUTINE Init
    INTEGER :: n,m,ij
    
    twrtpi=2.0d0/SQRT(pi)
    alphal = Ewald__Param % alpha
    
    n=SIZE(LennardJones__Par % Par_SE)
    ALLOCATE(Id_ij(n,n))
    ALLOCATE(ecc6(n*(n+1)/2),ecc12(n*(n+1)/2))
    DO n=1,SIZE(LennardJones__Par % Par_SE)
       DO m=n,SIZE(LennardJones__Par % Par_SE)
          ij=m*(m-1)/2+n
          Id_ij(n,m)=ij
          Id_ij(m,n)=ij
          ecc6(ij)=LennardJones__Par % c6_(ij)
          ecc12(ij)=LennardJones__Par % c12_(ij)
       END DO
    END DO
  END SUBROUTINE Init
END SUBROUTINE Int14
