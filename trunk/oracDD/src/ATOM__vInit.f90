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
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Jan 23 2009 -                                     *
!!$                                                                      *
!!$***********************************************************************
FUNCTION Atom__vInit_() RESULT(out)
  LOGICAL :: out
  REAL(8), ALLOCATABLE :: vx(:),vy(:),vz(:),xa(:),ya(:),za(:),mass(:)
  REAL(8) :: ekt,massa,tvel,sig,v1,v2,v3,tmass
  INTEGER :: n,Grp_No
  
  ekt=0.5*gascon*T % t0
  
  IF(T % Gauss_Init) THEN
     DO n=1,SIZE(Atoms)
        v1=cgauss()       
        v2=cgauss()       
        v3=cgauss()       
        massa=Atoms(n) % mass
        IF(Atoms(n) % knwn == 1 .AND. massa /=0.0_8) THEN
           tvel=Boltz*T % t0/(unite*massa)
           sig=SQRT(tvel)
           Atoms(n) % vx = v1*sig
           Atoms(n) % vy = v2*sig
           Atoms(n) % vz = v3*sig
        END IF
     END DO
  END IF
  
  ALLOCATE(vx(SIZE(IndBox_a_p)),vy(SIZE(IndBox_a_p))&
       &,vz(SIZE(IndBox_a_p)),mass(SIZE(IndBox_a_p)))
  ALLOCATE(xa(SIZE(IndBox_a_p)),ya(SIZE(IndBox_a_p))&
       &,za(SIZE(IndBox_a_p)))
  
  vx=Atoms(IndBox_a_p(:)) % vx
  vy=Atoms(IndBox_a_p(:)) % vy
  vz=Atoms(IndBox_a_p(:)) % vz
  xa=Atoms(IndBox_a_p(:)) % xa
  ya=Atoms(IndBox_a_p(:)) % ya
  za=Atoms(IndBox_a_p(:)) % za
  mass=Atoms(IndBox_a_p(:)) % mass
  xa=xa+_PBC(xa)
  ya=ya+_PBC(ya)
  za=za+_PBC(za)
  
  CALL Velocities
  
  Atoms(IndBox_a_p(:)) % vx=vx
  Atoms(IndBox_a_p(:)) % vy=vy
  Atoms(IndBox_a_p(:)) % vz=vz
  
  out=.TRUE.
CONTAINS
  SUBROUTINE Velocities
    REAL(8) :: lm(3),omega(3)
    REAL(8) :: lm_i(3),cm_i(3),tmass_i,iner_i(3,3),am_i(3),ek_i
    REAL(8) :: cm(3),iner(3,3),i_iner(3,3),xc,yc,zc&
         &,am(3),massa,vax,vay,vaz,ek,scale
    INTEGER :: n
    
    lm=0.0_8
    cm=0.0_8
    
    DO n=1,SIZE(mass)
       lm(1)=lm(1)+mass(n)*vx(n)
       lm(2)=lm(2)+mass(n)*vy(n)
       lm(3)=lm(3)+mass(n)*vz(n)
       xc=co(1,1)*xa(n)+co(1,2)*ya(n)+co(1,3)*za(n)
       yc=              co(2,2)*ya(n)+co(2,3)*za(n)
       zc=                            co(3,3)*za(n)
       
       cm(1)=cm(1)+mass(n)*xc
       cm(2)=cm(2)+mass(n)*yc
       cm(3)=cm(3)+mass(n)*zc
       tmass=tmass+mass(n)
    END DO
#ifdef HAVE_MPI
    IF(PI_Nprocs > 1) THEN
       tmass_i=tmass; lm_i=lm;cm_i=cm
       CALL MPI_ALLREDUCE(tmass_i,tmass,1,MPI_REAL8,MPI_SUM&
            &,PI_Comm_Cart,ierr)
       CALL MPI_ALLREDUCE(lm_i,lm,3,MPI_REAL8,MPI_SUM&
            &,PI_Comm_Cart,ierr)
       CALL MPI_ALLREDUCE(cm_i,cm,3,MPI_REAL8,MPI_SUM&
            &,PI_Comm_Cart,ierr)
    END IF
#endif
    cm=cm/tmass
    iner=0.0_8
    am=0.0_8
    DO n=1,SIZE(Mass)
       xc=co(1,1)*xa(n)+co(1,2)*ya(n)+co(1,3)*za(n)
       yc=              co(2,2)*ya(n)+co(2,3)*za(n)
       zc=                            co(3,3)*za(n)
       massa=mass(n)
       
       am(1)=am(1)+massa*(yc*vz(n)-zc*vy(n))
       am(2)=am(2)+massa*(zc*vx(n)-xc*vz(n))
       am(3)=am(3)+massa*(xc*vy(n)-yc*vx(n))
       iner(1,1)=iner(1,1)+massa*(yc**2+zc**2)
       iner(2,2)=iner(2,2)+massa*(zc**2+xc**2)
       iner(3,3)=iner(3,3)+massa*(xc**2+yc**2)
       iner(1,2)=iner(1,2)-massa*xc*yc
       iner(1,3)=iner(1,3)-massa*xc*zc
       iner(2,3)=iner(2,3)-massa*yc*zc
    END DO
    iner(2,1)=iner(1,2)
    iner(3,1)=iner(1,3)
    iner(3,2)=iner(2,3)
    
#ifdef HAVE_MPI
    IF(PI_Nprocs > 1) THEN
       am_i=am;iner_i=iner
       CALL MPI_ALLREDUCE(am_i,am,3,MPI_REAL8,MPI_SUM&
            &,PI_Comm_Cart,ierr)
       CALL MPI_ALLREDUCE(iner_i,iner,9,MPI_REAL8,MPI_SUM&
            &,PI_Comm_Cart,ierr)
    END IF
#endif
    am(1)=am(1)-(cm(2)*lm(3)-cm(3)*lm(2))
    am(2)=am(2)-(cm(3)*lm(1)-cm(1)*lm(3))
    am(3)=am(3)-(cm(1)*lm(2)-cm(2)*lm(1))
    
    IF(.NOT. Matinv_(iner,i_iner)) CALL Print_Errors()
    
    omega(1)=i_iner(1,1)*am(1)+i_iner(1,2)*am(2)+i_iner(1,3)*am(3)
    omega(2)=i_iner(2,1)*am(1)+i_iner(2,2)*am(2)+i_iner(2,3)*am(3)
    omega(3)=i_iner(3,1)*am(1)+i_iner(3,2)*am(2)+i_iner(3,3)*am(3)
    
    DO n=1,SIZE(Mass)
       xc=co(1,1)*xa(n)+co(1,2)*ya(n)+co(1,3)*za(n)
       yc=              co(2,2)*ya(n)+co(2,3)*za(n)
       zc=                            co(3,3)*za(n)
       vax=omega(2)*zc-omega(3)*yc
       vay=omega(3)*xc-omega(1)*zc
       vaz=omega(1)*yc-omega(2)*xc
       vx(n)=vx(n)-vax
       vy(n)=vy(n)-vay
       vz(n)=vz(n)-vaz
    END DO
    lm=0.0_8
    DO n=1,SIZE(mass)
       lm(1)=lm(1)+mass(n)*vx(n)
       lm(2)=lm(2)+mass(n)*vy(n)
       lm(3)=lm(3)+mass(n)*vz(n)
    END DO
#ifdef HAVE_MPI
    IF(PI_Nprocs > 1) THEN
       lm_i=lm
       CALL MPI_ALLREDUCE(lm_i,lm,3,MPI_REAL8,MPI_SUM&
            &,PI_Comm_Cart,ierr)
    END IF
#endif
    lm=lm/tmass
    vx=vx-lm(1)
    vy=vy-lm(2)
    vz=vz-lm(3)
    
    ek=SUM(mass*(vx**2+vy**2+vz**2))
    ek=0.5*ek*efact
#ifdef HAVE_MPI
    IF(PI_Nprocs > 1) THEN
       ek_i=ek
       CALL MPI_ALLREDUCE(ek_i,ek,1,MPI_REAL8,MPI_SUM&
            &,PI_Comm_Cart,ierr)
    END IF
#endif
    IF(ek .EQ. 0.0D0) THEN
       scale=0.0D0
    ELSE
       scale=SQRT(ekt*DBLE(3*SIZE(Atoms))/ek)
    END IF
    vx=vx*scale
    vy=vy*scale
    vz=vz*scale
    
  END SUBROUTINE Velocities
END FUNCTION Atom__vInit_
