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
SUBROUTINE Forces
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Mon Sep 29 2008 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This routine is part of the program  oracDD ----*
!!$---- This routine is contained in the Module Direct ----*

  INTEGER :: ii,n,p,q,r,o,iv,jv,kv,nx,ny,nz,numcell,count_g&
         &,count_gs,ngrp_j,ngrp_js,count0
  REAL(8) :: rcut1,rcut2,rcuts1,rcuts2
  REAL(8) :: xg,yg,zg,xd1,yd1,zd1,xd,yd,zd,xc,yc,zc,rsq,rsqi,r6,r12&
       &,furpar,qforce,conf,uconfa,st1,st2,st3,st4,st5,st6,st7,st8,st9&
       &,emvir,qfx,qfy,qfz,uconf(3),ucoul(3),xpi,ypi,zpi,rsp,xpgi&
       &,ypgi,zpgi,xgg,ygg,zgg,drj,massi,massj,xpgj,ypgj,zpgj,xa,ya&
       &,za,X_PBC,Y_PBC,Z_PBC,chrgei,ucoula,ssvir,rspi,rspqi,alphar&
       &,qt,erfcst,aux1,expcst,ucoul_o(3),uconf_o(3)
  
  INTEGER :: AtSt,AtEn,AtSt_i,AtEn_i,AtSt_j,AtEn_j,ig,j,k,l,i1,jj,j1&
       &,Slv_i,Slv_j,Slv_ij,Id_i,Id_j,Id_ij,nbti,p_mapa,p_j,lij
  

  rcut1=rcut_i(i_p)
  rcut2=rcut1**2
  rcuts1=rcut_o(i_p)
  rcuts2=rcuts1**2

  ucoul=0.0D0
  uconf=0.0D0
  DO ii=1,SIZE(IndBox_g_p)
     ig=IndBox_g_p(ii)
     n=IndBox_g_t(ig)
     xpgi=xpg(ig)
     ypgi=ypg(ig)
     zpgi=zpg(ig)
     p=Chain_xyz(ig) % i
     q=Chain_xyz(ig) % j
     r=Chain_xyz(ig) % k
     AtSt_i=grppt(1,ig)
     AtEn_i=grppt(2,ig)
     count_g=0
     count_gs=0
     DO o=1,SIZE(Ind_xyz)
        iv=Ind_xyz(o) % i
        jv=Ind_xyz(o) % j
        kv=Ind_xyz(o) % k
        nx=mod(mod(p+iv,ncx)+ncx,ncx)
        ny=mod(mod(q+jv,ncy)+ncy,ncy)
        nz=mod(mod(r+kv,ncz)+ncz,ncz)
        numcell=nz+ncz*(ny+ncy*nx)+1
        IF(numcell > ncx*ncy*ncz) STOP
        l=Head_xyz(numcell)
        DO WHILE(l > 0)
           xpgj=xpg(l)
           ypgj=ypg(l)
           zpgj=zpg(l)
           xa=xpgj-xpgi
           ya=ypgj-ypgi
           za=zpgj-zpgi
           X_PBC=PBC(xa)
           Y_PBC=PBC(ya)
           Z_PBC=PBC(za)
           xa=xa+X_PBC
           ya=ya+Y_PBC
           za=za+Z_PBC
           xc=co(1,1)*xa+co(1,2)*ya+co(1,3)*za
           yc=           co(2,2)*ya+co(2,3)*za
           zc=                      co(3,3)*za
           rsq=xc*xc+yc*yc+zc*zc
           IF(rsq <= rcut2) THEN
              count_g=count_g+1
              indGrp(count_g)=l
              Xg_PBC(count_g)=X_PBC
              Yg_PBC(count_g)=Y_PBC
              Zg_PBC(count_g)=Z_PBC
           ELSE IF(rsq <= rcuts2) THEN
              count_gs=count_gs+1
              indGrps(count_gs)=l
              Xgs_PBC(count_gs)=X_PBC
              Ygs_PBC(count_gs)=Y_PBC
              Zgs_PBC(count_gs)=Z_PBC
           END IF
           l=Chain_xyz(l) % p
        END DO
     END DO
     ngrp_j=count_g
     ngrp_js=count_gs

     DO i1=AtSt_i,AtEn_i
        xpi=xp0(i1)
        ypi=yp0(i1)
        zpi=zp0(i1)
        nbti=Id(i1)
        chrgei=chg(i1)
        maplg=.TRUE.
        maplg(i1)=.FALSE.
        IF(ALLOCATED(Maps(i1) % ex)) maplg(Maps(i1) % ex(:))=.FALSE.
        p_mapa=0        
        DO jj=1,ngrp_j
           j1=indGrp(jj)
           AtSt_j=grppt(1,j1)
           AtEn_j=grppt(2,j1)
           DO j=AtSt_j,AtEn_j
              IF(maplg(j)) THEN
                 p_mapa=p_mapa+1 
                 p_index_jj(p_mapa)=jj
                 p_index_j(p_mapa)=j
              END IF
           END DO
        END DO
        Slv_i=Slv(i1)
        Id_i=Id(i1)
        DO p_j=1,p_mapa
           jj=p_index_jj(p_j)
           j=p_index_j(p_j)
           j1=indGrp(jj)

           Slv_j=Slv(j)
           Slv_ij=Slv_j+Slv_i-1

           Id_j=Id(j)
           

           xd=Xg_PBC(jj)
           yd=Yg_PBC(jj)
           zd=Zg_PBC(jj)
           xd1=xpi+xd
           yd1=ypi+yd
           zd1=zpi+zd
           uconfa=0.0D0
           ucoula=0.0D0
!!$
!!$--- Possible Include
!!$           
           lij=Id_j*(Id_j-1)/2+Id_i
           xg=xd1-xp0(j)
           yg=yd1-yp0(j)
           zg=zd1-zp0(j)
           xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
           yc=           co(2,2)*yg+co(2,3)*zg
           zc=                      co(3,3)*zg
           rsq=xc*xc+yc*yc+zc*zc
           rsqi=1.0d0/rsq
           r6=rsqi*rsqi*rsqi
           r12=r6*r6
           ssvir=12.0d0*ecc12(lij)*r12-6.0d0*ecc6(lij)*r6
           qforce=ssvir*rsqi
           emvir = ssvir

           conf=ecc12(lij)*r12-ecc6(lij)*r6
           uconfa=uconfa+conf
                           

           rsp=DSQRT(rsq)
           rspi=1.0D0/rsp
           rspqi=rsqi*rspi
           alphar=alphal*rsp
           
           qt=1.0d0/(1.0d0+qp*alphar)
           expcst=EXP(-alphar*alphar)
           erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*expcst
           furpar=chrgei*chg(j)
           ucoula=ucoula+furpar*erfcst*rspi
           aux1  = furpar*(erfcst+twrtpi*alphar*expcst)*rspqi
           qforce=qforce+aux1
           emvir = emvir*rsqi + aux1
           

           fppx(i1)=fppx(i1)+qforce*xc
           fppy(i1)=fppy(i1)+qforce*yc
           fppz(i1)=fppz(i1)+qforce*zc
           fppx(j)=fppx(j)-qforce*xc
           fppy(j)=fppy(j)-qforce*yc
           fppz(j)=fppz(j)-qforce*zc

           qfx=emvir*xc
           qfy=emvir*yc
           qfz=emvir*zc
           st1 = st1+qfx*xg
           st2 = st2+qfx*yg
           st3 = st3+qfx*zg
           st4 = st4+qfy*xg
           st5 = st5+qfy*yg
           st6 = st6+qfy*zg
           st7 = st7+qfz*xg
           st8 = st8+qfz*yg
           st9 = st9+qfz*zg
!!$
!!$--- Possible Include
!!$
           ucoul(Slv_ij)=ucoul(Slv_ij)+ucoula
           uconf(Slv_ij)=uconf(Slv_ij)+uconfa
        END DO
        

!!$        DO l=1,ngrp_js
!!$           j1=IndGrps(l)
!!$           AtSt_j=grppt(1,j1)
!!$           AtEn_j=grppt(2,j1)
!!$           DO j=AtSt_j,AtEn_j
!!$              
!!$
!!$           END DO
!!$        END DO
     END DO
  END DO
  CALL MPI_ALLREDUCE(ucoul,ucoul_o,3,MPI_REAL8,MPI_SUM,PI_Comm_Cart,ierr)
  CALL MPI_ALLREDUCE(uconf,uconf_o,3,MPI_REAL8,MPI_SUM,PI_Comm_Cart,ierr)

  IF(PI_Node_Cart == 0) THEN
     WRITE(*,*) ucoul_o
     WRITE(*,*) uconf_o
  END IF
END SUBROUTINE Forces
FUNCTION PBC(x) RESULT(out)
  REAL(8) :: out
  REAL(8) :: x
  out=-2.0D0*ANINT(0.5D0*x)
END FUNCTION PBC
