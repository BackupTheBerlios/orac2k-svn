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
#include "forces.h"
#include "Erfc_Spline.h"
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
         &,count_gs,ngrp_j,ngrp_js,count0,i
  REAL(8) :: rcut1,rcut2,rcuts1,rcuts2
  Real(8), Save :: xa,xd1,xd,xc,ya,yd1,yd,yc,za,zd1,zd,zc
  REAL(8) :: rsq,rsqi,r6,r12&
       &,furpar,qforce,conf,uconfa,st1,st2,st3,st4,st5,st6,st7,st8,st9&
       &,emvir,qfx,qfy,qfz,uconf(3),ucoul(3),xpi,ypi,zpi,rsp,xpgi&
       &,ypgi,zpgi,xgg,ygg,zgg,drj,massi,massj,xpgci,ypgci,zpgci&
       &,X_PBC,Y_PBC,Z_PBC,chrgei,ucoula,ssvir,rspi,rspqi,alphar& 
       &,qt,erfcst,aux1,expcst,ucoul_o(3),uconf_o(3),rcutb,rcutb2
  REAL(8) ::  r2neigh,r2inn,r2out,rinn0,r2inn0,rout0,r2out0,arsout1&
       &,arsout2,arsinn1,arsinn2,rinn,rout,rtolout,rtolinn,auxa,auxb&
       &,h,MyErfc_,MyDErfc_,fpix,fpiy,fpiz

  INTEGER :: AtSt,AtEn,AtSt_i,AtEn_i,AtSt_j,AtEn_j,ig,j,k,l,i1,jj,j1&
       &,Slv_i,Slv_j,Slv_ij,Id_i,Id_j,nbti,p_mapa,p_j,lij,i_pb&
       &,count_b,count_c,p_mapb,nol,natom1,i_c
  REAL(8), SAVE :: Tol_q=1.0D-5

  LOGICAL :: lskip_ewald,ok_LJ,ok_Co
  LOGICAL, ALLOCATABLE :: ol(:)
  TYPE(Neigha__), DIMENSION(:), POINTER :: Neigha,Neighb
  Real(8), Allocatable :: xpc(:),ypc(:),zpc(:),xpgc(:),ypgc(:),zpgc(:)

  IF(i_p > SIZE(Radii)) RETURN

  Allocate(xpc(natom),ypc(natom),zpc(natom))
  Allocate(xpgc(ngroup),ypgc(ngroup),zpgc(ngroup))
  xpc=co(1,1)*xp0+co(1,2)*yp0+co(1,3)*zp0
  ypc=            co(2,2)*yp0+co(2,3)*zp0
  zpc=                        co(3,3)*zp0

  xpgc=co(1,1)*xpg+co(1,2)*ypg+co(1,3)*zpg
  ypgc=            co(2,2)*ypg+co(2,3)*zpg
  zpgc=                        co(3,3)*zpg
  

  Neigha=>List(i_p) % Neigh
  i_pb=i_p-1
  rcutb=0.0D0
  rcutb2=0.0D0
!!$          Translation scheme orac2k => oracDD
!!$
!!$
!!$         rinn    rinn0      rout     rout0
!!$       rcut_i(1) rcut_o(1) rcut_i(2) rcut_o(2)  
!!$            |      |         |        |
!!$            |      |         |        |
!!$            v      v         v        v
!!$------------xxxxxxxx----------xxxxxxxxx----
!!$           A       A          A       A
!!$           |       |          |       |
!!$           |       |          |       |
!!$           rcut_s(2)          rcut_s(2) 
!!$           rtolin             rtolout
!!$

  rcut1=Radii(i_p) % inn 
  rcut2=rcut1**2
  rcuts1=Radii(i_p) % out
  rcuts2=rcuts1**2

  rinn=0.0D0
  rinn0=0.0D0
  rtolinn=0.0D0

  rout=Radii(i_p) % inn
  rout0=Radii(i_p) % out
  rtolout=Radii(i_p) % shell

  IF(i_pb /= 0) THEN
     Neighb=>List(i_pb) % Neigh
     rcutb=Radii(i_pb) % out + Radii(i_pb) % update 
     rcutb2=rcutb**2

     rinn=Radii(i_pb) % inn
     rinn0=Radii(i_pb) % out
     rtolinn=Radii(i_pb) % shell
  END IF

  r2out=rout**2
  r2out0=rout0**2
  r2inn=rinn**2
  r2inn0=rinn0**2

  arsout1=rout0-3.0d0*rout
  arsout2=rtolout**3
  arsinn1=rinn0-3.0d0*rinn
  arsinn2=rtolinn**3
  ucoul=0.0D0
  uconf=0.0D0

  maplg=.TRUE.
  count_c=0
  st1=0.0D0
  st2=0.0D0
  st3=0.0D0
  st4=0.0D0
  st5=0.0D0
  st6=0.0D0
  st7=0.0D0
  st8=0.0D0
  st9=0.0D0
  qt=1.0d0/(1.0d0+qp*alphal*rinn)
  expcst=DEXP(-alphal*rinn*alphal*rinn)
  erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*expcst
  lskip_ewald = erfcst.lt.1.0d-4
  IF(Ewald__Param % NoSkip) lskip_ewald = .FALSE.
  fppx=0.0_8
  fppy=0.0_8
  fppz=0.0_8
  DO ig=1,ngroup
     xpgi=xpg(ig)
     ypgi=ypg(ig)
     zpgi=zpg(ig)

     xpgci=xpgc(ig)
     ypgci=ypgc(ig)
     zpgci=zpgc(ig)

     AtSt_i=grppt(1,ig)
     AtEn_i=grppt(2,ig)
     count_g=0
     count_gs=0
     IF(i_pb /= 0) THEN
        IF(ALLOCATED(Neighb(ig) % nb)) DEALLOCATE(Neighb(ig) % nb)
        Neighb(ig) % no = 0 
     END IF
     count_b=0
     DO o=1,Neigha(ig) % no
        l=Neigha(ig) % nb(o)

        xa=xpgi-xpg(l)
        xa=_PBC(xa)

        ya=ypgi-ypg(l)
        ya=_PBC(ya)

        za=zpgi-zpg(l)
        za=_PBC(za)

        xd=co(1,1)*xa+co(1,2)*ya+co(1,3)*za
        xc=xpgci-xpgc(l)+xd
        rsq=xc*xc

        yd=           co(2,2)*ya+co(2,3)*za
        yc=ypgci-ypgc(l)+yd
        rsq=rsq+yc*yc

        zd=                      co(3,3)*za
        zc=zpgci-zpgc(l)+zd
        rsq=rsq+zc*zc
        
        IF(rsq >= r2inn0 .AND. rsq <= r2out) THEN
           count_g=count_g+1
           indGrp(count_g)=l
           Xg_PBC(count_g)=xd
           Yg_PBC(count_g)=yd
           Zg_PBC(count_g)=zd
        ELSE IF(rsq > r2out .AND. rsq <= r2out0) THEN
           count_gs=count_gs+1
           indGrps(count_gs)=l
           Xgs_PBC(count_gs)=xd
           Ygs_PBC(count_gs)=yd
           Zgs_PBC(count_gs)=zd

           xcs(count_gs)=xc
           ycs(count_gs)=yc
           zcs(count_gs)=zc
           rsp=DSQRT(rsq)
           auxa=(arsout1+2.0d0*rsp)/arsout2
           auxb=rout0-rsp
           swrs(count_gs)=auxa*auxb**2
           dswrs(count_gs)=-2.0d0*auxa*auxb+2.0d0*auxb**2/arsout2
           dswrs(count_gs)=dswrs(count_gs)/rsp
        ELSE IF(rsq > r2inn .AND. rsq <= r2inn0) THEN
           count_gs=count_gs+1
           indGrps(count_gs)=l
           Xgs_PBC(count_gs)=xd
           Ygs_PBC(count_gs)=yd
           Zgs_PBC(count_gs)=zd
           xcs(count_gs)=xc
           ycs(count_gs)=yc
           zcs(count_gs)=zc
           rsp=DSQRT(rsq)
           auxa=(arsinn1+2.0d0*rsp)/arsinn2
           auxb=rinn0-rsp
           swrs(count_gs)=1.d0-auxa*auxb**2
           dswrs(count_gs)=2.0d0*auxa*auxb-2.0d0*auxb**2/arsinn2
           dswrs(count_gs)=dswrs(count_gs)/rsp
        END IF
        IF(i_pb /= 0) THEN
           IF(rsq <= rcutb2) THEN
              count_b=count_b+1
              neib(count_b)=l
           END IF
        END IF
     END DO
     IF(count_b /= 0) THEN
        Neighb(ig) % no=count_b
        ALLOCATE(Neighb(ig) % nb(count_b))
        Neighb(ig) % nb=neib(1:count_b)
     END IF

     ngrp_j=count_g
     ngrp_js=count_gs
     cmap2(1:ngrp_js)=0.d0
     count_c=count_c+count_gs+count_g
     DO i1=AtSt_i,AtEn_i
        xpi=xpc(i1)
        ypi=ypc(i1)
        zpi=zpc(i1)
        nbti=Id(i1)
        chrgei=chg(i1)
        Slv_i=Slv(i1)
        Id_i=Id(i1)
        maplg(i1)=.FALSE.
        IF(ALLOCATED(Maps(i1) % ex)) maplg(Maps(i1) % ex(:))=.FALSE.
        p_mapa=0        
        p_mapb=0        
        DO jj=1,ngrp_j
           j1=indGrp(jj)
           AtSt_j=grppt(1,j1)
           AtEn_j=grppt(2,j1)
           IF(ig == j1) AtSt_j=i1+1
           DO j=AtSt_j,AtEn_j
              IF(maplg(j)) THEN
                 p_mapa=p_mapa+1 
                 p_index_jj(p_mapa)=jj
                 p_index_j(p_mapa)=j
              END IF
           END DO
        END DO
        fpix=0.0_8; fpiy=0.0_8; fpiz=0.0_8
        IF(Erfc_Switch .AND. (.NOT. lskip_ewald)) THEN
           DO p_j=1,p_mapa
              jj=p_index_jj(p_j)
              j=p_index_j(p_j)
              j1=indGrp(jj)
              
              Slv_j=Slv(j)
              Slv_ij=Slv_j+Slv_i-1
              
              Id_j=Id(j)
              
              uconfa=0.0D0
              ucoula=0.0D0
!!$
!!$--- Possible Include
!!$           
              lij=Id_ij(Id_i,Id_j)

              xd=Xg_PBC(jj)
              xd1=xpi+xd
              xc=xd1-xpc(j)
              rsq=xc*xc

              yd=Yg_PBC(jj)
              yd1=ypi+yd
              yc=yd1-ypc(j)
              rsq=rsq+yc*yc

              zd=Zg_PBC(jj)
              zd1=zpi+zd
              zc=zd1-zpc(j)
              rsq=rsq+zc*zc

              qforce=0.0D0; emvir=0.0D0; aux1=0.0D0
              furpar=chrgei*chg(j)

              ok_LJ=rsq < eccc(lij)
              ok_Co=(rsq < Ewald_cut .AND. ABS(furpar) > Tol_q) .OR. Ewald__Param % noskip 
              

              IF((.NOT. ok_LJ) .AND. (.NOT. ok_Co)) CYCLE

              IF(ok_LJ) THEN
                 rsqi=1.0d0/rsq
                 r6=rsqi*rsqi*rsqi
                 r12=r6*r6
                 ssvir=12.0d0*ecc12(lij)*r12-6.0d0*ecc6(lij)*r6
                 qforce=ssvir*rsqi
                 emvir = ssvir*rsqi                 
                 conf=ecc12(lij)*r12-ecc6(lij)*r6
                 uconfa=uconfa+conf
              END IF
              
              IF(ok_Co) THEN
                 rsp=SQRT(rsq)
                 rspi=1.0_8/rsp
                 i_c=INT((rsp-E_xbeg)/E_dx)+1
                 h=rsp-DBLE(i_c-1)*E_dx-E_xbeg
                 MyErfc_=_Erfc(i_c,h)
                 MyDErfc_=_DErfc(i_c,h)                 
                 ucoula=ucoula+furpar*MyErfc_
                 aux1  = -rspi*furpar*MyDErfc_
                 qforce=qforce+aux1
              END IF

              emvir = emvir + aux1
              aux1=qforce*xc
              fpix=fpix+aux1
              fppx(j)=fppx(j)-aux1
              
              aux1=qforce*yc
              fpiy=fpiy+aux1
              fppy(j)=fppy(j)-aux1
              
              aux1=qforce*zc
              fpiz=fpiz+aux1
              fppz(j)=fppz(j)-aux1
              
              qfx=emvir*xc
              st1 = st1+qfx*xc
              st2 = st2+qfx*yc
              st3 = st3+qfx*zc
              qfy=emvir*yc
              st4 = st4+qfy*xc
              st5 = st5+qfy*yc
              st6 = st6+qfy*zc
              qfz=emvir*zc
              st7 = st7+qfz*xc
              st8 = st8+qfz*yc
              st9 = st9+qfz*zc
!!$
!!$--- Possible Include
!!$
              ucoul(Slv_ij)=ucoul(Slv_ij)+ucoula
              uconf(Slv_ij)=uconf(Slv_ij)+uconfa
           END DO
        ELSE IF((.NOT. Erfc_switch) .AND. (.NOT. lskip_ewald)) THEN
           DO p_j=1,p_mapa
              jj=p_index_jj(p_j)
              j=p_index_j(p_j)
              j1=indGrp(jj)
              
              Slv_j=Slv(j)
              Slv_ij=Slv_j+Slv_i-1
              
              Id_j=Id(j)
              
              uconfa=0.0D0
              ucoula=0.0D0
!!$
!!$--- Possible Include
!!$           
              lij=Id_ij(Id_i,Id_j)

              xd=Xg_PBC(jj)
              xd1=xpi+xd
              xc=xd1-xpc(j)
              rsq=xc*xc

              yd=Yg_PBC(jj)
              yd1=ypi+yd
              yc=yd1-ypc(j)
              rsq=rsq+yc*yc

              zd=Zg_PBC(jj)
              zd1=zpi+zd
              zc=zd1-zpc(j)
              rsq=rsq+zc*zc


              qforce=0.0D0; emvir=0.0D0; aux1=0.0D0
              rsqi=1.0d0/rsq
              furpar=chrgei*chg(j)

              ok_LJ=rsq < eccc(lij)
              ok_Co=(rsq < Ewald_cut .AND. ABS(furpar) > Tol_q) .OR. Ewald__Param % noskip 

              IF((.NOT. ok_LJ) .AND. (.NOT. ok_Co)) CYCLE
              IF(ok_LJ) THEN
                 r6=rsqi*rsqi*rsqi
                 r12=r6*r6
                 ssvir=12.0d0*ecc12(lij)*r12-6.0d0*ecc6(lij)*r6
                 qforce=ssvir*rsqi
                 emvir = ssvir*rsqi
                 
                 conf=ecc12(lij)*r12-ecc6(lij)*r6
                 uconfa=uconfa+conf
              END IF
              
              IF(ok_Co) THEN
                 rsp=SQRT(rsq)
                 rspi=1.0D0/rsp
                 rspqi=rsqi*rspi
                 alphar=alphal*rsp
                 
                 qt=1.0d0/(1.0d0+qp*alphar)
                 expcst=EXP(-alphar*alphar)
                 erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*expcst
                 ucoula=ucoula+furpar*erfcst*rspi                 
                 aux1  = furpar*(erfcst+twrtpi*alphar*expcst)*rspqi
                 qforce=qforce+aux1
              END IF

              emvir = emvir + aux1
              aux1=qforce*xc
              fpix=fpix+aux1
              fppx(j)=fppx(j)-aux1
              
              aux1=qforce*yc
              fpiy=fpiy+aux1
              fppy(j)=fppy(j)-aux1
              
              aux1=qforce*zc
              fpiz=fpiz+aux1
              fppz(j)=fppz(j)-aux1
              
              qfx=emvir*xc
              st1 = st1+qfx*xc
              st2 = st2+qfx*yc
              st3 = st3+qfx*zc
              qfy=emvir*yc
              st4 = st4+qfy*xc
              st5 = st5+qfy*yc
              st6 = st6+qfy*zc
              qfz=emvir*zc
              st7 = st7+qfz*xc
              st8 = st8+qfz*yc
              st9 = st9+qfz*zc
!!$
!!$--- Possible Include
!!$
              ucoul(Slv_ij)=ucoul(Slv_ij)+ucoula
              uconf(Slv_ij)=uconf(Slv_ij)+uconfa
           END DO
        ELSE IF(lskip_ewald) THEN
           DO p_j=1,p_mapa
              jj=p_index_jj(p_j)
              j=p_index_j(p_j)
              j1=indGrp(jj)
              
              
              uconfa=0.0D0
              ucoula=0.0D0
!!$
!!$--- Possible Include
!!$           
              xd=Xg_PBC(jj)
              xd1=xpi+xd
              xc=xd1-xpc(j)
              rsq=xc*xc

              yd=Yg_PBC(jj)
              yd1=ypi+yd
              yc=yd1-ypc(j)
              rsq=rsq+yc*yc

              zd=Zg_PBC(jj)
              zd1=zpi+zd
              zc=zd1-zpc(j)
              rsq=rsq+zc*zc

              Id_j=Id(j)
              lij=Id_ij(Id_i,Id_j)
              IF(rsq > eccc(lij)) CYCLE

              Slv_j=Slv(j)
              Slv_ij=Slv_j+Slv_i-1
              

              rsqi=1.0d0/rsq
              qforce=0.0D0; emvir=0.0D0
              r6=rsqi*rsqi*rsqi
              r12=r6*r6
              ssvir=12.0d0*ecc12(lij)*r12-6.0d0*ecc6(lij)*r6
              qforce=ssvir*rsqi
              emvir = ssvir*rsqi
              
              conf=ecc12(lij)*r12-ecc6(lij)*r6
              uconfa=uconfa+conf

              aux1=qforce*xc
              fpix=fpix+aux1
              fppx(j)=fppx(j)-aux1
              
              aux1=qforce*yc
              fpiy=fpiy+aux1
              fppy(j)=fppy(j)-aux1
              
              aux1=qforce*zc
              fpiz=fpiz+aux1
              fppz(j)=fppz(j)-aux1
              
              qfx=emvir*xc
              st1 = st1+qfx*xc
              st2 = st2+qfx*yc
              st3 = st3+qfx*zc
              qfy=emvir*yc
              st4 = st4+qfy*xc
              st5 = st5+qfy*yc
              st6 = st6+qfy*zc
              qfz=emvir*zc
              st7 = st7+qfz*xc
              st8 = st8+qfz*yc
              st9 = st9+qfz*zc
!!$
!!$--- Possible Include
!!$
              ucoul(Slv_ij)=ucoul(Slv_ij)+ucoula
              uconf(Slv_ij)=uconf(Slv_ij)+uconfa
           END DO

        END IF
        p_mapa=0        
        DO jj=1,ngrp_js
           j1=indGrps(jj)
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
        IF(Erfc_Switch .AND. (.NOT. lskip_ewald)) THEN
           DO p_j=1,p_mapa
              jj=p_index_jj(p_j)
              j=p_index_j(p_j)
              j1=indGrps(jj)
              
              Slv_j=Slv(j)
              Slv_ij=Slv_j+Slv_i-1
              
              Id_j=Id(j)
              
              
              uconfa=0.0D0
              ucoula=0.0D0
!!$
!!$--- Possible Include
!!$           
              lij=Id_ij(Id_i,Id_j)

              xd=Xgs_PBC(jj)
              xd1=xpi+xd
              xc=xd1-xpc(j)
              rsq=xc*xc

              yd=Ygs_PBC(jj)
              yd1=ypi+yd
              yc=yd1-ypc(j)
              rsq=rsq+yc*yc

              zd=Zgs_PBC(jj)
              zd1=zpi+zd
              zc=zd1-zpc(j)
              rsq=rsq+zc*zc

              rsqi=1.0d0/rsq
              
              qforce=0.0D0
              furpar=chrgei*chg(j)

              ok_LJ=rsq < eccc(lij)
              ok_Co=(rsq < Ewald_cut .AND. ABS(furpar) > Tol_q) .OR. Ewald__Param % noskip 
              IF((.NOT. ok_LJ) .AND. (.NOT. ok_Co)) CYCLE

              IF(ok_LJ) THEN
                 r6=rsqi*rsqi*rsqi
                 r12=r6*r6
                 ssvir=12.0d0*ecc12(lij)*r12-6.0d0*ecc6(lij)*r6
                 qforce=ssvir*rsqi*swrs(jj)
                 conf=ecc12(lij)*r12-ecc6(lij)*r6
                 cmap2(jj)=cmap2(jj)+conf
                 uconf(Slv_ij)=uconf(Slv_ij)+swrs(jj)*conf
              END IF
              
              IF(ok_Co) THEN
                 rsp=SQRT(rsq)
                 rspi=1.0_8/rsp
                 i_c=INT((rsp-E_xbeg)/E_dx)+1
                 h=rsp-tau(i_c)
                 
                 MyErfc_=_Erfc(i_c,h)
                 MyDErfc_=_DErfc(i_c,h)
                 
                 ucoul(Slv_ij)=ucoul(Slv_ij)+swrs(jj)*furpar*MyErfc_
                 cmap2(jj)=cmap2(jj)+furpar*MyErfc_
                 aux1  = -swrs(jj)*rspi*furpar*MyDErfc_
                 qforce=qforce+aux1
              END IF
              emvir=qforce

              fpix=fpix+qforce*xc
              fpiy=fpiy+qforce*yc
              fpiz=fpiz+qforce*zc
              fppx(j)=fppx(j)-qforce*xc
              fppy(j)=fppy(j)-qforce*yc
              fppz(j)=fppz(j)-qforce*zc
              st1 = st1+emvir*xc*xc
              st2 = st2+emvir*xc*yc
              st3 = st3+emvir*xc*zc
              st4 = st4+emvir*yc*xc
              st5 = st5+emvir*yc*yc
              st6 = st6+emvir*yc*zc
              st7 = st7+emvir*zc*xc
              st8 = st8+emvir*zc*yc
              st9 = st9+emvir*zc*zc
           END DO
        ELSE IF((.NOT. Erfc_switch) .AND. (.NOT. lskip_ewald)) THEN
           DO p_j=1,p_mapa
              jj=p_index_jj(p_j)
              j=p_index_j(p_j)
              j1=indGrps(jj)
              
              Slv_j=Slv(j)
              Slv_ij=Slv_j+Slv_i-1
              
              Id_j=Id(j)
              
              
              uconfa=0.0D0
              ucoula=0.0D0
!!$
!!$--- Possible Include
!!$           
              lij=Id_ij(Id_i,Id_j)

              xd=Xgs_PBC(jj)
              xd1=xpi+xd
              xc=xd1-xpc(j)
              rsq=xc*xc

              yd=Ygs_PBC(jj)
              yd1=ypi+yd
              yc=yd1-ypc(j)
              rsq=rsq+yc*yc

              zd=Zgs_PBC(jj)
              zd1=zpi+zd
              zc=zd1-zpc(j)
              rsq=rsq+zc*zc

              rsqi=1.0d0/rsq
              
              qforce=0.0D0;furpar=chrgei*chg(j)
              ok_LJ=rsq < eccc(lij)
              ok_Co=(rsq < Ewald_cut .AND. ABS(furpar) > Tol_q) .OR. Ewald__Param % noskip 
              IF((.NOT. ok_LJ) .AND. (.NOT. ok_Co)) CYCLE
              IF(ok_LJ) THEN
                 r6=rsqi*rsqi*rsqi
                 r12=r6*r6
                 ssvir=12.0d0*ecc12(lij)*r12-6.0d0*ecc6(lij)*r6
                 qforce=ssvir*rsqi*swrs(jj)
                 conf=ecc12(lij)*r12-ecc6(lij)*r6
                 cmap2(jj)=cmap2(jj)+conf
                 uconf(Slv_ij)=uconf(Slv_ij)+swrs(jj)*conf
              END IF
              
              IF(ok_Co) THEN
                 rsp=DSQRT(rsq)
                 rspi=1.0_8/rsp
                 rspqi=rsqi/rsp
                 alphar=alphal*rsp
                 qt=1.0d0/(1.0d0+qp*alphar)
                 expcst=EXP(-alphar*alphar)
                 erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*expcst
                 ucoul(Slv_ij)=ucoul(Slv_ij)+swrs(jj)*furpar*erfcst/rsp
                 cmap2(jj)=cmap2(jj)+furpar*erfcst/rsp
                 aux1=furpar*(erfcst+twrtpi*alphar*expcst)*rspqi*swrs(jj)
                 qforce=qforce+aux1
              END IF

              emvir=qforce
              
              fpix=fpix+qforce*xc
              fpiy=fpiy+qforce*yc
              fpiz=fpiz+qforce*zc
              fppx(j)=fppx(j)-qforce*xc
              fppy(j)=fppy(j)-qforce*yc
              fppz(j)=fppz(j)-qforce*zc
              st1 = st1+emvir*xc*xc
              st2 = st2+emvir*xc*yc
              st3 = st3+emvir*xc*zc
              st4 = st4+emvir*yc*xc
              st5 = st5+emvir*yc*yc
              st6 = st6+emvir*yc*zc
              st7 = st7+emvir*zc*xc
              st8 = st8+emvir*zc*yc
              st9 = st9+emvir*zc*zc
           END DO
        ELSE
           DO p_j=1,p_mapa
              jj=p_index_jj(p_j)
              j=p_index_j(p_j)
              j1=indGrps(jj)
              
              Slv_j=Slv(j)
              Slv_ij=Slv_j+Slv_i-1
              
              Id_j=Id(j)
              
              
              uconfa=0.0D0
              ucoula=0.0D0
!!$
!!$--- Possible Include
!!$           
              lij=Id_ij(Id_i,Id_j)

              xd=Xgs_PBC(jj)
              xd1=xpi+xd
              xc=xd1-xpc(j)
              rsq=xc*xc

              yd=Ygs_PBC(jj)
              yd1=ypi+yd
              yc=yd1-ypc(j)
              rsq=rsq+yc*yc

              zd=Zgs_PBC(jj)
              zd1=zpi+zd
              zc=zd1-zpc(j)
              rsq=rsq+zc*zc

              rsqi=1.0d0/rsq
              
              qforce=0.0D0
              IF(rsq > eccc(lij)) CYCLE
              r6=rsqi*rsqi*rsqi
              r12=r6*r6
              ssvir=12.0d0*ecc12(lij)*r12-6.0d0*ecc6(lij)*r6
              qforce=ssvir*rsqi*swrs(jj)
              conf=ecc12(lij)*r12-ecc6(lij)*r6
              cmap2(jj)=cmap2(jj)+conf
              uconf(Slv_ij)=uconf(Slv_ij)+swrs(jj)*conf
              
              emvir=qforce
              
              fpix=fpix+qforce*xc
              fpiy=fpiy+qforce*yc
              fpiz=fpiz+qforce*zc
              fppx(j)=fppx(j)-qforce*xc
              fppy(j)=fppy(j)-qforce*yc
              fppz(j)=fppz(j)-qforce*zc
              st1 = st1+emvir*xc*xc
              st2 = st2+emvir*xc*yc
              st3 = st3+emvir*xc*zc
              st4 = st4+emvir*yc*xc
              st5 = st5+emvir*yc*yc
              st6 = st6+emvir*yc*zc
              st7 = st7+emvir*zc*xc
              st8 = st8+emvir*zc*yc
              st9 = st9+emvir*zc*zc
           END DO
              
        END IF
        maplg(i1)=.TRUE.
        IF(ALLOCATED(Maps(i1) % ex)) maplg(Maps(i1) % ex(:))=.TRUE.
        fppx(i1)=fppx(i1)+fpix
        fppy(i1)=fppy(i1)+fpiy
        fppz(i1)=fppz(i1)+fpiz
     END DO
!!$
!!$===     add the S*V term to the atomic forces
!!$
     DO jj=1,ngrp_js
        xmap3(jj)=-dswrs(jj)*cmap2(jj)*xcs(jj)
        ymap3(jj)=-dswrs(jj)*cmap2(jj)*ycs(jj)
        zmap3(jj)=-dswrs(jj)*cmap2(jj)*zcs(jj)
     END DO
     
     IF(ngrp_js /= 0) THEN
        DO i1=grppt(1,ig),grppt(2,ig)
           massi=gmass(i1)
           xpi=xpc(i1)
           ypi=ypc(i1)
           zpi=zpc(i1)
!DEC$ IVDEP
           DO jj=1,ngrp_js
              fppx(i1)=fppx(i1)+massi*xmap3(jj)
              fppy(i1)=fppy(i1)+massi*ymap3(jj)
              fppz(i1)=fppz(i1)+massi*zmap3(jj)
              st1=st1+massi*xmap3(jj)*xpi
              st2=st2+massi*xmap3(jj)*ypi
              st3=st3+massi*xmap3(jj)*zpi
              st4=st4+massi*ymap3(jj)*xpi
              st5=st5+massi*ymap3(jj)*ypi
              st6=st6+massi*ymap3(jj)*zpi
              st7=st7+massi*zmap3(jj)*xpi
              st8=st8+massi*zmap3(jj)*ypi
              st9=st9+massi*zmap3(jj)*zpi
           END DO
        END DO
        DO jj=1,ngrp_js
           j=Indgrps(jj)
           xc=-Xgs_PBC(jj)
           yc=-Ygs_PBC(jj)
           zc=-Zgs_PBC(jj)
!DEC$ IVDEP
           DO j1=grppt(1,j),grppt(2,j)
              massj=gmass(j1)
              fppx(j1)=fppx(j1)-massj*xmap3(jj)
              fppy(j1)=fppy(j1)-massj*ymap3(jj)
              fppz(j1)=fppz(j1)-massj*zmap3(jj)
              xpi=-(xpc(j1)+xc)*massj
              ypi=-(ypc(j1)+yc)*massj
              zpi=-(zpc(j1)+zc)*massj
              st1=st1+xmap3(jj)*xpi
              st2=st2+xmap3(jj)*ypi
              st3=st3+xmap3(jj)*zpi
              st4=st4+ymap3(jj)*xpi
              st5=st5+ymap3(jj)*ypi
              st6=st6+ymap3(jj)*zpi
              st7=st7+zmap3(jj)*xpi
              st8=st8+zmap3(jj)*ypi
              st9=st9+zmap3(jj)*zpi
           END DO
        END DO
     END IF
  END DO

  fp(IndBox_a_t(:)) % x = fp(IndBox_a_t(:)) % x + fppx(:)
  fp(IndBox_a_t(:)) % y = fp(IndBox_a_t(:)) % y + fppy(:)
  fp(IndBox_a_t(:)) % z = fp(IndBox_a_t(:)) % z + fppz(:)

  CALL EN_LJ_(i_p,uconf(3),uconf(1),uconf(2))
  CALL EN_Coul_Dir_(i_p,ucoul(3),ucoul(1),ucoul(2))

END SUBROUTINE Forces
