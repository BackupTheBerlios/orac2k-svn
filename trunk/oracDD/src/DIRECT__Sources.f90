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
         &,count_gs,ngrp_j,ngrp_js,count0,i
  REAL(8) :: rcut1,rcut2,rcuts1,rcuts2
  REAL(8) :: xg,yg,zg,xd1,yd1,zd1,xd,yd,zd,xc,yc,zc,rsq,rsqi,r6,r12&
       &,furpar,qforce,conf,uconfa,st1,st2,st3,st4,st5,st6,st7,st8,st9&
       &,emvir,qfx,qfy,qfz,uconf(3),ucoul(3),xpi,ypi,zpi,rsp,xpgi&
       &,ypgi,zpgi,xgg,ygg,zgg,drj,massi,massj,xpgj,ypgj,zpgj,xa,ya&
       &,za,X_PBC,Y_PBC,Z_PBC,chrgei,ucoula,ssvir,rspi,rspqi,alphar&
       &,qt,erfcst,aux1,expcst,ucoul_o(3),uconf_o(3),rcutb,rcutb2
  REAL(8) ::  r2neigh,r2inn,r2out,rinn0,r2inn0,rout0,r2out0,arsout1&
       &,arsout2,arsinn1,arsinn2,rinn,rout,rtolout,rtolinn,auxa,auxb

  INTEGER :: AtSt,AtEn,AtSt_i,AtEn_i,AtSt_j,AtEn_j,ig,j,k,l,i1,jj,j1&
       &,Slv_i,Slv_j,Slv_ij,Id_i,Id_j,nbti,p_mapa,p_j,lij,i_pb&
       &,count_b,count_c,p_mapb
  TYPE(Neigha__), DIMENSION(:), POINTER :: Neigha,Neighb

  IF(i_p > SIZE(Radii)) RETURN

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
  WRITE(*,*) rinn,rinn0
  WRITE(*,*) rout,rout0
!!$  DO ii=1,natom
!!$     Id_i=Id(ii)
!!$     
!!$     lij=Id_ij(Id_i,Id_i)
!!$     WRITE(61,'(2i7,2x,e16.7,f14.6,2e14.6)') ii,Id_i,xp0(ii),chg(ii)&
!!$          &*SQRT(unitc)
!!$  END DO
!!$  
!!$  DO i=22325,22327
!!$     Id_i=Id(i)
!!$     DO j=i,22327
!!$        Id_j=Id(j)
!!$        lij=Id_ij(Id_i,Id_j)
!!$        WRITE(81,'(5i7,2x,2e14.6)') i,j,id_i,id_j,lij,ecc6(lij),ecc12(lij)
!!$     END DO
!!$  END DO

  maplg=.TRUE.
  count_c=0
  DO ii=1,SIZE(IndBox_g_p)
     ig=IndBox_g_p(ii)
     xpgi=xpg(ig)
     ypgi=ypg(ig)
     zpgi=zpg(ig)
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
        xpgj=xpg(l)
        ypgj=ypg(l)
        zpgj=zpg(l)
        xa=xpgi-xpgj
        ya=ypgi-ypgj
        za=zpgi-zpgj
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
        IF(rsq >= r2inn0 .AND. rsq <= r2out) THEN
           count_g=count_g+1
           indGrp(count_g)=l
           Xg_PBC(count_g)=X_PBC
           Yg_PBC(count_g)=Y_PBC
           Zg_PBC(count_g)=Z_PBC
        ELSE IF(rsq > r2out .AND. rsq <= r2out0) THEN
           count_gs=count_gs+1
           indGrps(count_gs)=l
           Xgs_PBC(count_gs)=X_PBC
           Ygs_PBC(count_gs)=Y_PBC
           Zgs_PBC(count_gs)=Z_PBC
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
           Xgs_PBC(count_gs)=X_PBC
           Ygs_PBC(count_gs)=Y_PBC
           Zgs_PBC(count_gs)=Z_PBC
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
        IF(i_pb /= 0 .AND. rsq <= rcutb2) THEN
           count_b=count_b+1
           nei(count_b)=l
        END IF
     END DO
     IF(count_b /= 0) THEN
        Neighb(ig) % no=count_b
        ALLOCATE(Neighb(ig) % nb(count_b))
        Neighb(ig) % nb=nei(1:count_b)
     END IF

     ngrp_j=count_g
     ngrp_js=count_gs
     cmap2(1:ngrp_js)=0.d0
     count_c=count_c+count_gs
     DO i1=AtSt_i,AtEn_i
        xpi=xp0(i1)
        ypi=yp0(i1)
        zpi=zp0(i1)
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
           lij=Id_ij(Id_i,Id_j)
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
                           
           furpar=chrgei*chg(j)
           rsp=DSQRT(rsq)
           rspi=1.0D0/rsp
           rspqi=rsqi*rspi
           alphar=alphal*rsp
              
           qt=1.0d0/(1.0d0+qp*alphar)
           expcst=EXP(-alphar*alphar)
           erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*expcst
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
        DO p_j=1,p_mapa
           jj=p_index_jj(p_j)
           j=p_index_j(p_j)
           j1=indGrps(jj)

           Slv_j=Slv(j)
           Slv_ij=Slv_j+Slv_i-1

           Id_j=Id(j)
           

           xd=Xgs_PBC(jj)
           yd=Ygs_PBC(jj)
           zd=Zgs_PBC(jj)
           xd1=xpi+xd
           yd1=ypi+yd
           zd1=zpi+zd
           uconfa=0.0D0
           ucoula=0.0D0
!!$
!!$--- Possible Include
!!$           
           lij=Id_ij(Id_i,Id_j)
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
           qforce=ssvir*rsqi*swrs(jj)
           conf=ecc12(lij)*r12-ecc6(lij)*r6
           cmap2(jj)=cmap2(jj)+conf
           uconf(Slv_ij)=uconf(Slv_ij)+swrs(jj)*conf
           
           rsp=DSQRT(rsq)
           rspqi=rsqi/rsp
           alphar=alphal*rsp
           qt=1.0d0/(1.0d0+qp*alphar)
           expcst=EXP(-alphar*alphar)
           erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*expcst
           furpar=chrgei*chg(j)
           ucoul(Slv_ij)=ucoul(Slv_ij)+swrs(jj)*furpar*erfcst/rsp
           cmap2(jj)=cmap2(jj)+furpar*erfcst/rsp
           aux1=furpar*(erfcst+twrtpi*alphar*expcst)*rspqi*swrs(jj)
           qforce=qforce+aux1
           emvir=qforce

           fppx(i1)=fppx(i1)+qforce*xc
           fppy(i1)=fppy(i1)+qforce*yc
           fppz(i1)=fppz(i1)+qforce*zc
           fppx(j)=fppx(j)-qforce*xc
           fppy(j)=fppy(j)-qforce*yc
           fppz(j)=fppz(j)-qforce*zc
           st1 = st1+emvir*xc*xg
           st2 = st2+emvir*xc*yg
           st3 = st3+emvir*xc*zg
           st4 = st4+emvir*yc*xg
           st5 = st5+emvir*yc*yg
           st6 = st6+emvir*yc*zg
           st7 = st7+emvir*zc*xg
           st8 = st8+emvir*zc*yg
           st9 = st9+emvir*zc*zg
        END DO
        maplg(i1)=.TRUE.
        IF(ALLOCATED(Maps(i1) % ex)) maplg(Maps(i1) % ex(:))=.TRUE.
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
           xpi=xp0(i1)
           ypi=yp0(i1)
           zpi=zp0(i1)
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
           xg=-Xgs_PBC(jj)
           yg=-Ygs_PBC(jj)
           zg=-Zgs_PBC(jj)
           DO j1=grppt(1,j),grppt(2,j)
              massj=gmass(j1)
              fppx(j1)=fppx(j1)-massj*xmap3(jj)
              fppy(j1)=fppy(j1)-massj*ymap3(jj)
              fppz(j1)=fppz(j1)-massj*zmap3(jj)
              xpi=-(xp0(j1)+xg)*massj
              ypi=-(yp0(j1)+yg)*massj
              zpi=-(zp0(j1)+zg)*massj
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
