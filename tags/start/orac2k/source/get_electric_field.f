      SUBROUTINE get_electric_field(xp0,yp0,zp0,xpa,ypa,zpa,xpga,ypga
     &     ,zpga,xpcma,ypcma,zpcma,mapnl,ingrpp_ef,ingrp_ef,ingrp_ef_x
     &     ,fpx,fpy,fpz,flj_x,flj_y,flj_z,Ex_rec,Ey_rec,Ez_rec,Edx
     &     ,Edy,Edz,Ex_Cor,Ey_Cor,Ez_Cor,worka,cpu,ncpu
     &     ,nstart,nend,nlocal,nstart_a,nend_a,nlocal_a,node,nodex,nodey
     &     ,nodez,ictxt,npy,npz,descQ,nprocs,ncube,tag_bndg,Ugrp,eer
     &     ,uself,uself_dip,Udirect,ene_ferrf,U_Thole,ELJ,ucoul_slt
     &     ,ucoul_slv,ucoul_ss,uconf_slt,uconf_slv,uconf_ss,charges
     &     ,dipole,do_LJ,skip_direct,Gauss_mode,Direct_Only,set_zero)

************************************************************************
*                                                                      *
*              Author:  Matteo Ceccarelli                              *
*              CECAM/ENS Lyon, FRANCE, may 99                          *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER node,nprocs,ncube,mapnl(*),ingrpp_ef,ingrp_ef(2,*)
      INTEGER nodex,nodey,nodez,ictxt,npy,npz,descQ(*),ingrp_ef_x(*)
      INTEGER tag_bndg(*)
      REAL*8  xpa(*),ypa(*),zpa(*),xpcma(*),ypcma(*),zpcma(*)
     &        ,xpga(*),ypga(*),zpga(*),xp0(*),yp0(*),zp0(*)
      REAL*8  Edx(*),Edy(*),Edz(*),dipole(3,*),charges(*)
      REAL*8  Ex_rec(*),Ey_rec(*),Ez_rec(*)
      REAL*8  Ex_Cor(*),Ey_Cor(*),Ez_Cor(*)
      REAL*8  fpx(*),fpy(*),fpz(*),cpu,Ugrp,Udirect,eer,uself,ene_ferrf
     &     ,uself_dip,U_Thole,ELJ,flj_x(*),flj_y(*),flj_z(*),ucoul_slt
     &     ,ucoul_slv,ucoul_ss,uconf_slt,uconf_slv,uconf_ss
      INTEGER nstart,nend,nlocal,nstart_a,nend_a,nlocal_a,worka(*)
      INTEGER npp_m,ncpu
      LOGICAL do_LJ,skip_direct,Gauss_mode,Direct_Only,set_zero


*----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'parst.h'
      INCLUDE 'cpropar.h'
      INCLUDE 'pme.h'

      REAL*8  work(m1),qt(4),fdx(m1),fdy(m1),fdz(m1)

      COMMON /rag1/ work,fdx,fdy,fdz

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8  urcs,urcp,urcsp,fudgec,ene_rec(m1)
      INTEGER nnlppaux(sitslu+1)
      INTEGER i,j,k,npp
      REAL*8  unb14_slt,cnb14_slt,ungrp_slt,cngrp_slt,uptors_slt
     &       ,unb14_slv,cnb14_slv,ungrp_slv,cngrp_slv,uptors_slv
     &     ,self_slt,self_slv,coulgrp_slt,coulgrp_slv
     &     ,confgrp_slt,confgrp_slv,Ucos,Ucop,Ucosp,ene,enedip

*----------------------- EXECUTABLE STATEMENTS ------------------------*
*=======================================================================
c------------ initialize some variables
*=======================================================================

      fudgec = 1.0d0 - fudge

      
      IF(do_LJ) THEN
         ELJ=0.0D0
         uconf_slt=0.0D0
         uconf_slv=0.0D0
         uconf_ss=0.0D0
         CALL zeroa(flj_x,flj_y,flj_z,ntap,1)
      END IF
      IF(set_zero) THEN
         U_Thole=0.0D0
         Udirect   = 0.0d0
         Ugrp = 0.0d0
         eer   = 0.0d0
         ene_ferrf = 0.0d0
         cngrp_slt = 0.0d0
         cngrp_slv = 0.0d0
         self_slt = 0.0d0
         self_slv = 0.0d0
         ucoul_slt=0.0D0
         ucoul_slv=0.0D0
         ucoul_ss=0.0D0
         CALL zeroa(fpx,fpy,fpz,ntap,1)
         CALL zeroa(Edx,Edy,Edz,ntap,1)
         CALL zeroa(Ex_rec,Ey_rec,Ez_rec,ntap,1)
         CALL zeroa(Ex_Cor,Ey_Cor,Ez_Cor,ntap,1)
      END IF

!=======================================================================
!     get forces: mts_forpp2 fnbgrp2
!=======================================================================
      IF(.NOT. Gauss_Mode) THEN
         CALL mts_forces2('m',xpa,ypa,zpa,xpga,ypga,zpga,xpcma
     &        ,ypcma,zpcma,mapnl,U_Thole,ELJ,fpx,fpy,fpz,do_LJ,flj_x
     &        ,flj_y,flj_z,worka,cpu,ncpu
     &        ,nstart,nend,nstart_a,nend_a,nlocal_a,node,nprocs,ncube
     &        ,' ',Udirect,ucoul_slt,ucoul_slv,ucoul_ss,uconf_slt
     &        ,uconf_slv,uconf_ss,Edx,Edy,Edz,charges,dipole,alphal
     &        ,skip_direct)
         CALL fnbgrp2(ss_index,xp0,yp0,zp0,charges,ecc12,ecc6,nbtype
     &        ,type_table,m6,alphal,ingrp_ef,ingrpp_ef,ingrp_ef_x,fpx
     &        ,fpy,fpz,do_LJ,flj_x,flj_y,flj_z,coulgrp_slt,coulgrp_slv
     &        ,confgrp_slt,confgrp_slv,plrzbij,Udirect,U_Thole,Edx,Edy
     &        ,Edz,dipole)
      END IF

#if defined PARALLEL
      IF(nprocs .GT. 1) THEN
         CALL P_merge_r8(Ucos,node,nprocs,ncube,rbyte)
         CALL P_merge_r8(Ucop,node,nprocs,ncube,rbyte)
         CALL P_merge_r8(Ucosp,node,nprocs,ncube,rbyte)
         CALL P_merge_r8(Udirect,node,nprocs,ncube,rbyte)
         CALL P_merge_r8(cngrp_slt,node,nprocs,ncube,rbyte)
         CALL P_merge_r8(cngrp_slv,node,nprocs,ncube,rbyte)

         CALL P_fold_r8(fpx,work,nstart_a,nend_a,nlocal_a,node,nprocs
     &        ,ncube,rbyte,nbyte)
         CALL P_fold_r8(fpy,work,nstart_a,nend_a,nlocal_a,node,nprocs
     &        ,ncube,rbyte,nbyte)
         CALL P_fold_r8(fpz,work,nstart_a,nend_a,nlocal_a,node,nprocs
     &        ,ncube,rbyte,nbyte)

         CALL P_fold_r8(Edx,work,nstart_a,nend_a,nlocal_a,node,nprocs
     &        ,ncube,rbyte,nbyte)
         CALL P_fold_r8(Edy,work,nstart_a,nend_a,nlocal_a,node,nprocs
     &        ,ncube,rbyte,nbyte)
         CALL P_fold_r8(Edz,work,nstart_a,nend_a,nlocal_a,node,nprocs
     &        ,ncube,rbyte,nbyte)

         CALL P_expand_r8(fpx,nstart_a,nend_a,nlocal_a,node,nprocs,ncube
     &        ,rbyte,nbyte)
         CALL P_expand_r8(fpy,nstart_a,nend_a,nlocal_a,node,nprocs,ncube
     &        ,rbyte,nbyte)
         CALL P_expand_r8(fpz,nstart_a,nend_a,nlocal_a,node,nprocs,ncube
     &        ,rbyte,nbyte)

         CALL P_expand_r8(Edx,nstart_a,nend_a,nlocal_a,node,nprocs,ncube
     &        ,rbyte,nbyte)
         CALL P_expand_r8(Edy,nstart_a,nend_a,nlocal_a,node,nprocs,ncube
     &        ,rbyte,nbyte)
         CALL P_expand_r8(Edz,nstart_a,nend_a,nlocal_a,node,nprocs,ncube
     &        ,rbyte,nbyte)

         CALL P_expand_r8(ene,nstart_a,nend_a,nlocal_a,node,nprocs,ncube
     &        ,rbyte,nbyte)
         CALL P_expand_r8(enedip,nstart_a,nend_a,nlocal_a,node,nprocs
     &        ,ncube,rbyte,nbyte)
      END IF
#endif

      Ugrp = coulgrp_slt+coulgrp_slv

      ucoul_slv=ucoul_slv+coulgrp_slv

      ucoul_slt=ucoul_slt+coulgrp_slt

      IF(do_LJ) THEN
         ELJ=ELJ+confgrp_slv+confgrp_slt
         uconf_slv=uconf_slv+confgrp_slv
         uconf_slt=uconf_slt+confgrp_slt
      END IF
!=======================================================================
!     get forces: pme + fnb14 + ferrf
!=======================================================================
! With pme for dipoles call always furier2

      IF(.NOT. Direct_Only) THEN
         CALL mts_furier2(node,nodex,nodey,nodez,ictxt,npy,npz,descQ
     &        ,nprocs,ncube,nstart_a,nend_a,nlocal_a,nstart_a,nend_a
     &        ,nlocal_a,xp0,yp0,zp0,xpa,ypa,zpa,xpcma,ypcma,zpcma,urcsp
     &        ,urcs,urcp,fpx,fpy,fpz,Ex_rec,Ey_rec,Ez_rec,Ex_Cor,Ey_Cor
     &        ,Ez_Cor,ene_rec,eer,fudgec,tag_bndg,ene_ferrf,charges
     &        ,dipole)
         eer = eer + ene_ferrf
         CALL cself_dipole(ntap,alphal,alphal,charges,dipole,volume
     &        ,uself,uself_dip)
c$$$      ELSE
c$$$         CALL zeroa(fdx,fdy,fdz,ntap,1)
c$$$         CALL mts_furier2(node,nodex,nodey,nodez,ictxt,npy,npz,descQ
c$$$     &        ,nprocs,ncube,nstart_a,nend_a,nlocal_a,nstart_a,nend_a
c$$$     &        ,nlocal_a,xp0,yp0,zp0,xpa,ypa,zpa,xpcma,ypcma,zpcma,urcsp
c$$$     &        ,urcs,urcp,fdx,fdy,fdz,Ex_rec,Ey_rec,Ez_rec,Ex_Cor,Ey_Cor
c$$$     &        ,Ez_Cor,ene_rec,eer,fudgec,tag_bndg,ene_ferrf,charges
c$$$     &        ,dipole)
c$$$         eer = eer + ene_ferrf
c$$$         CALL cself_dipole(ntap,alphal,alphal,charges,dipole,volume
c$$$     &        ,uself,uself_dip)
      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
