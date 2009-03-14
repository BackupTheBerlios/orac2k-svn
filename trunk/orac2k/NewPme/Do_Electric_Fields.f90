SUBROUTINE Do_Electric_Fields(xp0,yp0,zp0,xpa,ypa,zpa,xpga,ypga,zpga&
     &,xpcma,ypcma,zpcma,mapnl,ingrpp_ef,ingrp_ef,ingrp_ef_x,fpx,fpy&
     &,fpz,flj_x,flj_y,flj_z,Ex_rec,Ey_rec,Ez_rec,Edx,Edy,Edz&
     &,Ex_Cor,Ey_Cor,Ez_Cor,worka,cpu,ncpu,nstart,nend&
     &,nlocal,nstart_a,nend_a,nlocal_a,node,nodex,nodey,nodez,ictxt&
     &,npy,npz,descQ,nprocs,ncube,tag_bndg,Ugrp,eer,uself,uself_dip&
     &,Udirect,ene_ferrf,U_Thole,ELJ,ucoul_slt,ucoul_slv,ucoul_ss&
     &,uconf_slt,uconf_slv,uconf_ss,charges,dipole,ntap,What_To_Do)

!!$***********************************************************************
!!$   Time-stamp: <2002-09-27 10:49:30 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Mon Jun 28 2004 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*

  IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

  INTEGER :: node,nprocs,ncube,mapnl(*),ingrpp_ef,ingrp_ef(2,*)
  INTEGER :: nodex,nodey,nodez,ictxt,npy,npz,descQ(*),ingrp_ef_x(*)
  INTEGER :: tag_bndg(*)
  REAL(8)  xpa(*),ypa(*),zpa(*),xpcma(*),ypcma(*),zpcma(*),xpga(*)&
       &,ypga(*),zpga(*),xp0(*),yp0(*),zp0(*)
  REAL(8)  Edx(*),Edy(*),Edz(*),dipole(3,*),charges(*)
  REAL(8)  Ex_rec(*),Ey_rec(*),Ez_rec(*)
  REAL(8)  Ex_Cor(*),Ey_Cor(*),Ez_Cor(*)
  REAL(8)  fpx(*),fpy(*),fpz(*),cpu,Ugrp,Udirect,eer,uself,ene_ferrf&
       &,uself_dip,U_Thole,ELJ,flj_x(*),flj_y(*),flj_z(*),ucoul_slt&
       &,ucoul_slv,ucoul_ss,uconf_slt,uconf_slv,uconf_ss 
  INTEGER :: nstart,nend,nlocal,nstart_a,nend_a,nlocal_a,worka(*)
  INTEGER :: npp_m,ncpu,ntap
  CHARACTER(80) :: What_To_Do
  

!!$----------------------- VARIABLES IN COMMON --------------------------*
  INTERFACE
     SUBROUTINE get_electric_field(xp0,yp0,zp0,xpa,ypa,zpa,xpga,ypga&
          &,zpga,xpcma,ypcma,zpcma,mapnl,ingrpp_ef,ingrp_ef&
          &,ingrp_ef_x,fpx,fpy,fpz,flj_x,flj_y,flj_z,Ex_rec&
          &,Ey_rec,Ez_rec,Edx,Edy,Edz,Ex_Cor,Ey_Cor,Ez_Cor,worka,cpu&
          &,ncpu,nstart,nend,nlocal,nstart_a,nend_a&
          &,nlocal_a,node,nodex,nodey,nodez,ictxt,npy,npz,descQ&
          &,nprocs,ncube,tag_bndg,Ugrp,eer,uself,uself_dip,Udirect&
          &,ene_ferrf,U_Thole,ELJ,ucoul_slt,ucoul_slv,ucoul_ss&
          &,uconf_slt,uconf_slv,uconf_ss,charges,dipole,do_LJ&
          &,skip_direct,Gauss,Direct_Only,set_zero)
       IMPLICIT none

       INTEGER :: node,nprocs,ncube,mapnl(*),ingrpp_ef,ingrp_ef(2,*)
       INTEGER :: nodex,nodey,nodez,ictxt,npy,npz,descQ(*)&
            &,ingrp_ef_x(*) 
       INTEGER :: tag_bndg(*) 
       REAL(8) ::  xpa(*),ypa(*),zpa(*),xpcma(*),ypcma(*),zpcma(*)&
            &,xpga(*),ypga(*),zpga(*),xp0(*),yp0(*),zp0(*) 
       REAL(8) ::  Edx(*),Edy(*),Edz(*),dipole(3,*),charges(*)
       REAL(8) ::  Ex_rec(*),Ey_rec(*),Ez_rec(*)
       REAL(8) ::  Ex_Cor(*),Ey_Cor(*),Ez_Cor(*)
       REAL(8) ::  fpx(*),fpy(*),fpz(*),cpu,Ugrp,Udirect,eer,uself&
            &,ene_ferrf,uself_dip,U_Thole,ELJ,flj_x(*),flj_y(*),flj_z(&
            &*),ucoul_slt,ucoul_slv,ucoul_ss,uconf_slt,uconf_slv&
            &,uconf_ss 
       INTEGER :: nstart,nend,nlocal,nstart_a,nend_a,nlocal_a&
            &,worka(*) 
       INTEGER :: npp_m,ncpu 
       LOGICAL :: do_LJ,skip_direct,do_CD,Gauss,Direct_Only,set_zero
     END SUBROUTINE get_electric_field
  END INTERFACE

!!$------------------------- LOCAL VARIABLES ----------------------------*

  LOGICAL :: do_LJ,skip_direct,Gauss,Direct_Only,set_zero
  REAL(8), DIMENSION (:,:), ALLOCATABLE, SAVE :: dipole_1
  REAL(8), DIMENSION (:), ALLOCATABLE, SAVE :: charges_1
  INTEGER, SAVE :: First_Time=0 

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

  IF(First_Time == 0) THEN
     First_Time=1
     ALLOCATE(dipole_1(3,ntap))
     ALLOCATE(charges_1(ntap))
  END IF
  
  charges_1=charges(1:ntap)
  dipole_1=dipole(1:3,1:ntap)
  set_zero=.TRUE.
  Gauss=.FALSE.
  Direct_only=.FALSE.

  SELECT CASE(What_to_Do)
  CASE('Charges Direct and LJ')
     do_LJ=.TRUE.
     skip_direct=.FALSE.
     Direct_Only=.TRUE.
     dipole_1=0.0D0
     CALL get_electric_field(xp0,yp0,zp0,xpa,ypa,zpa,xpga,ypga&
          &,zpga,xpcma,ypcma,zpcma,mapnl,ingrpp_ef,ingrp_ef&
          &,ingrp_ef_x,fpx,fpy,fpz,flj_x,flj_y,flj_z,Ex_rec&
          &,Ey_rec,Ez_rec,Edx,Edy,Edz,Ex_Cor,Ey_Cor,Ez_Cor,worka,cpu&
          &,ncpu,nstart,nend,nlocal,nstart_a,nend_a& 
          &,nlocal_a,node,nodex,nodey,nodez,ictxt,npy,npz,descQ&
          &,nprocs,ncube,tag_bndg,Ugrp,eer,uself,uself_dip,Udirect&
          &,ene_ferrf,U_Thole,ELJ,ucoul_slt,ucoul_slv,ucoul_ss&
          &,uconf_slt,uconf_slv,uconf_ss,charges_1,dipole_1,do_LJ&
          &,skip_direct,Gauss,Direct_Only,set_zero)
  CASE('Charges Reciprocal no LJ')
     do_LJ=.FALSE.
     skip_direct=.TRUE.
     Direct_Only=.FALSE.
     dipole_1=0.0D0
     CALL get_electric_field(xp0,yp0,zp0,xpa,ypa,zpa,xpga,ypga&
          &,zpga,xpcma,ypcma,zpcma,mapnl,ingrpp_ef,ingrp_ef&
          &,ingrp_ef_x,fpx,fpy,fpz,flj_x,flj_y,flj_z,Ex_rec&
          &,Ey_rec,Ez_rec,Edx,Edy,Edz,Ex_Cor,Ey_Cor,Ez_Cor,worka,cpu&
          &,ncpu,nstart,nend,nlocal,nstart_a,nend_a&
          &,nlocal_a,node,nodex,nodey,nodez,ictxt,npy,npz,descQ&
          &,nprocs,ncube,tag_bndg,Ugrp,eer,uself,uself_dip,Udirect&
          &,ene_ferrf,U_Thole,ELJ,ucoul_slt,ucoul_slv,ucoul_ss&
          &,uconf_slt,uconf_slv,uconf_ss,charges_1,dipole_1,do_LJ&
          &,skip_direct,Gauss,Direct_Only,set_zero)
  CASE('Charges Full and LJ')
     do_LJ=.TRUE.
     skip_direct=.FALSE.
     Direct_Only=.FALSE.
     dipole_1=0.0D0
     CALL get_electric_field(xp0,yp0,zp0,xpa,ypa,zpa,xpga,ypga&
          &,zpga,xpcma,ypcma,zpcma,mapnl,ingrpp_ef,ingrp_ef&
          &,ingrp_ef_x,fpx,fpy,fpz,flj_x,flj_y,flj_z,Ex_rec&
          &,Ey_rec,Ez_rec,Edx,Edy,Edz,Ex_Cor,Ey_Cor,Ez_Cor,worka,cpu&
          &,ncpu,nstart,nend,nlocal,nstart_a,nend_a&
          &,nlocal_a,node,nodex,nodey,nodez,ictxt,npy,npz,descQ&
          &,nprocs,ncube,tag_bndg,Ugrp,eer,uself,uself_dip,Udirect&
          &,ene_ferrf,U_Thole,ELJ,ucoul_slt,ucoul_slv,ucoul_ss&
          &,uconf_slt,uconf_slv,uconf_ss,charges_1,dipole_1,do_LJ&
          &,skip_direct,Gauss,Direct_Only,set_zero)
  CASE('Full and LJ')
     do_LJ=.TRUE.
     skip_direct=.FALSE.
     Direct_Only=.FALSE.
     CALL get_electric_field(xp0,yp0,zp0,xpa,ypa,zpa,xpga,ypga&
          &,zpga,xpcma,ypcma,zpcma,mapnl,ingrpp_ef,ingrp_ef&
          &,ingrp_ef_x,fpx,fpy,fpz,flj_x,flj_y,flj_z,Ex_rec&
          &,Ey_rec,Ez_rec,Edx,Edy,Edz,Ex_Cor,Ey_Cor,Ez_Cor,worka,cpu&
          &,ncpu,nstart,nend,nlocal,nstart_a,nend_a&
          &,nlocal_a,node,nodex,nodey,nodez,ictxt,npy,npz,descQ&
          &,nprocs,ncube,tag_bndg,Ugrp,eer,uself,uself_dip,Udirect&
          &,ene_ferrf,U_Thole,ELJ,ucoul_slt,ucoul_slv,ucoul_ss&
          &,uconf_slt,uconf_slv,uconf_ss,charges_1,dipole_1,do_LJ&
          &,skip_direct,Gauss,Direct_Only,set_zero)
  CASE('Full no LJ')
     do_LJ=.FALSE.
     skip_direct=.FALSE.
     Direct_Only=.FALSE.
     
     CALL get_electric_field(xp0,yp0,zp0,xpa,ypa,zpa,xpga,ypga&
          &,zpga,xpcma,ypcma,zpcma,mapnl,ingrpp_ef,ingrp_ef&
          &,ingrp_ef_x,fpx,fpy,fpz,flj_x,flj_y,flj_z,Ex_rec&
          &,Ey_rec,Ez_rec,Edx,Edy,Edz,Ex_Cor,Ey_Cor,Ez_Cor,worka,cpu&
          &,ncpu,nstart,nend,nlocal,nstart_a,nend_a&
          &,nlocal_a,node,nodex,nodey,nodez,ictxt,npy,npz,descQ&
          &,nprocs,ncube,tag_bndg,Ugrp,eer,uself,uself_dip,Udirect&
          &,ene_ferrf,U_Thole,ELJ,ucoul_slt,ucoul_slv,ucoul_ss&
          &,uconf_slt,uconf_slv,uconf_ss,charges_1,dipole_1,do_LJ&
          &,skip_direct,Gauss,Direct_Only,set_zero)
  CASE('Gauss no LJ')
     do_LJ=.FALSE.
     skip_direct=.TRUE.
     Gauss=.TRUE.

     CALL get_electric_field(xp0,yp0,zp0,xpa,ypa,zpa,xpga,ypga&
          &,zpga,xpcma,ypcma,zpcma,mapnl,ingrpp_ef,ingrp_ef&
          &,ingrp_ef_x,fpx,fpy,fpz,flj_x,flj_y,flj_z,Ex_rec&
          &,Ey_rec,Ez_rec,Edx,Edy,Edz,Ex_Cor,Ey_Cor,Ez_Cor,worka,cpu&
          &,ncpu,nstart,nend,nlocal,nstart_a,nend_a&
          &,nlocal_a,node,nodex,nodey,nodez,ictxt,npy,npz,descQ&
          &,nprocs,ncube,tag_bndg,Ugrp,eer,uself,uself_dip,Udirect&
          &,ene_ferrf,U_Thole,ELJ,ucoul_slt,ucoul_slv,ucoul_ss&
          &,uconf_slt,uconf_slv,uconf_ss,charges_1,dipole_1,do_LJ&
          &,skip_direct,Gauss,Direct_Only,set_zero)  
     
  CASE('Gauss and LJ')
     do_LJ=.TRUE.
     skip_direct=.TRUE.
     Gauss=.TRUE.

     CALL get_electric_field(xp0,yp0,zp0,xpa,ypa,zpa,xpga,ypga&
          &,zpga,xpcma,ypcma,zpcma,mapnl,ingrpp_ef,ingrp_ef&
          &,ingrp_ef_x,fpx,fpy,fpz,flj_x,flj_y,flj_z,Ex_rec&
          &,Ey_rec,Ez_rec,Edx,Edy,Edz,Ex_Cor,Ey_Cor,Ez_Cor,worka,cpu&
          &,ncpu,nstart,nend,nlocal,nstart_a,nend_a&
          &,nlocal_a,node,nodex,nodey,nodez,ictxt,npy,npz,descQ&
          &,nprocs,ncube,tag_bndg,Ugrp,eer,uself,uself_dip,Udirect&
          &,ene_ferrf,U_Thole,ELJ,ucoul_slt,ucoul_slv,ucoul_ss&
          &,uconf_slt,uconf_slv,uconf_ss,charges_1,dipole_1,do_LJ&
          &,skip_direct,Gauss,Direct_Only,set_zero)  
     
  CASE default
     WRITE(*,*) 'Don''t know what to do!'
     STOP
  END SELECT

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE Do_Electric_Fields
