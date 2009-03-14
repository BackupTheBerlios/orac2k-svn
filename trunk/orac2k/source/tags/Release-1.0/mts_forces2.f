      SUBROUTINE  mts_forces2(rshell,xp0,yp0,zp0,xpg,ypg,zpg,xpcm,ypcm
     &     ,zpcm,mapnl,U_Thole,ELJ,fpx,fpy,fpz,do_LJ,flj_x,flj_y,flj_z
     &     ,worka,cpu_h,ncpu_h,nstart,nend
     &     ,nstart_a,nend_a,nlocal_a,node,nprocs,ncube,P_shell,Utotal
     &     ,ucoul_slt,ucoul_slv,ucoul_ss,uconf_slt,uconf_slv,uconf_ss
     &     ,Edx,Edy,Edz,charges,dipole,aalphal,skip_direct)

************************************************************************
*                                                                      *
*     FORCES calls MTS_FORSLV,MTS_FORPW,MTS_FORPP,FINTRAPS             *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     Written by Piero Procacci, CECAM-ENS Lyon 1995                   *
*                                                                      *
*                                                                      *
************************************************************************


*======================= DECLARATIONS ==================================

      IMPLICIT none

*-------------------- VARIABLES IN COMMONS -----------------------------

      include 'parst.h'
      include 'cpropar.h'
      INCLUDE 'unit.h'
      REAL*8  work
      COMMON /rag1/ work(m1)

*----------------------- ARGUMENTS -------------------------------------

      CHARACTER*1  rshell,P_shell
      REAL*8  xp0(*),yp0(*),zp0(*),fpx(*),fpy(*),fpz(*),xpg(*),ypg(*)
     &     ,zpg(*),xpcm(*),ypcm(*),zpcm(*),cpu_h(*),flj_x(*),flj_y(*)
     &     ,flj_z(*)
      REAL*8  U_Thole,ELJ
      REAL*8 Utotal,Edx(*),Edy(*),Edz(*),dipole(3,*),charges(*)
     &     ,ucoul_slt,ucoul_slv,ucoul_ss,uconf_slt,uconf_slv,uconf_ss
     &     ,aalphal
      INTEGER ncpu_h,nstart,nend,node,nprocs,ncube,worka(*),npp_o
     &     ,npp_loc
      INTEGER mapnl(*),nstart_a,nend_a,nlocal_a
      LOGICAL do_LJ,skip_direct

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER m,n,i,npp,iret
      REAL*8 vfcp_ff,tfcp_ff,elapse,tdelta_ff,gcpu_ff
      REAL*8 vfcp_ff1,tfcp_ff1,tdelta_ff1,gcpu_ff1
      CHARACTER*80 errmsg

*==================== EXECUTABLE STATEMENTS ============================

*=======================================================================
*     select cutoffs according to current value of SHELL
*=======================================================================

      gcpu_ff1=0.0D0
      CALL timer(vfcp_ff1,tfcp_ff1,elapse)
      tdelta_ff1=tfcp_ff1

      U_Thole=0.0D0

      gcpu_ff=0.0D0
      CALL timer(vfcp_ff,tfcp_ff,elapse)
      tdelta_ff=tfcp_ff
 
      CALL mts_forpp2(ss_index,xp0,yp0,zp0,xpg,ypg,zpg,charges,nbtype
     &     ,type_table,m6,ntap,atomg,xpcm,ypcm,zpcm,groupp,atomp,co
     &     ,aalphal,mapnl,ngrp,grppt,ucoul_slt,ucoul_slv,ucoul_ss
     &     ,uconf_slt,uconf_slv,uconf_ss,U_Thole,ELJ,fpx,fpy,fpz,do_LJ
     &     ,flj_x,flj_y,flj_z,ecc12,ecc6,Rcut_El,Rcut_LJ,worka,nstart
     &     ,nend,iret,errmsg,plrzbij,Edx,Edy,Edz,dipole,skip_direct)

      Utotal=Utotal+ucoul_slt+ucoul_slv+ucoul_ss
      CALL timer(vfcp_ff,tfcp_ff,elapse)
      tdelta_ff=tfcp_ff-tdelta_ff
      gcpu_ff=gcpu_ff + tdelta_ff

#ifdef PARALLEL
      CALL P_get_errmsg(iret,errmsg,80,node,nprocs,ncube,nbyte)
      CALL P_merge_r8(ucoul_slt)
      CALL P_merge_r8(ucoul_slv)
      CALL P_merge_r8(ucoul_ss)
      CALL P_merge_r8(uconf_slt)
      CALL P_merge_r8(uconf_slv)
      CALL P_merge_r8(uconf_ss)
      CALL P_merge_r8(U_Thole)
      CALL P_merge_r8(ELJ)
#endif
      npp=SUM(worka(nstart:nend))+nend-nstart+1
      IF(rshell .EQ. P_shell) THEN
         DO i=1,nprocs
            cpu_h(i)=0.0D0
         END DO
         cpu_h(node+1)=gcpu_ff
#ifdef PARALLEL
         CALL P_merge_vecr8(cpu_h,nprocs,node,nprocs,ncube,rbyte)
#endif
      END IF

      CALL timer(vfcp_ff1,tfcp_ff1,elapse)
      tdelta_ff1=tfcp_ff1-tdelta_ff1
      gcpu_ff1=gcpu_ff1 + tdelta_ff1

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
