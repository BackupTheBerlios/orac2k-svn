      SUBROUTINE  mts_forces(rshell,xp0,yp0,zp0,xpg,ypg,zpg,xpcm,ypcm
     &     ,zpcm,mapnl,mapdn,nmapdn,ucns,ucos,virs,virsp,ucnp,ucop,ucnsp
     &     ,ucosp,fpx,fpy,fpz,phi,stress,worka
     &     ,cpu_h,ncpu_h,nstart,nend,nstart_a,nend_a,nlocal_a,node
     &     ,nprocs,ncube,P_shell)


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

*----------------------- ARGUMENTS -------------------------------------

      CHARACTER*1  rshell,P_shell
      REAL*8  xp0(*),yp0(*),zp0(*),fpx(*),fpy(*),fpz(*),xpg(*),ypg(*)
     &     ,zpg(*),xpcm(*),ypcm(*),zpcm(*),phi(*),cpu_h(*),stress(3,3)
      REAL*8  ucnp,ucop,ucns,ucos,ucnsp,ucosp,virs,virsp
      INTEGER ncpu_h,nstart,nend,node,nprocs,ncube,worka(*)
      INTEGER mapnl(*),mapdn(2,*),nmapdn(*),nstart_a,nend_a,nlocal_a
      

*-------------------- DEFINITION OF AN EXTERNAL FUNCTION ---------------

#ifdef PARALLEL
      INCLUDE 'mpif.h'
      include 'mpi_size.h'
#endif
      EXTERNAL  near0
      LOGICAL near0

*-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8 ro,rto,ri,rti,rni
      INTEGER iz1,m,n,i,npp,iret,nlocal_b,nend_b,nstart_b,ierr
      REAL*8 vfcp_ff,tfcp_ff,elapse,tdelta_ff,gcpu_ff
      REAL*8 vfcp_ff1,tfcp_ff1,tdelta_ff1,gcpu_ff1
      CHARACTER*80 errmsg
      LOGICAL flag

*==================== EXECUTABLE STATEMENTS ============================

*=======================================================================
*     select cutoffs according to current value of SHELL
*=======================================================================

      gcpu_ff1=0.0D0
      CALL timer(vfcp_ff1,tfcp_ff1,elapse)
      tdelta_ff1=tfcp_ff1
      flag=.TRUE.
      if(rshell.eq.'z') then
         iz1=0
         ro=rcuth
         ri=rcuth
         rto=0.d0
         rti=rtolh
         rni=rneih

         CALL zero0(worka,ngrp)
      endif
      if(rshell.eq.'u') then
         iz1=0
         ro=rcuth
         ri=rcuth
         rto=0.d0
         rti=rtolh
         rni=rneih

         CALL zero0(worka,ngrp)
      endif
      if(rshell.eq.'v') then
         iz1=0
         ro=rcuth
         ri=rcuth
         rto=0.d0
         rti=rtolh
         rni=rneih
         CALL zero0(worka,ngrp)
         flag=.FALSE.
      endif
      if(rshell.eq.'h') then
         iz1=1
         if(rneih.lt.1.d-10) iz1=0
         ro=rcuth
         ri=rcutl
         rto=rtolh
         rti=rtoll
         rni=rneil
      endif
      if(rshell.eq.'l') then
         iz1=1
         ro=rcutl
         rto=rtoll
         ri=rcutm
         rti=rtolm
         rni=rneim
      end if   
      if(rshell.eq.'m') then
         iz1=1
         ro=rcutm
         rto=rtolm
         ri=0.d0
         rti=0.d0
         rni=0.d0
      end if   
      ucnp=0.d0
      ucop=0.d0
      ucns=0.d0
      ucos=0.d0
      virs=0.d0
      ucnsp=0.d0
      ucosp=0.d0
      DO m=1,3
         DO n=1,3
            stress(n,m)=0.0D0
         END DO
      END DO

      gcpu_ff=0.0D0
      CALL timer(vfcp_ff,tfcp_ff,elapse)
      tdelta_ff=tfcp_ff
      CALL mts_forpp(rshell,ss_index,xp0,yp0,zp0,xpg,ypg,zpg,chrge
     &     ,nbtype,type_table,m6,ntap,atomg,xpcm,ypcm,zpcm,groupp,atomp
     &     ,co,ecc12,ecc6,clewld,alphal,erfc_spline,erfc_bin,erfc_arr
     &     ,mapnl,ngrp,grppt,ucnp,ucns,ucnsp,ucop,ucos,ucosp,fpx,fpy,fpz
     &     ,phi,stress,rni,ri,ro,rti,rto,pmass,iz1,worka,nstart,nend
     &     ,Boltz_Group,flag,lpnbd,iret,errmsg)

      CALL timer(vfcp_ff,tfcp_ff,elapse)
      tdelta_ff=tfcp_ff-tdelta_ff
      gcpu_ff=gcpu_ff + tdelta_ff      

#ifdef PARALLEL
      CALL P_get_errmsg(iret,errmsg,80,node,nprocs,ncube,nbyte)
      CALL MPI_ALLGATHER(gcpu_ff,1,MPI_REAL8,cpu_h,1,MPI_REAL8
     &     ,MPI_COMM_WORLD,ierr)
      ncpu_h=ncpu_h+1
#endif

#ifdef PARALLEL
      IF(iz1 .NE. 0 .AND. ri .NE. ro .AND. nprocs .GT. 1) THEN
         CALL P_merge_r8(ucnp)
         CALL P_merge_r8(ucns)
         CALL P_merge_r8(ucnsp)
         CALL P_merge_r8(ucop)
         CALL P_merge_r8(ucos)
         CALL P_merge_r8(ucosp)
         CALL P_merge_vecr8(stress,9)
      END IF
#endif
      npp=SUM(worka(nstart:nend))+nend-nstart+1
c$$$      CALL P_Write_Neighbors(rshell,npp,kprint,node,nprocs)
      IF(iz1 .EQ. 0) THEN
#ifdef PARALLEL
         IF(nprocs .GT. 1) CALL P_merge_i(npp)
#endif
         WRITE(kprint,10000) npp
      END IF

      CALL timer(vfcp_ff1,tfcp_ff1,elapse)
      tdelta_ff1=tfcp_ff1-tdelta_ff1
      gcpu_ff1=gcpu_ff1 + tdelta_ff1

*================= END OF EXECUTABLE STATEMENTS ========================

10000 FORMAT('Neighbor Lists Dimensions     *neighbor(',i7,')* '/)
      RETURN
      END
