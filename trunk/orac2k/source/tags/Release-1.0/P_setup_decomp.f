      SUBROUTINE P_setup_decomp(node,nprocs,ncube,rbyte,nbyte,nstart_h
     &     ,nend_h,nlocal_h,nstart_ah,nend_ah,nlocal_ah,nstart_1,nend_1
     &     ,nlocal_1,nstart_ex0,nend_ex0,nlocal_ex0,nstart_ex,nend_ex
     &     ,nlocal_ex,nstart_2,nend_2,nlocal_2,nstart_cmi,nend_cmi
     &     ,nlocal_cmi,nstart_cm,nend_cm,nlocal_cm,nstart_cme,nend_cme
     &     ,nlocal_cme,nstart_cme0,nend_cme0,nlocal_cme0,nstart_g1
     &     ,nend_g1,nlocal_g1,nstart_gex,nend_gex,nlocal_gex,nstart_gex0
     &     ,nend_gex0,nlocal_gex0,ptr_ex,ptr_cme,ptr_gex,ntot_1,ntot_cmi
     &     ,ntot_ex,ntot_cme,nprot,protl,cnst_protp_1,cnst_protl_1
     &     ,cnst_protp_ex,cnst_protl_ex,cnstp,cnst_protp,cnst_protl
     &     ,lstrtch,lstretch,lbndg,lbend,cpu_h,ncpu_h,atomp,atomg
     &     ,ntap,ngrp,grppt,worka,xpa,ypa,zpa,xpga,ypga,zpga)

************************************************************************
*   Time-stamp: <95/01/07 00:47:53 marchi>                             *
*                                                                      *
*                                                                      *
*   Set up the atomic decomposition boundaries                         *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Feb 25 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER node,nprocs,ncube,rbyte,nbyte,nstart_h,nend_h,nlocal_h
     &     ,nstart_ah,nend_ah,nlocal_ah,nstart_1,nend_1,nlocal_1
     &     ,nstart_ex0,nend_ex0,nlocal_ex0,nstart_ex,nend_ex,nlocal_ex
     &     ,nstart_2,nend_2,nlocal_2,nstart_cmi,nend_cmi,nlocal_cmi
     &     ,nstart_cm,nend_cm,nlocal_cm,nstart_cme,nend_cme,nlocal_cme
     &     ,nstart_cme0,nend_cme0,nlocal_cme0,nstart_g1
     &     ,nend_g1,nlocal_g1,nstart_gex,nend_gex,nlocal_gex,nstart_gex0
     &     ,nend_gex0,nlocal_gex0,ptr_ex,ptr_cme,ptr_gex,ntot_1
     &     ,ntot_cmi,ntot_ex,ntot_cme,ntap,ngrp,nprot,ncpu_h,mpp

      INTEGER protl(*),cnst_protp_1,cnst_protl_1(*),cnst_protp_ex
     &     ,cnst_protl_ex(*),cnstp,cnst_protp,cnst_protl(*),lstrtch(2,*)
     &     ,lstretch,lbndg(3,*),lbend,atomp(*),atomg(*),grppt(2,*)
     &     ,worka(*)

      REAL*8 xpa(*),ypa(*),zpa(*),xpga(*),ypga(*),zpga(*),cpu_h(*)
      CHARACTER*1 P_shell

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER mapnl,mapdn(2),nmapdn,iret,mpp_a,M_get_length
      REAL*8  ucns_p,ucos_p,virs_p,virsp_p,ucnp_p,ucop_p,ucnsp_p,ucosp_p
     &     ,fpx_p,fpy_p,fpz_p,stressd_p(3,3),xpcma,ypcma,zpcma
      CHARACTER*80 errmsg
      REAL(8), DIMENSION (:), POINTER ::  phi

*----------------------- EXECUTABLE STATEMENTS ------------------------*


*=======================================================================
*---  Steve Plimpton's decomposition boundaries                     ----
*---   Atomic decomposition:                                        ----
*---   nstart_h: Used by the mts_force routine                      ----
*---   nstart_1: Used by the mts_intra routines                     ----
*---   nstart_2: Used in the h,l,m loops (forces from mts_forces are----
*---           expanded!)                                           ----
*---   nstart_3: Used by mts_furier to expand the bspline           ----
*---           coefficients of the PME calculation                  ----
*---   nstart_ex: Used for ballistic atoms in the n0, n1 loops      ----
*---            (no force from loop acting except constraints)      ----
*---                                                                ----
*---   Group decomposition:                                         ----
*---   nstart_cm : Used in the h,l,m loops                          ----
*---   nstart_cmi: Used in the n0,n1 loops                          ----
*---   nstart_cme: Used in the n0,n1 loops for ballistic groups     ----
*---                                                                ----
*=======================================================================

      iret=0
      ALLOCATE(phi(ntap))
#if !defined PARALLEL

      nstart_h=1
      nend_h=ngrp
      nlocal_h=ngrp
      nstart_ah=1
      nend_ah=ntap
      nlocal_ah=ntap
#else

*----- Fake split and setup
      CALL P_split_scalar0(node,nprocs,nstart_h,nend_h,nlocal_h
     &     ,nstart_ah,nend_ah,nlocal_ah,ngrp,grppt)
*----- Compute neighbors

      P_shell='l'
      CALL mts_forces('z',xpa,ypa,zpa,xpga,ypga,zpga,xpcma,ypcma,zpcma
     &     ,mapnl,mapdn,nmapdn,ucns_p,ucos_p,virs_p,virsp_p,ucnp_p
     &     ,ucop_p,ucnsp_p,ucosp_p,fpx_p,fpy_p,fpz_p,phi,stressd_p,worka
     &     ,cpu_h,ncpu_h,nstart_h,nend_h,nstart_ah,nend_ah,nlocal_ah
     &     ,node,nprocs,ncube,P_shell)
      
*----- Now can split and set up correctly (nstart_h)

      cpu_h(1:nprocs)=0.0D0
      CALL P_split_scalar(node,nprocs,ncube,nbyte,nstart_h,nend_h
     &     ,nlocal_h,nstart_ah,nend_ah,nlocal_ah,ngrp,grppt,worka,cpu_h
     &     ,ncpu_h)
#endif

*----- Split bonded interaction and constraints (nstart_1, nstart_ex)

      CALL P_atoms_split_intra(node,nprocs,ncube,nbyte,nstart_1,nend_1
     &     ,nlocal_1,nstart_ex0,nend_ex0,nlocal_ex0,nstart_ex,nend_ex
     &     ,nlocal_ex,nprot,protl,cnstp,cnst_protp,cnst_protl
     &     ,cnst_protp_1,cnst_protl_1,cnst_protp_ex,cnst_protl_ex,ntap
     &     ,lstrtch,lstretch,lbndg,lbend,atomp,iret,errmsg)

*----- Split intramolecular interaction arrays

      CALL P_split_intra(node,nprocs,ncube,nstart_1,nend_1,worka,iret
     &     ,errmsg)
      IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)

*----- Split nonbonded interaction loop and do set up (nstart_2) 

      CALL P_atoms_split_inter(node,nprocs,ncube,nbyte,ntap,nstart_2
     &     ,nend_2,nlocal_2,cnstp,cnst_protp,cnst_protl,iret,errmsg)

*----- Split constraint interaction arrays according to nstart_2

      CALL P_split_constr(node,nprocs,ncube,nstart_2,nend_2,worka,iret
     &     ,errmsg)
      IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)

*----- Split nonbonded loop, bonded loop and ballistic groups 

      CALL P_cm_split(node,nstart_2,nend_2,nlocal_2,nstart_cm,nend_cm
     &     ,nlocal_cm,atomp)
      CALL P_cm_split(node,nstart_1,nend_1,nlocal_1,nstart_cmi,nend_cmi
     &     ,nlocal_cmi,atomp)
      CALL P_cm_split(node,nstart_ex,nend_ex,nlocal_ex,nstart_cme
     &     ,nend_cme,nlocal_cme,atomp)

*----- Split atom groups cm for nonbonded, intramolecular and ballistic

      IF(nstart_ex .EQ. ntap+1) nstart_cme=nprot+1
      nstart_cme0=nstart_cme
      nend_cme0=nend_cme
      nlocal_cme0=nlocal_cme

*----- Broadcast pointer of the first ballistic atom

      ptr_ex=nstart_ex
#if defined PARALLEL
      CALL P_broadcast_i(ptr_ex,1)
#endif      
      IF(nstart_cme .NE. nprot+1) THEN
         ptr_cme=atomp(ptr_ex)
         ntot_1=ptr_ex-1
         ntot_cmi=ptr_cme-1
         ntot_ex=ntap-ptr_ex+1
         ntot_cme=nprot-ptr_cme+1
         nstart_gex0=atomg(nstart_ex0)
         nend_gex0=atomg(nend_ex0)
         nlocal_gex0=nend_gex0-nstart_gex0+1
         nstart_gex=atomg(nstart_ex)
         nend_gex=atomg(nend_ex)
         nlocal_gex=nend_gex-nstart_gex+1
         ptr_gex=atomg(ptr_ex)
#ifdef PARALLEL
         CALL P_adjust_decomp(nstart_gex,nend_gex,nlocal_gex,node,nbyte
     &        ,nprocs)
         CALL P_adjust_decomp(nstart_gex0,nend_gex0,nlocal_gex0,node
     &        ,nbyte,nprocs)
#endif
      ELSE
         ptr_cme=nprot+1
         ntot_1=ntap
         ntot_cmi=nprot
         ntot_ex=0
         ntot_cme=0
      END IF
      IF(nstart_1 .NE. 0) THEN
         nstart_g1=atomg(nstart_1)
         nend_g1=atomg(nend_1)
      ELSE
         nstart_g1=1
         nend_g1=0
      END IF
      nlocal_g1=nend_g1-nstart_g1+1

#ifdef PARALLEL
      IF(nstart_g1 .NE. 0) CALL P_adjust_decomp(nstart_g1,nend_g1
     &     ,nlocal_g1,node,nbyte,nprocs)
#endif

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
