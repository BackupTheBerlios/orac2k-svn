      SUBROUTINE mts_furier(node,nodex,nodey,nodez,ictxt,npy,npz,descQ
     &     ,nprocs,ncube,nstart,nend,nlocal,nstart_2,nend_2,nlocal_2,xpb
     &     ,ypb,zpb,xp0,yp0,zp0,xpcm,ypcm,zpcm,urcsp,urcs,urcp,virsp
     &     ,virs,virp,fpx,fpy,fpz,phi,fsin14,fsbend,fsbond,fscnstr_slt
     &     ,fscnstr_slv,coul_bnd_slt,coul_bnd_slv,rsh,rshk,eer,virial
     &     ,fudgec,tag_bndg)

************************************************************************
*                                                                      *
*     Transit routine : calls MTS_FURIWW, MTS_FURIPP, MTS_FURIPW,      *
*     FERFF, WATSELF, PME_FFT                                          *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     Written by Piero Procacci, CECAM-ENS Lyon 1995                   *
*     Modified for parallel machines by                                *
*              Massimo Marchi                                          *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed Feb 10 1999 -                                     *
*                                                                      *
*                                                                      *
************************************************************************


*======================= DECLARATIONS ==================================

      Use Pme_Save
      USE REDUCE, ONLY: P_Reduce_Forces=>Reduce_Forces      
      IMPLICIT none
      INTERFACE
         SUBROUTINE P_Change_Decomposition(Decmp_name,nato,vpx,vpy,vpz
     &        ,nstart,nend,nstart1,nend1,node,nprocs)
         REAL(8) :: vpx(*),vpy(*),vpz(*)
         INTEGER :: nstart,nend,nstart1,nend1
         INTEGER :: node,nprocs,nato
         CHARACTER(*) :: Decmp_name
         END SUBROUTINE P_Change_Decomposition         
      END INTERFACE

*----------------------- ARGUMENTS -------------------------------------

      INTEGER tag_bndg(*),node,nodex,nodey,nodez,ictxt,npy
     &     ,npz,ncube,nprocs,descQ(*),ierr
      INTEGER nstart,nend,nlocal,nstart_2,nend_2,nlocal_2
      REAL*8     xp0(*),yp0(*),zp0(*),fpx(*),fpy(*),fpz(*),xpb(*),ypb(*)
     &     ,zpb(*),phi(*)
      REAL*8     urcsp,urcs,urcp,virsp,virs,virp,fsin14,gsin14,fsbend
     &     ,gsbend,fsbond,gsbond,fudgec,fscnstr_slt,fscnstr_slv
     &     ,coul_bnd_slt,coul_bnd_slv
      character*1 rsh,rshk
      REAL*8  xpcm(*),ypcm(*),zpcm(*)
      REAL*8   virial(3,3),eer,cpu_l(128)

c---  stuff for b-spline interpolation and FFT 

*---- VARIABLES IN INCLUDE --------------------------------------------*

      INCLUDE 'parst.h'
      INCLUDE 'cpropar.h'
#if defined PARALLEL
      INCLUDE 'mpif.h'
#endif
      REAL*8  fxx(m1,2),fyy(m1,2),fzz(m1,2),phiw(m1)
      SAVE fxx,fyy,fzz,phiw

      REAL*8  fsbond_slt,fsbend_slt,fsin14_slt,fsbond_slv,fsbend_slv
     &     ,fsin14_slv

*-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8  time1_avg,time2_avg,time_avg,time_a,time_b,dd,ela
      INTEGER ntime_avg
      INTEGER  i,j

#define _AXPY_  daxpy
      INTEGER n_loc,one,n_loc_2

      DATA one/1/
      DATA ntime_avg/0/
      DATA time_avg/0.0D0/

*==== EXECUTABLE STATEMENTS: ==========================================*

      n_loc=nlocal
      n_loc_2=nlocal_2

c---------------------------------------------------------------------
c     PROTEIN alone or PROTEIN + SOLVENT
c---------------------------------------------------------------------

      
      CALL timer(dd,time_a,ela)
      DO i=1,3
         DO j=1,3
            virial(i,j)=0.0D0
         END DO
      END DO
      eer=0.0D0
      IF(rsh .NE. rshk) RETURN

      IF(.NOT. pme) THEN
         CALL furipp(ss_index,oc,xp0,yp0,zp0,chrge,ntap,atomp,grppt
     &        ,alphal,rkcut,volume,urcp,urcs,urcsp,xpcm,ypcm,zpcm,fpx
     &        ,fpy,fpz,co,virial)
      ELSE

*=======================================================================
*--- Call fft_pme                                                   ----
*=======================================================================
      
         CALL fft_pme(node,nodex,nodey,nodez,ictxt,npy,npz,descQ,nprocs
     &        ,ncube,nbyte,rbyte,nstart_2,nend_2,nlocal_2,ntap,xp0,yp0
     &        ,zp0,xpcm,ypcm,zpcm,chrge,co,oc,volume,alphal,pme_order
     &        ,nfft1,nfft2,nfft3,nfft3_start,nfft3_local,nfft2_start
     &        ,nfft2_local,eer,fpx,fpy,fpz,phi,virial,atomp,grppt
     &        ,pressure,rkcut)
      END IF
#ifdef PARALLEL
      CALL P_Reduce_Forces(fpx,fpy,fpz)
#endif
         
      IF(remove_momentum) THEN
c$$$         CALL remove_mv(fpx,fpy,fpz,mass,ntap)
      END IF
      CALL zeroa(fxx(nstart,1),fyy(nstart,1),fzz(nstart,1),nlocal,1)
      CALL zeroa(fxx(nstart,2),fyy(nstart,2),fzz(nstart,2),nlocal,1)
      phiw=0.0D0

*=======================================================================
*      subtract "intramolecular term" in the ZERO cell:   BONDS     ----
*=======================================================================

      fsbond=0.0D0
      coul_bnd_slv=0.0D0
      coul_bnd_slt=0.0D0
      IF(lstretch .NE. 0) THEN
         CALL ferrf(ss_index,alphal,chrge,1.0D0,xpb,ypb,zpb,1,lstrtch
     &        ,lstretch,lbnd_x,fsbond_slt,fsbond_slv,fxx(1,2),fyy(1,2)
     &        ,fzz(1,2),phiw,erf_corr,erf_arr_corr,delew,rlew)
         fsbond=fsbond_slt
         coul_bnd_slv=fsbond_slv
         coul_bnd_slt=fsbond_slt
      END IF

*=======================================================================
*----- Subtract "intramolecular term" in the ZERO cell:   BENDS     ----
*=======================================================================

      fsbend=0.0D0
      IF(int13p .NE. 0) THEN
         CALL ferrf_tag(ss_index,alphal,chrge,xpb,ypb,zpb,int13
     &        ,int13p,int13_x,fsbend_slt,fsbend_slv,fxx,fyy,fzz,phiw,m1
     &        ,tag_bndg,erf_corr,erf_arr_corr,delew,rlew)
         coul_bnd_slt=coul_bnd_slt+fsbend_slt
         coul_bnd_slv=coul_bnd_slv+fsbend_slv
         fsbend=fsbend_slt
      END IF

*=======================================================================
*----- Subtract "intramolecular term" in the ZERO cell:   1-4 Fudged ---
*=======================================================================

      fsin14=0.0D0
      IF(int14p .NE. 0) THEN
         CALL ferrf(ss_index,alphal,chrge,fudgec,xpb,ypb,zpb,1,int14
     &        ,int14p,int14_x,fsin14_slt,fsin14_slv,fxx(1,2),fyy(1,2)
     &        ,fzz(1,2),phiw,erf_corr,erf_arr_corr,delew,rlew)
         fsin14=fsin14_slt
         coul_bnd_slv=coul_bnd_slv+fsin14_slv
         coul_bnd_slt=coul_bnd_slt+fsin14_slt
      END IF

*=======================================================================
*---- Compute stress tensor if coupling is by group --------------------
*---- Must be changed for parallel machines !!!!                     ---
*=======================================================================

      IF(pressure .OR. cpress .AND. coupl_grp) THEN
         CALL comp_stress_bnd(nstart,nend,atomp,virial,co,xpcm,ypcm,zpcm
     &        ,fxx(1,2),fyy(1,2),fzz(1,2))
      END IF

#if defined PARALLEL
      CALL P_Change_Decomposition("N1-N2",ntap,fxx(1,1),fyy(1,1),fzz(1,1
     &     ),nstart,nend,nstart_2,nend_2,node,nprocs)
      CALL P_Change_Decomposition("N1-N2",ntap,fxx(1,2),fyy(1,2),fzz(1,2
     &     ),nstart,nend,nstart_2,nend_2,node,nprocs)
#ifdef __Phi
      CALL P_fold_r8(ntap,phiw,nstart,nend,nlocal,node,nprocs)
      CALL P_expand_r8(phiw,nstart,nend,nlocal,node,nprocs)
#endif
#endif      
      
#ifdef __Phi
      phi(nstart_2:nend_2)=phi(nstart_2:nend_2)+phiw(nstart_2:nend_2)
#endif

*=======================================================================
*---- Add up all contributions from tagged ferrf and then to         ---
*---- full force                                                     ---
*=======================================================================

      CALL _AXPY_(n_loc_2,1.0D0,fxx(nstart_2,2),one,fpx(nstart_2),one)
      CALL _AXPY_(n_loc_2,1.0D0,fyy(nstart_2,2),one,fpy(nstart_2),one)
      CALL _AXPY_(n_loc_2,1.0D0,fzz(nstart_2,2),one,fpz(nstart_2),one)
      CALL _AXPY_(n_loc_2,1.0D0,fxx(nstart_2,1),one,fpx(nstart_2),one)
      CALL _AXPY_(n_loc_2,1.0D0,fyy(nstart_2,1),one,fpy(nstart_2),one)
      CALL _AXPY_(n_loc_2,1.0D0,fzz(nstart_2,1),one,fpz(nstart_2),one)


*=======================================================================
*---- Merge energies and expands forces in parallel runs             ---
*=======================================================================

#if defined PARALLEL
      IF(nprocs .GT. 1) THEN
         CALL P_merge_r8(coul_bnd_slt)
         CALL P_merge_r8(coul_bnd_slv)
         CALL P_merge_r8(fsbond)
         CALL P_merge_r8(fsbend)
         CALL P_merge_r8(fsin14)
         CALL P_merge_vecr8(virial,9)
      END IF
#endif         

*=======================================================================
*---- Finally add contributions from constraints                     ---
*=======================================================================

      fsbond=fsbond+fscnstr_slt
      coul_bnd_slt=coul_bnd_slt+fscnstr_slt
      coul_bnd_slv=coul_bnd_slv+fscnstr_slv
      CALL timer(dd,time_b,ela)
      ela=(time_b-time_a)
#if defined PARALLEL
      IF(nprocs .GT. 1) CALL P_merge_r8(ela)
#endif
      time_avg=time_avg+ela/DFLOAT(nprocs)
      ntime_avg=ntime_avg+1

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
