
*-------------------- VARIABLES IN COMMON ------------------------------

*=======================================================================
*---  Steve Plimpton's decomposition vectors                        ----
*---   Atomic decomposition:                                        ----
*---   vecx_0: Used by the mts_force routine                        ----
*---   vecx_1: Used by the mts_intra routines                       ----
*---   vecx_2: Used in the h,l,m loops (forces from mts_forces are  ----
*---           expanded!)                                           ----
*---   vecx_3: Used by mts_furier to expand the bspline coefficients----
*---           of the PME calculation                               ----
*---   vecx_ex: Used for ballistic atoms in the n0, n1 loops        ----
*---            (no force from loop acting except constraints)      ----
*---                                                                ----
*---   Group decomposition:                                         ----
*---   vecx_cm : Used in the h,l,m loops                            ----
*---   vecx_cmi: Used in the n0,n1 loops                            ----
*---   vecx_cme: Used in the n0,n1 loops for ballistic groups       ----
*---                                                                ----
*=======================================================================

      INTEGER n_h
      PARAMETER (n_h = 128)

      INTEGER node,nodex,nodey,nodez,ictxt,npy,npz,nprocs,ncube
      INTEGER nstart_cm,nend_cm,nlocal_cm,nstart_cmi,nend_cmi,nlocal_cmi
     &     ,nstart_cme,nend_cme,nlocal_cme,nstart_cme0,nend_cme0
     &     ,nlocal_cme0,nstart_1,nend_1,nlocal_1,nstart_2,nend_2
     &     ,nlocal_2,nstart_ex0,nend_ex0,nlocal_ex0,nstart_ex,nend_ex
     &     ,nlocal_ex,ptr_ex,ptr_cme,ptr,ntot_1,ntot_cmi,ntot_ex
     &     ,ntot_cme,nstart_g1,nend_g1,nlocal_g1,nstart_gex,nend_gex
     &     ,nlocal_gex,nstart_gex0,nend_gex0,nlocal_gex0,ptr_gex
     &     ,descQ(12)

      INTEGER ncpu_h,ncpu_l,ncpu_m,nstart_h,nend_h,nstart_ah,nend_ah
     &     ,nlocal_h,nlocal_ah
      INTEGER nat_listp,nat_list(2,nores)
      REAL*8  cpu_h(n_h),cpu_l(n_h),cpu_m(n_h)
      COMMON /parallel/ cpu_h,cpu_l,cpu_m,nstart_cm,nend_cm,nlocal_cm
     &     ,nstart_cmi,nend_cmi,nlocal_cmi,nstart_cme,nend_cme
     &     ,nlocal_cme,nstart_cme0,nend_cme0,nlocal_cme0,nstart_1,nend_1
     &     ,nlocal_1,nstart_2,nend_2,nlocal_2,nstart_ex0,nend_ex0
     &     ,nlocal_ex0,nstart_ex,nend_ex,nlocal_ex,ptr_ex,ptr_cme,ptr
     &     ,ntot_1,ntot_cmi,ntot_ex,ntot_cme,nstart_g1,nend_g1,nlocal_g1
     &     ,nstart_gex,nend_gex,nlocal_gex,nstart_gex0,nend_gex0
     &     ,nlocal_gex0,ptr_gex,ncpu_h,ncpu_l,ncpu_m,nstart_h,nend_h
     &     ,nstart_ah,nend_ah,nlocal_h,nlocal_ah,nat_listp,nat_list,node
     &     ,nodex,nodey,nodez,ictxt,npy,npz,nprocs,ncube,descQ
