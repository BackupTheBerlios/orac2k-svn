<MODULE FOURIER_Mod

!!$***********************************************************************
!!$   Time-stamp: <2005-10-25 18:55:26 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Oct 20 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*

  USE INPUT_Mod, ONLY: Read_String, Parser, err_open,err_end,err_unr&
       &,err_fnf,err_args
  USE xerror_mod

  LOGICAL, save :: Electrostatics=.FALSE.,Gaussian=.TRUE.
  INTEGER, SAVE :: n_write,kfourier,numatoms
  CHARACTER(80), SAVE :: filename
  
  INTEGER, SAVE :: nfft1=0,nfft2=0,nfft3=0,nfft3_start=0,nfft3_local&
       &=0,nfft2_start=0,nfft2_local=0,order=0
  REAL(8), SAVE :: alpha,oca(3,3)
  INTEGER, SAVE :: sizfftab,sizffwrk,siztheta,siz_Q,sizheap,sizstack
  INTEGER, SAVE :: mfft1,mfft2,mfft3,maxord,maxt,mth,max_atm,maxn,maxtw
  INTEGER, SAVE :: nax,nay,naz
  
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: Ex0,Ey0,Ez0
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: Ex,Ey,Ez
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: phi,rho
  REAL(8), DIMENSION(:), ALLOCATABLE, SAVE::  bsp_mod1,bsp_mod2,bsp_mod3
  INTEGER, SAVE :: node=0,nprocs=1,nodex=0,nodey=0,nodez=0,npy=0,npz=0

CONTAINS

!!$======================== DECLARATIONS ================================*

  SUBROUTINE Initialize(A_node,A_nprocs,A_nodex,A_nodey,A_nodez,A_npy&
       &,A_npz,ntap)
    
    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

    INTEGER :: A_node,A_nprocs,A_nodex,A_nodey,A_nodez,A_npy,A_npz,ntap

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: i,ictxt,descQ(2),fftable(2),iret=0,nfft2a
    CHARACTER(80) :: errmsg=' '

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*


    numatoms=ntap
    A_node=node
    A_nprocs=nprocs
    A_nodex=nodex
    A_nodey=nodey
    A_nodez=nodez
    A_npy=npy
    A_npz=npz
    nax=(nfft1/2+1)*2
    nfft2a=nfft2
    npya=npy
    npza=npz
    nay=nfft2a
    naz=nfft3_local

    CALL load_bsp_moduli(bsp_mod1,bsp_mod2,bsp_mod3,nfft1,nfft2,nfft3&
         &,order)

    CALL Pme_init(node,nprocs,nodex,nodey,nodez,npy,npz,ictxt,descQ&
         &,fftable,nfft1,nfft2,nfft3,nfft3_start,nfft3_local&
         &,nfft2_start,nfft2_local,iret,errmsg)

    IF(iret /= 0) CALL abort_now(errmsg)

  END SUBROUTINE Initialize
  SUBROUTINE Compute(xp0,yp0,zp0,co,oc)
    IMPLICIT NONE 

!!$----------------------------- ARGUMENTS ------------------------------*

    INTEGER :: xp0(*),yp0(*),zp0(*),co(3,3),oc(3,3)

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: naay,naaz
    REAL(8), DIMENSION (:), ALLOCATABLE :: Q,x,y,z,fx,fy,fz,ffwork
    REAL(8), DIMENSION (:,:), ALLOCATABLE :: theta1,dtheta1,theta2&
         &,dtheta2,theta3,dtheta3,d2theta1,d2theta2,d2theta3
    INTEGER, DIMENSION (:,:), ALLOCATABLE :: indk1,indk2,indj1,indj2
    INTEGER, DIMENSION (:),   ALLOCATABLE :: mk,mj


!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    ALLOCATE(phi(nax,nay,naz))
    ALLOCATE(theta1(order,numatoms),theta2(order,numatoms)&
         &,theta3(order,numatoms),dtheta1(order,numatoms)&
         &,dtheta2(order,numatoms),dtheta3(order,numatoms))
    ALLOCATE(indk1(order,numatoms),indk2(order,numatoms)&
         &,indj1(order,numatoms),indj2(order,numatoms),mk(order&
         &*numatoms),mj(order*numatoms))
    ALLOCATE(x(numatoms),y(numatoms),z(numatoms),fx(numatoms)&
         &,fy(numatoms),fz(numatoms))

    fx=0.0D0
    fy=0.0D0
    fz=0.0D0
    oca = 0.5*oc

    naay=nfft2
    naaz=nfft3
    naaz=nfft3_local
    
    kstart=nfft3_start
    kend  =nfft3_start+nfft3_local-1
    jstart=nodey*naay+1
    jend  =(nodey+1)*naay

    CALL get_scaled_fractionals(numatoms,x,y,z,recip,nfft1&
         &,nfft2,nfft3,x,y,z)

    CALL get_bspline_coeffs(numatoms,x,y,z,order,theta1,theta2&
         &,theta3,dtheta1,dtheta2,dtheta3,nfft1,nfft2,nfft3,kstart&
         &,kend,jstart,jend,indk1,indk2,indj1,indj2,mk,mj)

    CALL fill_charge_grid(node,nodey,nodez,jstart,kstart,numatoms&
         &,charge,theta1,theta2,theta3,x,y,z,order,nfft1,nfft2&
         &,nfft3,nax,nay,naz,phi,indk1,indk2,indj1,indj2,mk,mj)


!!$=======================================================================
!!$--- Use FFTW_TRANSPOSED_ORDER
!!$========================================================================

      CALL do_rfft3d(idir,Q)
      CALL timer(vfcp,time3,elapse)
      Time4=time3-time2

#if defined PARALLEL 

!!$=======================================================================
!!$--- Use transposed Q matrix with FFTW_TRANSPOSED_ORDER
!!$========================================================================

      nc1=nd1
      nc2=nfft3
      nc3=nfft2_local
      nb1=nfft1
      nb2=nfft3
      nb3=nfft2
      nb3_start=nfft2_start
      nb3_local=nfft2_local
      CALL scalar_sum_transp(node,nb3_start,nb3_local,Q,ewald_coeff
     &     ,volume,recip,bsp_mod1,bsp_mod2,bsp_mod3,nb1,nb2,nb3,nc1,nc2
     &     ,nc3,eer,virial,rkcut)
#else
      nc1=nd1
      nc2=nfft2
      nc3=nfft3_local
      nb1=nfft1
      nb2=nfft2
      nb3=nfft3
      nb3_start=nfft3_start
      nb3_local=nfft3_local
      CALL scalar_sum_normal(node,nb3_start,nb3_local,Q,ewald_coeff
     &     ,volume,recip,bsp_mod1,bsp_mod2,bsp_mod3,nb1,nb2,nb3,nc1,nc2
     &     ,nc3,eer,virial,rkcut)      
#endif
      idir=-1
      CALL timer(vfcp,time2,elapse)

!!$=======================================================================
!!$--- Use FFTW_TRANSPOSED_ORDER
!!$========================================================================

      CALL do_rfft3d(idir,Q)


    


  END SUBROUTINE Compute

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

  INCLUDE 'FOURIER_Read.f90'
END MODULE FOURIER_Mod
