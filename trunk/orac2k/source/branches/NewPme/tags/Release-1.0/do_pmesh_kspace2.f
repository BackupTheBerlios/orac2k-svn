      SUBROUTINE do_pmesh_kspace2(node,nodex,nodey,nodez,ictxt,npy,npz
     &     ,descQ,nprocs,ncube,nbyte,rbyte,numatoms,fr1,fr2,fr3,charge
     &     ,recip,volume,ewald_coeff,order,nfft1,nfft2,nfft3,eer,dx,dy
     &     ,dz,ene,virial,sizfftab,sizffwrk,siztheta,siz_Q,bsp_mod1
     &     ,bsp_mod2,bsp_mod3,fftable,Q,isys,ffwork,theta1,theta2,theta3
     &     ,dtheta1,dtheta2,dtheta3,indk1,indk2,indj1,indj2,mk,mj,rkcut)

************************************************************************
*   Time-stamp: <95/01/07 00:47:53 marchi>                             *
*                                                                      *
c INPUT                                                               
c       numatoms:  number of atoms
c       x,y,z:   atomic coords
c       charge  atomic charges
c       recip: 3x3 array of reciprocal unit cell vectors (stored as columns)
c       volume: the volume of the unit cell
c       ewald_coeff:   ewald convergence parameter
c       order: the order of Bspline interpolation. E.g. cubic is order 4
c          fifth degree is order 6 etc. The order must be an even number 
c          and at least 4.
c       nfft1,nfft2,nfft3: the dimensions of the charge grid array
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Tom Darden NIH                                 *
*              Modifications: Massimo Marchi                           *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Feb  9 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER numatoms,order,nfft1,nfft2,nfft3,isys(*)
      INTEGER indk1(*),indk2(*),indj1(*),indj2(*),mk(*),mj(*)
      INTEGER node,nprocs,ncube,nbyte,rbyte,nodex,nodey,nodez,ictxt,npy
     &     ,npz,descQ(*)
      REAL*8  charge(*)
     &     ,recip(3,3),volume,ewald_coeff,rkcut
c OUTPUT
c       eer:  ewald reciprocal or k-space  energy
c       dx,dy,dz: forces incremented by k-space sum
c       virial:  virial due to k-space sum (valid for atomic scaling;
c                rigid molecule virial needs a correction term not
c                computed here
      REAL*8  ene(*)
      REAL*8  eer,dx(numatoms),dy(numatoms),dz(numatoms),virial(3,3)

c SIZES of some arrays
      INTEGER    sizfftab,sizffwrk,siztheta,siz_Q


c HEAP STORAGE:  These arrays need to be preserved throughout simulation
      REAL*8  bsp_mod1(nfft1),bsp_mod2(nfft2),bsp_mod3(nfft3)
     &     ,fftable(sizfftab)
c STACK STORAGE: These arrays can be tossed after leaving this routine

      REAL*8  Q(*),ffwork(*),theta1(*),theta2(*),theta3(*),dtheta1(*)
     &     ,dtheta2(*),dtheta3(*),fr1(*),fr2(*),fr3(*)

      INTEGER  nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,kstart,kend,jstart,jend,nay,naz,p,gp_dim3,gp_n
      REAL*8  gcpu_ll,vfcp_ll,tfcp_ll,tdelta_ll,elapse
#if defined T3E
      INTEGER*8 iceil,npya,npza,nfft2a,nfft3a
#else
      INTEGER   iceil,npya,npza,nfft2a,nfft3a
#endif

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      nfft2a=nfft2
      nfft3a=nfft3
      npya=npy
      npza=npz
      nay=nfft2
      naz=nfft3
#ifdef PARALLEL
      nay=iceil(nfft2a,npya)
      naz=iceil(nfft3a,npza)
#endif
      kstart=nodez*naz+1
      kend  =(nodez+1)*naz
      jstart=nodey*nay+1
      jend  =(nodey+1)*nay

c  get some integer array dimensions
      CALL get_fftdims(nfft1,nay,naz,nfftdim1,nfftdim2,nfftdim3
     &     ,nfftable,nffwork,sfft,sffw)
      gp_dim3=nfftdim3
#ifdef _GPFFT_
      gp_dim3 = nfft3
      gp_n = nfft3/2
      if ( nfft3 .EQ. 2*gp_n ) gp_dim3 = nfft3+1
#endif

      CALL get_scaled_fractionals(numatoms,fr1,fr2,fr3,recip,nfft1,nfft2
     &     ,nfft3,fr1,fr2,fr3)

      CALL get_bspline_coeffs(numatoms,fr1,fr2,fr3,order,theta1,theta2
     &     ,theta3,dtheta1,dtheta2,dtheta3,nfft1,nfft2,nfft3,kstart,kend
     &     ,jstart,jend,indk1,indk2,indj1,indj2,mk,mj)

      CALL fill_charge_grid(node,nodey,nodez,nay,naz,numatoms
     &     ,charge,theta1,theta2,theta3,fr1,fr2,fr3,order,nfft1,nfft2
     &     ,nfft3,nfftdim1,nfftdim2,nfftdim3,Q,indk1,indk2,indj1,indj2
     &     ,mk,mj)

      CALL fft_back(nprocs,Q,descQ,fftable,ffwork,nfft1,nfft2,nfft3
     &     ,nfftdim1,nfftdim2,gp_dim3,nfftable,nffwork)

      CALL scalar_sum(node,nodey,nodez,nay,naz,Q,ewald_coeff,volume
     &     ,recip,bsp_mod1,bsp_mod2,bsp_mod3,nfft1,nfft2,nfft3,nfftdim1
     &     ,nfftdim2,nfftdim3,eer,virial,rkcut)

      CALL fft_forward(nprocs,Q,descQ,fftable,ffwork,nfft1,nfft2,nfft3
     &     ,nfftdim1,nfftdim2,gp_dim3,nfftable,nffwork)

      CALL grad_sum2(node,nodey,nodez,nay,naz,numatoms,charge,recip
     &     ,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3,dx,dy,dz,ene
     &     ,fr1,fr2,fr3,order,nfft1,nfft2,nfft3,nfftdim1,nfftdim2
     &     ,nfftdim3,Q,indk1,indk2,indj1,indj2,mk,mj)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
