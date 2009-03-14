      SUBROUTINE do_pmesh_kspace_dipole(node,nodex,nodey,nodez,ictxt,npy
     &     ,npz,descQ,nprocs,ncube,nbyte,rbyte,numatoms,fr1,fr2,fr3
     &     ,charge,dipole,recip,volume,ewald_coeff,order,nfft1,nfft2
     &     ,nfft3,nfft3_start,nfft3_local,nfft2_start,nfft2_local,nd1
     &     ,nd2,nd3,eer,fpx,fpy,fpz,phi,field_x,field_y,field_z,virial
     &     ,sizfftab,sizffwrk,siztheta,siz_Q,bsp_mod1,bsp_mod2,bsp_mod3
     &     ,fftable,Q,isys,ffwork,theta1,theta2,theta3,dtheta1,dtheta2
     &     ,dtheta3,d2theta1,d2theta2,d2theta3,indk1,indk2,indj1,indj2
     &     ,mk,mj,rkcut)

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

      USE RFFT3D
      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER numatoms,order,nfft1,nfft2,nfft3,isys(*)
      INTEGER indk1(*),indk2(*),indj1(*),indj2(*),mk(*),mj(*)
      INTEGER node,nprocs,ncube,nbyte,rbyte,nodex,nodey,nodez,ictxt,npy
     &     ,npz,descQ(*),nfft3_start,nfft3_local,nfft2_start,nfft2_local
     &     ,nd1,nd2,nd3
      REAL*8  charge(*),dipole(3,*),recip(3,3),volume,ewald_coeff,rkcut
     &     ,fpx(*),fpy(*),fpz(*)
c OUTPUT
c       eer:  ewald reciprocal or k-space  energy
c       dx,dy,dz: forces incremented by k-space sum
c       virial:  virial due to k-space sum (valid for atomic scaling;
c                rigid molecule virial needs a correction term not
c                computed here
      REAL*8  eer,field_x(numatoms),field_y(numatoms),field_z(numatoms)
     &     ,phi(numatoms),virial(3,3)

c SIZES of some arrays
      INTEGER    sizfftab,sizffwrk,siztheta,siz_Q


c HEAP STORAGE:  These arrays need to be preserved throughout simulation
      REAL*8  bsp_mod1(nfft1),bsp_mod2(nfft2),bsp_mod3(nfft3)
     &     ,fftable(sizfftab)
c STACK STORAGE: These arrays can be tossed after leaving this routine

      REAL*8  Q(nd1,nd2,nd3),ffwork(*),theta1(*),theta2(*),theta3(*)
     &     ,dtheta1(*),dtheta2(*),dtheta3(*),d2theta1(*),d2theta2(*)
     &     ,d2theta3(*),fr1(*),fr2(*),fr3(*)

      INTEGER  nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,kstart,kend,jstart,jend,nay,naz,p,gp_dim3,gp_n
      INTEGER idir,nb1,nb2,nb3,nc1,nc2,nc3,nb3_local,nb3_start
      REAL*8  gcpu_ll,vfcp_ll,tfcp_ll,tdelta_ll,elapse,eer_2
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
      nay=nfft2a
      naz=nfft3_local

      kstart=nfft3_start
      kend  =nfft3_start+nfft3_local-1
      jstart=nodey*nay+1
      jend  =(nodey+1)*nay
c  get some integer array dimensions
      nfftable=2
      nffwork=2
      sfft=2
      sffw=2
      gp_dim3=nd3

#ifdef _GPFFT_
      gp_dim3 = nfft3
      gp_n = nfft3/2
      if ( nfft3 .EQ. 2*gp_n ) gp_dim3 = nfft3+1
#endif

      CALL get_scaled_fractionals(numatoms,fr1,fr2,fr3,recip,nfft1,nfft2
     &     ,nfft3,fr1,fr2,fr3)

      CALL get_bspline_coeffs2(numatoms,fr1,fr2,fr3,order,theta1,theta2
     &     ,theta3,dtheta1,dtheta2,dtheta3,d2theta1,d2theta2,d2theta3
     &     ,nfft1,nfft2,nfft3,kstart,kend,jstart,jend,indk1,indk2,indj1
     &     ,indj2,mk,mj)

      CALL fill_dipole_grid(node,nodey,nodez,nay,naz,numatoms
     &     ,charge,dipole,recip,dtheta1,dtheta2,dtheta3,theta1,theta2
     &     ,theta3,fr1,fr2,fr3,order,nfft1,nfft2,nfft3,nd1,nd2,nd3,Q
     &     ,indk1,indk2,indj1,indj2,mk,mj)

c$$$=======================================================================
c$$$--- Use FFTW_TRANSPOSED_ORDER
c$$$========================================================================

      idir=1
      CALL do_rfft3d(idir,Q)

#if defined PARALLEL 
c$$$=======================================================================
c$$$--- Use transposed Q matrix with FFTW_TRANSPOSED_ORDER
c$$$========================================================================

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
c$$$=======================================================================
c$$$--- Use FFTW_TRANSPOSED_ORDER
c$$$========================================================================

      idir=-1
      CALL do_rfft3d(idir,Q)
      CALL grad_sum_dipole(node,nodey,nodez,nay,naz,numatoms,charge
     &     ,dipole,recip,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3
     &     ,d2theta1,d2theta2,d2theta3,fpx,fpy,fpz,eer_2,phi,field_x
     &     ,field_y,field_z,fr1,fr2,fr3,order,nfft1,nfft2,nfft3,nd1,nd2
     &     ,nd3,Q,indk1,indk2,indj1,indj2,mk,mj)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END

