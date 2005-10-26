      SUBROUTINE do_pmesh_kspace(node,nodex,nodey,nodez,ictxt,npy,npz
     &     ,descQ,nprocs,ncube,nbyte,rbyte,numatoms,fr1,fr2,fr3,charge
     &     ,recip,volume,ewald_coeff,order,nfft1,nfft2,nfft3,nfft3_start
     &     ,nfft3_local,nfft2_start,nfft2_local,nd1,nd2,nd3,eer,dx,dy,dz
     &     ,phi,virial,sizfftab,sizffwrk,siztheta,siz_Q,bsp_mod1
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
      USE RFFT3D
      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER numatoms,order,nfft1,nfft2,nfft3,nfft3_start,nfft3_local
     &     ,nfft2_start,nfft2_local,isys(*)
      INTEGER indk1(*),indk2(*),indj1(*),indj2(*),mk(*),mj(*)
      INTEGER node,nprocs,ncube,nbyte,rbyte,nodex,nodey,nodez,ictxt,npy
     &     ,npz,descQ(*),nd1,nd2,nd3
      REAL*8  charge(*)
     &     ,recip(3,3),volume,ewald_coeff,rkcut
c OUTPUT
c       eer:  ewald reciprocal or k-space  energy
c       dx,dy,dz: forces incremented by k-space sum
c       virial:  virial due to k-space sum (valid for atomic scaling;
c                rigid molecule virial needs a correction term not
c                computed here
      REAL*8  eer,dx(numatoms),dy(numatoms),dz(numatoms),virial(3,3)

c SIZES of some arrays
      INTEGER    sizfftab,sizffwrk,siztheta,siz_Q


c HEAP STORAGE:  These arrays need to be preserved throughout simulation
      REAL*8  bsp_mod1(nfft1),bsp_mod2(nfft2),bsp_mod3(nfft3)
     &     ,fftable(sizfftab)
c STACK STORAGE: These arrays can be tossed after leaving this routine

      REAL*8  Q(nd1,nd2,nd3),ffwork(*),theta1(*),theta2(*),theta3(*)
     &     ,dtheta1(*),dtheta2(*),dtheta3(*),fr1(*),fr2(*),fr3(*),phi(*)

      INTEGER  nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw
#if defined PARALLEL
      INCLUDE 'mpif.h'
#endif

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,kstart,kend,jstart,jend,nay,naz,p,gp_dim3,gp_n,ja,ka
      REAL*8  gcpu_ll,vfcp_ll,tfcp_ll,tdelta_ll,elapse,dum1,dum2,time_a
     &     ,time_b
#if defined T3E
      INTEGER*8 iceil,npya,npza,nfft2a,nfft3a
#else
      INTEGER   iceil,npya,npza,nfft2a,nfft3a
#endif
      REAL*8  time1_avg,time2_avg,time_avg,time,time1,time2,time3,time4
     &     ,time5,time_oth,time_fft,vfcp
      INTEGER ntime_avg,idir,nb1,nb2,nb3,nc1,nc2,nc3,nb3_local,nb3_start
      DATA ntime_avg/0/
      DATA time_avg,time_fft,time_oth/0.0D0,0.0D0,0.0D0/

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      nay=nfft2
      naz=nfft3
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

      CALL timer(vfcp,time1,elapse)
      
      CALL get_scaled_fractionals(numatoms,fr1,fr2,fr3,recip,nfft1,nfft2
     &     ,nfft3,fr1,fr2,fr3)

      CALL get_bspline_coeffs(numatoms,fr1,fr2,fr3,order,theta1,theta2
     &     ,theta3,dtheta1,dtheta2,dtheta3,nfft1,nfft2,nfft3,kstart,kend
     &     ,jstart,jend,indk1,indk2,indj1,indj2,mk,mj)

      CALL fill_charge_grid(node,nodey,nodez,jstart,kstart,numatoms
     &     ,charge,theta1,theta2,theta3,fr1,fr2,fr3,order,nfft1,nfft2
     &     ,nfft3,nd1,nd2,nd3,Q,indk1,indk2,indj1,indj2
     &     ,mk,mj)

      idir=1
      CALL timer(vfcp,time2,elapse)

c$$$=======================================================================
c$$$--- Use FFTW_TRANSPOSED_ORDER
c$$$========================================================================

      CALL do_rfft3d(idir,Q)
      CALL timer(vfcp,time3,elapse)
      Time4=time3-time2

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
      idir=-1
      CALL timer(vfcp,time2,elapse)

c$$$=======================================================================
c$$$--- Use FFTW_TRANSPOSED_ORDER
c$$$========================================================================

      CALL do_rfft3d(idir,Q)
      CALL timer(vfcp,time3,elapse)
      Time4=Time4+time3-time2

      CALL grad_sum(node,jstart,kstart,numatoms,charge,recip
     &     ,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3,dx,dy,dz,phi
     &     ,fr1,fr2,fr3,order,nfft1,nfft2,nfft3,nd1,nd2,nd3,Q,indk1
     &     ,indk2,indj1,indj2,mk,mj)
      CALL timer(vfcp,time5,elapse)
      time_oth=time_oth+time5-time1-time4
      time_fft=time_fft+time4
      ntime_avg=ntime_avg+1
      time1=time_oth
      time2=time_fft
#if defined PARALLLEL
      CALL P_merge_r8(time1,node,nprocs,ncube,rbyte)
      CALL P_merge_r8(time2,node,nprocs,ncube,rbyte)
#endif
c$$$      WRITE(*,*) 'other =',time1/DBLE(ntime_avg*nprocs)
c$$$      WRITE(*,*) 'fft =',time2/DBLE(ntime_avg*nprocs)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
