#if defined _CRAY_ | defined T3E
#define _DOT_     sdot
#else
#define _DOT_  ddot
#endif
      subroutine fft_pme_dipole(node,nodex,nodey,nodez,ictxt,npy,npz
     &     ,descQ,nprocs,ncube,nbyte,rbyte,nstart_2,nend_2,nlocal_2,ntap
     &     ,xp,yp,zp,xpcm,ypcm,zpcm,pmechg,dipole,co,oc,volume,alphal
     &     ,order,nfft1,nfft2,nfft3,nfft3_start,nfft3_local,nfft2_start
     &     ,nfft2_local,eer,fpx,fpy,fpz,phi,dphix,dphiy
     &     ,dphiz,rkcut)

************************************************************************
*   Time-stamp: <1999-10-12 14:00:30 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Piero Procacci LSM U. Firenze                  *
*              Modified by: Massimo Marchi                             *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Feb  9 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER   NFFT1,NFFT2,NFFT3,NUMATOMS,ORDER,node,nprocs,ncube
     &     ,nbyte,rbyte,nodex,nodey,nodez,ictxt,npy,npz,descQ(*)
     &     ,nfft3_start,nfft3_local,nfft2_start,nfft2_local
      INTEGER  ntap,count,mia,ia
      INTEGER nstart_2,nend_2,nlocal_2
      REAL*8   dphix(*),dphiy(*),dphiz(*),phi(*)
      REAL*8   xp(*),yp(*),zp(*),xpcm(*),ypcm(*),zpcm(*),pmechg(*)
     &     ,dipole(3,*),fpx(*),fpy(*),fpz(*)
      REAL*8   co(3,3),oc(3,3),volume,alphal,stressc(3,3),eer
      LOGICAL pressure

*----------------------- DYNAMIC ALLOCATION ---------------------------*

      REAL(8), DIMENSION (:), ALLOCATABLE :: Q,x,y,z,fx,fy,fz,ffwork
      REAL(8), DIMENSION (:,:), ALLOCATABLE :: theta1,dtheta1,d2theta1
     &     ,theta2,dtheta2,d2theta2,theta3,dtheta3,d2theta3
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: indk1,indk2,indj1,indj2
      INTEGER, DIMENSION (:),   ALLOCATABLE :: mk,mj

*----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'pme.h'

      REAL*8  rkcut
      REAL*8 recip(3,3),vir(3,3)
      REAL*8 sx,sy,sz,stc1,stc2,stc3,stc4,stc5,stc6,stc7,stc8,stc9,du,t1
     &     ,t0
#if defined T3E | defined _CRAY_
      INTEGER*8 n_loc,one
#else
      INTEGER n_loc,one
#endif
      REAL*8 _DOT_
      DATA one/1/

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,k,nax,nay,naz,lena,dim_Q,mx,ierr1,M_get_length,len
     &     ,mlimit(3)
      REAL*8  gcpu_ll,vfcp_ll,tfcp_ll,tdelta_ll,elapse
#if defined T3E
      INTEGER*8 iceil,ihplen,npya,npza,nfft2a,nfft3a
#else
      INTEGER   iceil,ihplen,npya,npza,nfft2a,nfft3a
#endif

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      n_loc=nlocal_2
      numatoms=ntap
      dim_Q=siz_Q
      nax=nfft1+1
      nay=nfft2+1
      naz=nfft3+1

      nax=(nfft1/2+1)*2
      nfft2a=nfft2
      nfft3a=nfft3
      npya=npy
      npza=npz
      nay=nfft2a
      naz=nfft3_local
      len=nax*nay*naz
      dim_Q=len

      ALLOCATE(Q(len))

#if defined _GPFFT_
      len=numatoms
      ALLOCATE(ffwork(len))
#else
      ALLOCATE(ffwork(len))
#endif

      ALLOCATE(theta1(order,numatoms),theta2(order,numatoms)
     &     ,theta3(order,numatoms),dtheta1(order,numatoms),dtheta2(order
     &     ,numatoms),dtheta3(order,numatoms),d2theta1(order,numatoms)
     &     ,d2theta2(order,numatoms),d2theta3(order,numatoms))
      ALLOCATE(indk1(order,numatoms),indk2(order,numatoms)
     &     ,indj1(order,numatoms),indj2(order,numatoms),mk(order
     &     *numatoms),mj(order*numatoms))
      ALLOCATE(x(numatoms),y(numatoms),z(numatoms),fx(numatoms)
     &     ,fy(numatoms),fz(numatoms))
         
c---  transform all to cartesian coordinates
      
      CALL change_frame(co,oc,1,ntap,xp,yp,zp,x,y,z)
      CALL zeroa(dphix,dphiy,dphiz,ntap,1)
      CALL zeroa(fx,fy,fz,ntap,1)

c--   call stand-alone Darden's routines 

      DO i=1,3
         DO j=1,3
            recip(i,j) = 0.5*oc(j,i)
         END DO
      END DO

      CALL timer(du,t0,du)
      CALL do_pmesh_kspace_dipole(node,nodex,nodey,nodez,ictxt,npy,npz
     &     ,descQ,nprocs,ncube,nbyte,rbyte,numatoms,x,y,z,pmechg,dipole
     &     ,recip,volume,alphal,order,nfft1,nfft2,nfft3,nfft3_start
     &     ,nfft3_local,nfft2_start,nfft2_local,nax,nay,naz,eer,fx,fy,fz
     &     ,phi,dphix,dphiy,dphiz,vir,sizfftab,sizffwrk,siztheta,dim_Q
     &     ,bsp_mod1,bsp_mod2,bsp_mod3,fftable,Q,isys,ffwork,theta1
     &     ,theta2,theta3,dtheta1,dtheta2,dtheta3,d2theta1,d2theta2
     &     ,d2theta3,indk1,indk2,indj1,indj2,mk,mj,rkcut)

#ifdef PARALLEL
      IF(nprocs .GT. 1) THEN
         CALL P_merge_r8(eer,node,nprocs,ncube,rbyte)
      END IF
#endif

      DO i=1,ntap
         fpx(i) = fpx(i)+fx(i)
         fpy(i) = fpy(i)+fy(i)
         fpz(i) = fpz(i)+fz(i)
      END DO

      DEALLOCATE(Q,ffwork,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3
     &     ,d2theta1,d2theta2,d2theta3,indk1,indk2,indj1,indj2,mk,mj,x,y
     &     ,z,fx,fy,fz)

      CALL timer(du,t1,du)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
