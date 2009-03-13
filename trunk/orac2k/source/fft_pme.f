#if defined _CRAY_ | defined T3E
#define _DOT_     sdot
#else
#define _DOT_  ddot
#endif
#if defined _CRAY_ & !defined T3E
#define _INT2_   8
#else
#define _INT2_   4
#endif
      subroutine fft_pme(node,nodex,nodey,nodez,ictxt,npy,npz,descQ
     &     ,nprocs,ncube,nbyte,rbyte,nstart_2,nend_2,nlocal_2,ntap,xp,yp
     &     ,zp,xpcm,ypcm,zpcm,pmechg,co,oc,volume,alphal,order,nfft1
     &     ,nfft2,nfft3,nfft3_start,nfft3_local,nfft2_start,nfft2_local
     &     ,eer,fpx,fpy,fpz,phi,stressc,atomp,grppt,pressure,rkcut)

************************************************************************
*   Time-stamp: <1999-10-12 14:00:30 marchi>                           *
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

      Use Pme_Save
      Use RFFT3D
      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER   NFFT1,NFFT2,NFFT3,nfft3_start,nfft3_local,nfft2_start
     &     ,nfft2_local,NUMATOMS,ORDER,node,nprocs,ncube,nbyte,rbyte
     &     ,nodex,nodey,nodez,ictxt,npy,npz,descQ(*)
      INTEGER  ntap,atomp(*),grppt(2,*),count,mia,ia
      INTEGER nstart_2,nend_2,nlocal_2
      REAL*8   fpx(*),fpy(*),fpz(*),phi(*)
      REAL*8   xp(*),yp(*),zp(*),xpcm(*),ypcm(*),zpcm(*),pmechg(*)
      REAL*8   co(3,3),oc(3,3),volume,alphal,stressc(3,3),eer
      LOGICAL pressure

*----------------------- DYNAMIC ALLOCATION ---------------------------*

      REAL(8), DIMENSION (:), ALLOCATABLE :: x,y,z,fx,fy,fz,ffwork
      REAL(8), DIMENSION (:,:), ALLOCATABLE :: theta1,dtheta1,theta2
     &     ,dtheta2,theta3,dtheta3,d2theta1,d2theta2,d2theta3

      INTEGER, DIMENSION (:,:), ALLOCATABLE :: indk1,indk2,indj1,indj2
      INTEGER, DIMENSION (:),   ALLOCATABLE :: mk,mj
      Real(8), Dimension(:,:,:), Allocatable :: Q

*----------------------- VARIABLES IN COMMON --------------------------*

#if defined PARALLEL
      INCLUDE 'mpif.h'
#endif
      REAL*8  rkcut
      REAL*8 recip(3,3),vir(3,3)
      REAL*8 sx,sy,sz,stc1,stc2,stc3,stc4,stc5,stc6,stc7,stc8,stc9
      Integer :: k1,k2,k3
      INTEGER n_loc,one
      INTEGER ier
      REAL*8 _DOT_
      DATA one/1/

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,k,lena,dim_Q,mx,ierr1,M_get_length,len
      REAL*8  gcpu_ll,vfcp_ll,tfcp_ll,tdelta_ll,elapse,pi
      INTEGER   iceil,ihplen,npya,npza,nfft2a,nfft3a
      LOGICAL already_allocated
      DATA already_allocated/.FALSE./
      Integer, save  :: MyFirst=0,nax,nay,naz
      

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      nax=naax
      nay=naay
      naz=naaz
      pi=4.0D0*DATAN(1.0D0)
      n_loc=nlocal_2
      numatoms=ntap
      dim_Q=siz_Q

      nfft2a=nfft2
      nfft3a=nfft3
      npya=npy
      npza=npz

      len=nax*nay*naz

      dim_Q=len
      
      ALLOCATE(Q(nax,nay,nax))

#if defined _FFT_CRAY_
      mx=max(nfft1,nfft2)
      mx=max(nfft3,mx)+1
      len=512*mx
      ALLOCATE(ffwork(len))
#elif defined _GPFFT_
      len=numatoms
      ALLOCATE(ffwork(len))
#else
      len=1
      ALLOCATE(ffwork(len))
#endif
      ALLOCATE(theta1(order,numatoms)
     &     ,theta2(order,numatoms),theta3(order,numatoms)
     &     ,dtheta1(order,numatoms),dtheta2(order,numatoms)
     &     ,dtheta3(order,numatoms))
      ALLOCATE(indk1(order,numatoms),indk2(order,numatoms)
     &     ,indj1(order,numatoms),indj2(order,numatoms),mk(order
     &     *numatoms),mj(order*numatoms))
      ALLOCATE(x(numatoms),y(numatoms),z(numatoms),fx(numatoms)
     &     ,fy(numatoms),fz(numatoms))
         
c---  transform all to cartesian coordinates
      
      CALL change_frame(co,oc,1,ntap,xp,yp,zp,x,y,z)
      CALL zeroa(fx,fy,fz,ntap,1)

c--   call stand-alone Darden's routines 

      DO i=1,3
         DO j=1,3
            recip(i,j) = 0.5*oc(j,i)
         END DO
      END DO

      CALL do_pmesh_kspace(node,nodex,nodey,nodez,ictxt,npy,npz,descQ
     &     ,nprocs,ncube,nbyte,rbyte,numatoms,x,y,z,pmechg,recip,volume
     &     ,alphal,order,nfft1,nfft2,nfft3,nfft3_start,nfft3_local
     &     ,nfft2_start,nfft2_local,nax,nay,naz,eer,fx,fy,fz,phi,vir
     &     ,sizfftab,sizffwrk,siztheta,dim_Q,bsp_mod1,bsp_mod2,bsp_mod3
     &     ,fftable,Q,isys,ffwork,theta1,theta2,theta3,dtheta1,dtheta2
     &     ,dtheta3,indk1,indk2,indj1,indj2,mk,mj,rkcut)

#ifdef PARALLEL
      IF(nprocs .GT. 1) THEN
         CALL P_merge_r8(eer,node,nprocs,ncube,rbyte)
      END IF
#endif
      IF(pressure) THEN
         DO i=1,ntap
            j=atomp(i)
            sx=xp(i)-xpcm(j)
            sy=yp(i)-ypcm(j)
            sz=zp(i)-zpcm(j)
            x(i)=sx*co(1,1)+sy*co(1,2)+sz*co(1,3)
            y(i)=sx*co(2,1)+sy*co(2,2)+sz*co(2,3)
            z(i)=sx*co(3,1)+sy*co(3,2)+sz*co(3,3)
         END DO

         stc1=-_DOT_(ntap,x,one,fx,one)
         stc2=-_DOT_(ntap,y,one,fx,one)
         stc3=-_DOT_(ntap,z,one,fx,one)
         stc4=-_DOT_(ntap,x,one,fy,one)
         stc5=-_DOT_(ntap,y,one,fy,one)
         stc6=-_DOT_(ntap,z,one,fy,one)
         stc7=-_DOT_(ntap,x,one,fz,one)
         stc8=-_DOT_(ntap,y,one,fz,one)
         stc9=-_DOT_(ntap,z,one,fz,one)

         stressc(1,1)=stressc(1,1)-vir(1,1)+stc1
         stressc(1,2)=stressc(1,2)-vir(1,2)+stc2
         stressc(1,3)=stressc(1,3)-vir(1,3)+stc3
         stressc(2,1)=stressc(2,1)-vir(2,1)+stc4
         stressc(2,2)=stressc(2,2)-vir(2,2)+stc5
         stressc(2,3)=stressc(2,3)-vir(2,3)+stc6
         stressc(3,1)=stressc(3,1)-vir(3,1)+stc7
         stressc(3,2)=stressc(3,2)-vir(3,2)+stc8
         stressc(3,3)=stressc(3,3)-vir(3,3)+stc9
      END IF

      DO i=1,ntap
         fpx(i)=fpx(i)+fx(i)
         fpy(i)=fpy(i)+fy(i)
         fpz(i)=fpz(i)+fz(i)
      END DO

      DEALLOCATE(Q,ffwork,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3
     &     ,indk1,indk2,indj1,indj2,mk,mj,x,y,z,fx,fy,fz)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
