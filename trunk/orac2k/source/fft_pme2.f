#if defined _CRAY_ | defined T3E
#define _DOT_     sdot
#else
#define _DOT_  ddot
#endif
#if defined _CRAY_ & !defined T3E
#define _INT2_   8
#elif defined T3E
#define _INT2_   4
#else
#define _INT2_   2
#endif
      subroutine fft_pme2(node,nodex,nodey,nodez,ictxt,npy,npz,descQ
     &     ,nprocs,ncube,nbyte,rbyte,nstart_2,nend_2,nlocal_2,ntap,xp,yp
     &     ,zp,xpcm,ypcm,zpcm,pmechg,co,oc,volume,alphal,order,nfft1
     &     ,nfft2,nfft3,eer,fx,fy,fz,ene,atomp
     &     ,grppt,rkcut)

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

      Use Pme_Save
      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER   NFFT1,NFFT2,NFFT3,NUMATOMS,ORDER,node,nprocs,ncube
     &     ,nbyte,rbyte,nodex,nodey,nodez,ictxt,npy,npz,descQ(*)
      INTEGER  ntap,atomp(*),grppt(2,*),count,mia,ia
      INTEGER nstart_2,nend_2,nlocal_2
      REAL*8   ene(*),fx(*),fy(*),fz(*)
      REAL*8   xp(*),yp(*),zp(*),xpcm(*),ypcm(*),zpcm(*),pmechg(*)
      REAL*8   co(3,3),oc(3,3),volume,alphal,eer

*----------------------- DYNAMIC ALLOCATION ---------------------------*

#if defined DYNAMIC_MEM
      REAL*8  Q(*),ffwork(*),theta1(order,*),dtheta1(order,*)
     &     ,theta2(order,*),dtheta2(order,*),theta3(order,*)
     &     ,dtheta3(order,*),x(*),y(*),z(*)
      INTEGER indk1(order,*),indk2(order,*),indj1(order,*),indj2(order
     &     ,*),mk(*),mj(*)
      POINTER (ip_Q,Q),(ip_ffwork,ffwork),(ip_theta1,theta1),(ip_theta2
     &     ,theta2),(ip_theta3,theta3),(ip_dtheta1,dtheta1),(ip_dtheta2
     &     ,dtheta2),(ip_dtheta3,dtheta3),(ip_indk1,indk1),(ip_indk2
     &     ,indk2),(ip_indj1,indj1),(ip_indj2,indj2),(ip_mk,mk),(ip_mj
     &     ,mj),(ip_x,x),(ip_y,y),(ip_z,z)
#endif

*----------------------- VARIABLES IN COMMON --------------------------*

      REAL*8  rkcut
      REAL*8 recip(3,3),vir(3,3)
      REAL*8 sx,sy,sz,stc1,stc2,stc3,stc4,stc5,stc6,stc7,stc8,stc9
#if !defined DYNAMIC_MEM
      INTEGER indk1(mth),indk2(mth),indj1(mth),indj2(mth),mk(max_atm)
     &     ,mj(max_atm)
      REAL*8 x(max_atm),y(max_atm),z(max_atm)
      REAL*8 Q(maxt),ffwork(maxtw),theta1(mth),dtheta1(mth),theta2(mth)
     &     ,dtheta2(mth),theta3(mth),dtheta3(mth)
      COMMON /rag1/ Q,ffwork,theta1,dtheta1,theta2,dtheta2,theta3
     &     ,dtheta3,x,y,z,indk1,indk2,indj1,indj2,mk,mj
#endif
#if defined T3E | defined _CRAY_
      INTEGER*8 n_loc,one
#else
      INTEGER n_loc,one
#endif
      REAL*8 _DOT_
      DATA one/1/

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,k,nax,nay,naz,lena,dim_Q,mx,ierr1,M_get_length,len
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

#if defined _FFT_CRAY_
      nax=nfft1+1
      nay=iceil(nfft2,npy)+1
      naz=iceil(nfft3,npz)+1
#elif defined _FFT_T3E_
      nax=nfft1+1
      nfft2a=nfft2
      nfft3a=nfft3
      npya=npy
      npza=npz
      nay=iceil(nfft2a,npya)+1
      naz=iceil(nfft3a,npza)+1
#elif defined _GPFFT_
      nax=nfft1+1
      nfft2a=nfft2
      nfft3a=nfft3
      npya=npy
      npza=npz
      nay=iceil(nfft2a,npya)+1
      naz=iceil(nfft3a,npza)+1
#endif

      len=nax*nay*naz*2
      dim_Q=len
      
#if defined DYNAMIC_MEM
      CALL memory(ip_q,len,8)
#if defined _FFT_CRAY_
      mx=max(nfft1,nfft2)
      mx=max(nfft3,mx)+1
      len=512*mx
      CALL memory(ip_ffwork,len,8)
#else
      CALL memory(ip_ffwork,len,8)
#endif
      len=order*numatoms
      CALL memory(ip_theta1,len,8)
      CALL memory(ip_theta2,len,8)
      CALL memory(ip_theta3,len,8)
      CALL memory(ip_dtheta1,len,8)
      CALL memory(ip_dtheta2,len,8)
      CALL memory(ip_dtheta3,len,8)
      
      len=order*numatoms
      
      CALL memory(ip_indk1,len,_INT2_)
      CALL memory(ip_indk2,len,_INT2_)
      CALL memory(ip_indj1,len,_INT2_)
      CALL memory(ip_indj2,len,_INT2_)
      CALL memory(ip_mk,len,_INT2_)
      CALL memory(ip_mj,len,_INT2_)
      
      len=numatoms
      CALL memory(ip_x,len,8)
      CALL memory(ip_y,len,8)
      CALL memory(ip_z,len,8)
#endif
         
c---  transform all to cartesian coordinates
      
      CALL change_frame(co,oc,1,ntap,xp,yp,zp,x,y,z)
      CALL zeroa(fx,fy,fz,ntap,1)
      CALL zero0(ene,ntap)

c--   call stand-alone Darden's routines 

      DO i=1,3
         DO j=1,3
            recip(i,j) = 0.5*oc(j,i)
         END DO
      END DO

      CALL do_pmesh_kspace2(node,nodex,nodey,nodez,ictxt,npy,npz,descQ
     &     ,nprocs,ncube,nbyte,rbyte,numatoms,x,y,z,pmechg,recip,volume
     &     ,alphal,order,nfft1,nfft2,nfft3,eer,fx,fy,fz,ene,vir,sizfftab
     &     ,sizffwrk,siztheta,dim_Q,bsp_mod1,bsp_mod2,bsp_mod3,fftable,Q
     &     ,isys,ffwork,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3
     &     ,indk1,indk2,indj1,indj2,mk,mj,rkcut)

#ifdef PARALLEL
      IF(nprocs .GT. 1) THEN
         CALL P_fold_r8(ntap,fx,nstart_2,nend_2,nlocal_2,node,nprocs)
         CALL P_fold_r8(ntap,fy,nstart_2,nend_2,nlocal_2,node,nprocs)
         CALL P_fold_r8(ntap,fz,nstart_2,nend_2,nlocal_2,node,nprocs)
         CALL P_fold_r8(ntap,ene,nstart_2,nend_2,nlocal_2,node,nprocs)

         CALL P_merge_r8(eer,node,nprocs,ncube,rbyte)

         CALL P_expand_r8(fx,nstart_2,nend_2,nlocal_2,node,nprocs)
         CALL P_expand_r8(fy,nstart_2,nend_2,nlocal_2,node,nprocs)
         CALL P_expand_r8(fz,nstart_2,nend_2,nlocal_2,node,nprocs)
         CALL P_expand_r8(ene,nstart_2,nend_2,nlocal_2,node,nprocs)
      END IF
#endif

#if defined DYNAMIC_MEM
      CALL free_mem(ip_q)
      CALL free_mem(ip_ffwork)
      CALL free_mem(ip_theta1)
      CALL free_mem(ip_theta2)
      CALL free_mem(ip_theta3)
      CALL free_mem(ip_dtheta1)
      CALL free_mem(ip_dtheta2)
      CALL free_mem(ip_dtheta3)
      CALL free_mem(ip_indk1)
      CALL free_mem(ip_indk2)
      CALL free_mem(ip_indj1)
      CALL free_mem(ip_indj2)
      CALL free_mem(ip_mk)
      CALL free_mem(ip_mj)
      CALL free_mem(ip_x)
      CALL free_mem(ip_y)
      CALL free_mem(ip_z)
#endif

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
