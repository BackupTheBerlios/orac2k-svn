      subroutine fft_pme_init(numatoms,nfft1,nfft2,nfft3,order,sizfftab
     &     ,sizffwrk,siztheta,siz_Q,sizheap,sizstack,bsp_mod1,bsp_mod2
     &     ,bsp_mod3,fftable,ffwork)

************************************************************************
*
*     To be called from MTSMD
*---  ON INPUT
*     nfft1,nfft2,nfft3: grid points in the k1,k2,k3 directions
*     order            : order of B-spline interpolation
*---  ON OUTPUT 
*     sizfftab is permanent 3d fft table storage
*     sizffwrk is temporary 3d fft work storage
*     siztheta is size of arrays theta1-3 dtheta1-3
*     sizheap is total size of permanent storage
*     sizstack is total size of temporary storage
*     bsp_mod1-3 hold the moduli of the inverse DFT of the B splines
*
************************************************************************

      IMPLICIT NONE

      INTEGER   NFFT1,NFFT2,NFFT3,NUMATOMS,ORDER
      INTEGER   sizfftab,sizffwrk,siztheta,siz_Q,sizheap,sizstack
      DOUBLE PRECISION bsp_mod1(*),bsp_mod2(*),bsp_mod3(*)
      DOUBLE PRECISION fftable(*),ffwork(*)


      call pmesh_kspace_get_sizes(nfft1,nfft2,nfft3,numatoms,order
     &     ,sizfftab,sizffwrk,siztheta,siz_Q,sizheap,sizstack)
      call pmesh_kspace_setup(
     $    bsp_mod1,bsp_mod2,bsp_mod3,fftable,ffwork,
     $    nfft1,nfft2,nfft3,order,sizfftab,sizffwrk)

      return
      end

