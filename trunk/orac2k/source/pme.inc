#include "config.h"	
      INTEGER mfft1,mfft2,mfft3,maxord,maxt,mth,max_atm,maxn,maxtw
      PARAMETER(mfft1=_FFT1_+1, mfft2=_FFT2_+1, mfft3=_FFT3_+1)
      PARAMETER(maxn=_MAX_GRID_+1,maxord=_MORD_,
     &     maxt=2*mfft1*mfft2*mfft3,max_atm=_SIT_SOLU_,mth=maxord
     &     *max_atm)
      REAL*8  bsp_mod1(mfft1),bsp_mod2(mfft2),bsp_mod3(mfft3)
#if  defined _FFT_CRAY_
      PARAMETER (maxtw=512*maxn)
      REAL*8    fftable(2*(mfft1+mfft2+mfft3+50))
      INTEGER isys(4)
#elif defined _FFT_T3E_ &  !defined PARALLEL
      PARAMETER (maxtw=512*maxn)
      REAL*8  fftable(2*(mfft1+mfft2+mfft3+50))
      INTEGER isys(4)
#elif defined _FFT_T3E_ &  defined PARALLEL
      PARAMETER (maxtw=maxt)
      REAL*8  fftable(2*12*(mfft1+mfft2+mfft3))
      INTEGER isys(4)
#elif  defined _GPFFT_
      PARAMETER (maxtw=2*maxn)
      REAL*8  fftable(4*1024)
      INTEGER isys(4)
#else
      PARAMETER (maxtw=2*maxn)
      REAL*8  fftable(3*(maxn*4+15))
      INTEGER isys(4)
#endif
      INTEGER sizfftab,sizffwrk,siztheta,siz_Q,sizheap,sizstack
      COMMON /pme_data/ bsp_mod1,bsp_mod2,bsp_mod3,fftable,sizfftab
     &     ,sizffwrk,siztheta,siz_Q,sizheap,sizstack,isys








