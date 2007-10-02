      SUBROUTINE fft_setup(array,fftable,ffwork,nfft1,nfft2,nfft3
     &     ,nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork)

************************************************************************
*   Time-stamp: <01/02/24 17:07:41 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Tom Darden NIHS                                *
*              Modified by: Massimo MARCHI                             *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Feb 13 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8  array(*),fftable(*),ffwork(*)
      INTEGER nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
      INTEGER nfftable,nffwork,isys(4)

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8  scale
#if defined T3E
      INTEGER*8 nfft1a,nfft2a,nfft3a,nfftdim1a,nfftdim2a,nfftdim3a
     &     ,isysa(4),isign
#else
      INTEGER isign
#endif

*----------------------- EXECUTABLE STATEMENTS ------------------------*

#ifdef SGIFFT
      call ZFFT3DI(nfft1,nfft2,nfft3,fftable)
#endif
#if defined _FFT_CRAY_
      write(6,*)'using cray fft code'
      isign = 0
      scale = 1.d0
      isys(1)=3
      isys(2)=0
      isys(3)=0
      isys(4)=0
      CALL ccfft3d(isign,nfft1,nfft2,nfft3,scale,array,nfftdim1,nfftdim2
     &     ,array,nfftdim1,nfftdim2,fftable,ffwork,isys)
#elif defined  _FFT_T3E_ & MODE != PARALLEL
      write(6,*)'using cray fft code'
      isign = 0
      scale = 1.d0
      isys(1)=3
      isys(2)=0
      isys(3)=0
      isys(4)=0
      isysa(1)=3
      isysa(2)=0
      isysa(3)=0
      isysa(4)=0
      nfft1a=nfft1
      nfft2a=nfft2
      nfft3a=nfft3
      nfftdim1a=nfftdim1
      nfftdim2a=nfftdim2
      nfftdim3a=nfftdim3
      CALL ccfft3d(isign,nfft1a,nfft2a,nfft3a,scale,array,nfftdim1a
     &     ,nfftdim2a,array,nfftdim1a,nfftdim2a,fftable,ffwork,isysa)
#elif defined _PDFFT_ 
      write(6,*)'using public domain fft code'
      call pubz3di(nfft1,nfft2,nfft3,fftable,nfftable)
#endif
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
