      SUBROUTINE fft_back(nprocs,array,descQ,fftable,ffwork,nfft1,nfft2
     &     ,nfft3,nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork)

************************************************************************
*   Time-stamp: <2005-01-29 11:48:40 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
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
      INTEGER nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,descQ(*)
     &     ,nfftable,nffwork,nprocs

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER isign,isys(4),one
      REAL*8  scale
      INTEGER isignb
      DATA one/1/

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      isign = -1

      CALL pubz3d(isign,nfft1,nfft2,nfft3,array,
     $   nfftdim1,nfftdim2,fftable,nfftable,ffwork,nffwork)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END

