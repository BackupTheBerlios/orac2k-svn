      SUBROUTINE M_free_s(ip_point)

************************************************************************
*   Time-stamp: <99/02/14 15:37:03 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed Feb 10 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

#if defined _CRAY_ | defined T3E
      INTEGER*8 ip_point
#else
      INTEGER ip_point
#endif

*------------------------- LOCAL VARIABLES ----------------------------*

#if defined _CRAY_ | defined T3E
      INTEGER*8 ierr1,one
#else
      INTEGER ierr1,one
#endif
      DATA one/1/

*----------------------- EXECUTABLE STATEMENTS ------------------------*

#if defined _CRAY_ | defined T3E
      CALL shpdeallc(ip_point,ierr1,one)
#endif

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
