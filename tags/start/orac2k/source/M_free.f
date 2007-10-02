      SUBROUTINE M_free(ip_point)

************************************************************************
*   Time-stamp: <1999-11-22 11:52:48 marchi>                             *
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

#if defined _CRAY_ | defined T3E | defined AXP
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
      CALL hpdeallc(ip_point,ierr1,one)
#elif defined AXP  | defined OSF1
      CALL free(ip_point)
#elif AIX
      CALL free(%val(ip_point))
#endif

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
