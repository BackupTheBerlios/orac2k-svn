      SUBROUTINE free_mem(ip_point)

************************************************************************
*   Time-stamp: <1999-11-02 10:32:28 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Nov  2 1999 -                                     *
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

*----------------------- EXECUTABLE STATEMENTS ------------------------*

#if defined T3E & defined PARALLEL
      CALL M_free_s(ip_point)
#else
      CALL M_free(ip_point)
#endif

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
