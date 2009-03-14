      INTEGER FUNCTION M_get_length(n,nbyte)

************************************************************************
*   Time-stamp: <1999-10-12 17:28:20 marchi>                             *
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
      
      INTEGER n,nbyte

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER bytes_x_unit
#if defined _CRAY_ | defined T3E
      PARAMETER (bytes_x_unit=8)
#else
      PARAMETER (bytes_x_unit=1)
#endif

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      M_get_length=n*nbyte/bytes_x_unit

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
