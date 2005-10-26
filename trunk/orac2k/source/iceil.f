      INTEGER FUNCTION iceil(j,k)

************************************************************************
*   Time-stamp: <1999-10-29 17:49:52 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Fri Oct 29 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER j,k

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i
      REAL*8  a,b

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      a=DBLE(j)/DBLE(k)
      b=DINT(a)
      iceil=IDINT(a)
      IF(DABS(b-a) .GT. 1.0D-5) THEN
         iceil=iceil+1
      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
