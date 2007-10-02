      SUBROUTINE zeroi(a,n)

************************************************************************
*   Time-stamp: <95/05/23 17:48:51 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue May 23 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER a(*),n

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO i=1,n
         a(i)=0
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
