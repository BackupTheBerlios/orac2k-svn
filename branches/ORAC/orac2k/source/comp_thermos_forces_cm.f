      SUBROUTINE comp_thermos_forces_cm(fth,vco,masspr)

************************************************************************
*   Time-stamp: <99/02/18 18:42:26 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Feb 18 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      REAL*8  fth,vco(3,3),masspr

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      DO i=1,3
         DO j=i,3
            fth=fth+masspr*vco(i,j)**2
         END DO
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
