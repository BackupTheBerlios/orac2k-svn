      SUBROUTINE correc3x3(vp,fp,mass,dt)

************************************************************************
*   Time-stamp: <05/02/28 10:59:01 gmarchet>                             *
*                                                                      *
*                                                                      *
*     Advances velocities of co matrix for half time step and          *
*     correct the same velocities at full time step                    *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Mar  9 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*
      USE Module_Stress
      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8       fp(3,3)
      REAL*8       vp(3,3)
      REAL*8       dt,mass(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER      i,j
      REAL*8       tfact

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO i=1,3
         tfact=0.5*dt/mass(i)
         DO j=i,3
            IF(Orthogonal_Stress) THEN
               IF(i /= j) THEN
                  tfact=0.0D0
               END IF
            END IF
            vp(i,j)=vp(i,j) + fp(i,j)*tfact
         END DO
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
