      SUBROUTINE set_tempp(masspr,vco,temppra,t)

************************************************************************
*   Time-stamp: <2009-03-12 16:21:45 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Mar  8 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      USE Module_Stress
      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8  masspr,vco(3,3),temppra,t

*----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'unit.h'

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j
      REAL*8  sig,tvel,t1,u1,u2,tto,ranf

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO i=1,3
         DO j=1,3
            vco(i,j)=0.0D0
         END DO
      END DO
      
      DO i=1,3
         DO j=i,3
            tvel=boltz*t/(unite*masspr)
            sig=DSQRT(tvel)
            u1=ranf()
            u2=ranf()
            t1=DSQRT(-2.0d0*DLOG(u1))*DCOS(2.0d0*pi*u2)
            IF(Orthogonal_Stress) THEN
               IF(i /= j) THEN
                  sig=0.0D0
                  t1=0.0D0
               END IF
            END IF
            vco(i,j)=t1*sig
         END DO
      END DO
      tto=0.0D0
      DO i=1,3
         DO j=i,3
            tto=tto+0.50D0*masspr*vco(i,j)*vco(i,j)
         END DO
      END DO
      tto=tto*efact
      temppra=(2.0D0/6.0D0)*tto/gascon

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
