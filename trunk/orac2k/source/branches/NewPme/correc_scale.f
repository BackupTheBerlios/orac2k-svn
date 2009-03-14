      SUBROUTINE correc_scale(noff,n,vx,vy,vz,etap,dt)

************************************************************************
*   Time-stamp: <99/02/18 15:42:54 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Apr  6 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER n,noff
      REAL*8  etap,dt,vx(*),vy(*),vz(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i
      REAL*8  scale

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      scale=DEXP(-dt*etap)
      DO i=noff,n
         vx(i)=vx(i)*scale
         vy(i)=vy(i)*scale
         vz(i)=vz(i)*scale
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
