      SUBROUTINE correc_scale_2(nstart,nend,pointer,vx,vy,vz,etap,dt)

************************************************************************
*   Time-stamp: <99/02/18 13:27:08 marchi>                             *
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

      INTEGER pointer(*),nstart,nend
      REAL*8  etap(*),dt,vx(*),vy(*),vz(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,m
      REAL*8  scale(2)

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      scale(1)=DEXP(-dt*etap(1))
      scale(2)=DEXP(-dt*etap(2))
      DO i=nstart,nend
         j=pointer(i)
         vx(i)=vx(i)*scale(j)
         vy(i)=vy(i)*scale(j)
         vz(i)=vz(i)*scale(j)
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
