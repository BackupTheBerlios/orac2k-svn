      SUBROUTINE remove_mv(fpx,fpy,fpz,mass,ntap)

************************************************************************
*   Time-stamp: <98/02/10 16:09:11 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Feb 10 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8  fpx(*),fpy(*),fpz(*),mass(*)
      INTEGER ntap

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i
      REAL*8  totmass,fxtot,fytot,fztot,xm

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      fxtot=0.0D0
      fytot=0.0D0
      fztot=0.0D0
      totmass=0.0D0
      DO i=1,ntap
         fxtot=fxtot+fpx(i)
         fytot=fytot+fpy(i)
         fztot=fztot+fpz(i)
         totmass=totmass+mass(i)
      END DO
      DO i=1,ntap
         xm=mass(i)/totmass
         fpx(i)=fpx(i)-fxtot*xm
         fpy(i)=fpy(i)-fytot*xm
         fpz(i)=fpz(i)-fztot*xm
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
