      SUBROUTINE comp_noneq_free(vpx,vpy,vpz,fpx,fpy,fpz,nstart,nend
     &     ,dfree)

************************************************************************
*   Time-stamp: <99/04/12 13:29:37 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Apr  8 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nstart,nend
      REAL*8  dfree,vpx(*),vpy(*),vpz(*),fpx(*),fpy(*),fpz(*)

*----------------------- VARIABLES IN COMMON --------------------------*

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      dfree=0.0D0
      DO i=nstart,nend
         dfree=dfree+vpx(i)*fpx(i)+vpy(i)*fpy(i)+vpz(i)*fpz(i)
      END DO
      dfree=-dfree

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
