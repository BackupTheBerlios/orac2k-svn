      SUBROUTINE dcopy3(n,ax,ay,az,na,bx,by,bz,nb)

************************************************************************
*   Time-stamp: <98/01/09 18:04:40 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Fri Jan  9 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER n,na,nb
      REAL*8  ax(*),ay(*),az(*),bx(*),by(*),bz(*)

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      CALL dcopy(n,ax,na,bx,nb)
      CALL dcopy(n,ay,na,by,nb)
      CALL dcopy(n,az,na,bz,nb)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
