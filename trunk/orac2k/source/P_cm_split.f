      SUBROUTINE P_cm_split(node,nstart,nend,nlocal,nstart_cm,nend_cm
     &     ,nlocal_cm,atomp)

************************************************************************
*   Time-stamp: <1999-10-18 14:29:38 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Jan 23 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER node,nstart,nend,nlocal,nstart_cm,nend_cm,nlocal_cm,atomp(
     &     *)

*------------------------- LOCAL VARIABLES ----------------------------*


*----------------------- EXECUTABLE STATEMENTS ------------------------*


*=======================================================================
*---- Split groups center of mass with the above choice of nstart ------
*---- and nend ---------------------------------------------------------
*=======================================================================

      IF(nlocal .NE. 0) THEN
         nstart_cm=atomp(nstart)
         nend_cm=atomp(nend)
         nlocal_cm=nend_cm-nstart_cm+1
      ELSE
         nstart_cm=1
         nend_cm=0
         nlocal_cm=0
      END IF
      
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
