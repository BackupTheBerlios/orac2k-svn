      FUNCTION get_pointer_protl(nprot,protl)

************************************************************************
*   Time-stamp: <99/01/30 16:28:32 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Jan 30 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER get_pointer_protl,nprot,protl(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,count,m

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      count=0
      DO i=1,nprot
         m=protl(count+1)
         count=count+m+1
      END DO
      get_pointer_protl=count+1

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
