      SUBROUTINE P_get_intra(intra,p1,nintra,intrax,nstart,nend,iret)

************************************************************************
*   Time-stamp: <99/02/23 18:27:45 marchi>                             *
*                                                                      *
*   intrax: The intramolecular interactions having at least an atom    *
*           between nstart and nend                                    *
*   intrah: The intramolecular interactions with at least an atom      *
*           outside the bondary                                        *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Fri Jul 10 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER p1,intra(p1,*),intrax(*),nstart,nend,nintra,iret
      
*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,count
      LOGICAL ok

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      count=0
      DO i=1,nintra
         ok=.FALSE.
         DO j=1,p1
            IF(intra(j,i) .GE. nstart .AND. intra(j,i) .LE. nend) ok=
     &           .TRUE.
         END DO
         IF(ok) THEN
            count=count+1
            intrax(1+count)=i
         END IF
      END DO
      intrax(1)=count

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
