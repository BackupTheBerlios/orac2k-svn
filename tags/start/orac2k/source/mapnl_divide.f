      SUBROUTINE mapnl_divide(node,nstart,nend,grppt,mapnl)

************************************************************************
*   Time-stamp: <99/01/31 17:54:15 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Jul 15 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nstart,nend,node
      INTEGER grppt(2,*),mapnl(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,l,m,map,map1

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      map1=0
      DO i=1,nstart-1
         DO j=grppt(1,i),grppt(2,i)
            m=mapnl(map1+1)
            map1=map1+m+1
         END DO
      END DO
      map=0
      DO i=nstart,nend
         DO j=grppt(1,i),grppt(2,i)
            m=mapnl(map1+1)
            mapnl(map+1)=m
            DO l=1,m
               mapnl(map+1+l)=mapnl(map1+1+l)
            END DO
            map=map+m+1
            map1=map1+m+1
         END DO
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
