      SUBROUTINE gen_stress_tag_bnd(na,nb,index,ns,atomp,tag)

************************************************************************
*   Time-stamp: <98/12/02 13:14:16 marchi>                             *
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

      INTEGER na,nb
      INTEGER index(nb,*),ns(*),atomp(*),tag(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,l,ia,ib,ii,ma
      LOGICAL ok

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      ma=ns(1)
      DO ii=1,ma
         i=IABS(ns(ii+1))
         ia=atomp(index(1,i))
         ok=.TRUE.
         DO l=2,nb
            ib=atomp(index(l,i))
            IF(ib .NE. ia) ok=.FALSE.
            ia=ib
         END DO
         IF(ok) THEN
            tag(i)=1
         ELSE
            tag(i)=2
         END IF
      END DO
      
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
