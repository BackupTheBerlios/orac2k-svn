      SUBROUTINE setup_fake_map(ngrp,grppt,nnlppf,pfix,mass)

************************************************************************
*   Time-stamp: <04/12/16 15:05:27 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Feb 28 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER ngrp,grppt(2,*)
      INTEGER nnlppf(*)
      REAL*8  mass(*)
      LOGICAL pfix
      
*----------------------- VARIABLES IN COMMON --------------------------*

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER count,j,i
      LOGICAL ok

*----------------------- EXECUTABLE STATEMENTS ------------------------*
      
      count=0
      DO j=1,ngrp
         ok=.TRUE.
         IF(pfix) THEN
            DO i=grppt(1,j),grppt(2,j)
               IF(mass(i) .GT. 1.0D5) ok=.FALSE.
            END DO
         END IF
         IF(ok) THEN
            count=count+1
            nnlppf(1+count)=j
         END IF
      END DO
      nnlppf(1)=count      

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
