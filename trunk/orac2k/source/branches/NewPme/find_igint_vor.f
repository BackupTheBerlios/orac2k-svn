      SUBROUTINE find_igint_vor(ngrp,grppt,iret,errmsg)

************************************************************************
*   Time-stamp: <99/05/14 13:41:56 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed Jun 18 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INCLUDE 'parst.h'
      INCLUDE 'voronoi.h'
      INTEGER ngrp,grppt(2,*),iret
      CHARACTER*80 errmsg
      
*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i1,j1,n,map

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO n=1,ngrp
         DO i1=grppt(1,n),grppt(2,n)
            map=0
            DO j1=grppt(1,n),grppt(2,n)
               map=map+1
               IF(map+1 .GT. pig_nnl) THEN
                  iret=1
                  errmsg=
     &   'While counting bonded neighbors, _MAX_IG_NNL_ in config.h'
     &   / /' is too small. Abort.'
                  RETURN
               END IF
               ig_nnl(map+1,i1)=j1
            END DO
            ig_nnl(1,i1)=map
         END DO
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
