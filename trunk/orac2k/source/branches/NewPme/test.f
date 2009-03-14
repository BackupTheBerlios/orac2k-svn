      PROGRAM test

************************************************************************
*   Time-stamp: <2008-03-14 13:07:58 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Fri Mar 14 2008 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program  ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

*----------------------- VARIABLES IN COMMON --------------------------*

*------------------------- LOCAL VARIABLES ----------------------------*
      REAL(8), ALLOCATABLE :: pla(:,:),vrt(:,:,:),area(:),d2(:)
      INTEGER, ALLOCATABLE ::  nver(:)


*----------------------- EXECUTABLE STATEMENTS ------------------------*


      CALL comp_voronoi


*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      END
