      SUBROUTINE zero_voronoi

************************************************************************
*   Time-stamp: <2005-01-27 17:52:05 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Jun 21 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

*----------------------- VARIABLES IN SCRATCH COMMON ------------------*
      
      INCLUDE 'voronoi.h'

*------------------------- LOCAL VARIABLES ----------------------------*


*----------------------- EXECUTABLE STATEMENTS ------------------------*

      volume_vor=0.0D0
      area_vor=0.0D0
      nnlpp_vor=0

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
