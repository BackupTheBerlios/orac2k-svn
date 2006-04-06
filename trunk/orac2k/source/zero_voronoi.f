      SUBROUTINE zero_voronoi

************************************************************************
*   Time-stamp: <2006-04-05 14:31:18 marchi>                             *
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

      USE VORONOI_Mod, ONLY: volume_vor,area_vor,nnlpp_vor
      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

*----------------------- VARIABLES IN SCRATCH COMMON ------------------*
      

*------------------------- LOCAL VARIABLES ----------------------------*


*----------------------- EXECUTABLE STATEMENTS ------------------------*

      volume_vor=0.0D0
      area_vor=0.0D0
      nnlpp_vor=0

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
