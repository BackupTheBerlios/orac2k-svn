!!$***********************************************************************
!!$   Time-stamp: <2005-02-25 11:59:44 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Feb 25 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This MODULE is part of the program ORAC ----*
MODULE Extra_Forces
  TYPE extra_stretch
     REAL(8) :: K, r0
  END TYPE extra_stretch
  INCLUDE 'unit.h'
  TYPE(extra_stretch) :: pop 
CONTAINS
  SUBROUTINE Init
!!$======================== DECLARATIONS ================================*

  IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

!!$----------------------- VARIABLES IN COMMON --------------------------*

!!$------------------------- LOCAL VARIABLES ----------------------------*


!!$----------------------- EXECUTABLE STATEMENTS ------------------------*





!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE Extra_Forces
