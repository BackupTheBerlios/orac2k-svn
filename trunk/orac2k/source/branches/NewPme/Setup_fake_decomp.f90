SUBROUTINE Setup_fake_decomp()

!!$***********************************************************************
!!$   Time-stamp: <2005-01-29 12:44:18 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Sat Jan 29 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*

  IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

!!$----------------------- VARIABLES IN COMMON --------------------------*

!!$------------------------- LOCAL VARIABLES ----------------------------*


!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

  nstart_h=1
  nend_h=ngrp
  nlocal_h=ngrp
  nstart_ah=1
  nend_ah=ntap
  nlocal_ah=ntap
  


!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

  
END SUBROUTINE Setup_fake_decomp
