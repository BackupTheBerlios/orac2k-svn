MODULE ORAC_PROCESS_Mod

!!$***********************************************************************
!!$   Time-stamp: <2006-12-15 13:58:50 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Tue Nov 21 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program  ----*

  USE PARAMETERS_Mod, ONLY: Scan_Parameters
  USE ERROR_Mod, ONLY: Print_Errors, Add_Errors=>Add,Setup_Errors
  USE MYPARSE_Mod, ONLY: max_pars
  USE CLASS_Tree
  IMPLICIT none
  PRIVATE
  PUBLIC Inputs,Grammar,Process_Commands
  TYPE(tree), POINTER, SAVE :: Inputs,Grammar
  CHARACTER(len=max_pars), DIMENSION(:), ALLOCATABLE, PRIVATE :: strngs
CONTAINS
  SUBROUTINE Process_Commands
    TYPE(Branch), POINTER :: checks,heldo
    INTEGER :: o
    CALL Setup_Errors
    CALL Get_Tree(Inputs)
    CALL Scan_Parameters
    CALL Print_Errors()
  END SUBROUTINE Process_Commands
END MODULE ORAC_PROCESS_Mod
