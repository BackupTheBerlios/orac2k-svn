MODULE Process

!!$***********************************************************************
!!$   Time-stamp: <2007-01-04 18:04:11 marchi>                           *
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

  USE Parameters
  USE Setup
  USE Grammars
  USE Errors, ONLY: Print_Errors, Add_Errors=>Add,Setup_Errors
  USE Tree
  IMPLICIT none
  PRIVATE
  PUBLIC Process_
CONTAINS
  SUBROUTINE Process_
    INTEGER :: o
    CALL Setup_Errors
    CALL Tree__Get_Tree(Grammars__Inputs)
    CALL Setups__Scan
    CALL Parameters__Scan
    CALL Print_Errors()
  END SUBROUTINE Process_
END MODULE Process
