MODULE SYSTEM_Mod

!!$***********************************************************************
!!$   Time-stamp: <2006-12-20 21:53:59 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Wed Dec 20 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*

  USE ERROR_Mod, ONLY: Add_Errors=>Add, Print_Errors, error_args, errmsg_f
  USE PDB_Utils, ONLY: PDB__Coordinates=>Coordinates, ResiduePDB
  USE Setup__solute, ONLY: PDB_Solute, PDB_Template
  USE Class_Atom, ONLY: Atom__Init=>Init, Atom  
  IMPLICIT none
  TYPE(ResiduePDB), DIMENSION(:), ALLOCATABLE, SAVE :: Res_Solute, Res_Template&
       &, Res_Solvent
CONTAINS
  SUBROUTINE Setup
    CALL Atom__Init
    IF(ALLOCATED(PDB_Solute))&
         & CALL PDB__Coordinates('Solute',PDB_Solute, Res_Solute)
    CALL Print_Errors()
  END SUBROUTINE Setup

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE SYSTEM_Mod
