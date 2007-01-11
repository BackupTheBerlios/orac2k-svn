MODULE SystemPrm

!!$***********************************************************************
!!$   Time-stamp: <2007-01-11 19:18:18 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Tue Nov 28 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program ORAC ----*
  
  USE Constants
  USE PrmUtilities
  USE LennardJones
  USE TypesPrm
  USE BondsPrm
  USE AnglesPrm
  USE TorsionsPrm
  USE ImpropersPrm
  USE SystemTpg
  USE AtomCnt
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f
  IMPLICIT none
  PRIVATE
  PUBLIC :: SystemPrm_, SystemPrm__Type,SystemPrm__Chain, Prm
  TYPE :: SystemPrm__Type
     TYPE(SystemPrm__Chain), POINTER :: bonds(:)=>NULL(),angles(:)=>NULL(),&
          &dihed(:)=>NULL(),imph(:)=>NULL()
!!$     ,int14(:)=>NULL()
     TYPE(LennardJones__Type), POINTER :: LJ=>NULL()
     CHARACTER(len=max_atm), POINTER :: Types(:)=>NULL()
  END TYPE SystemPrm__Type

  TYPE(SystemPrm__Type), SAVE :: Prm

  REAL(8), DIMENSION(:), ALLOCATABLE  :: Params
  LOGICAL, DIMENSION(:), ALLOCATABLE  :: oks
CONTAINS
  SUBROUTINE SystemPrm_
    INTEGER :: n,m
    Prm % Types=>TypesPrm_()
    IF(.NOT. ASSOCIATED(Prm % Types)) CALL Print_errors()
    Prm % LJ=>LennardJones_() 
    IF(.NOT. ASSOCIATED(Prm % LJ)) CALL Print_errors()
    Prm % bonds=> BondsPrm_()
    IF(.NOT. ASSOCIATED(Prm % bonds)) CALL Print_errors()
    Prm % angles=>AnglesPrm_()
    IF(.NOT. ASSOCIATED(Prm % angles)) CALL Print_errors()
    Prm % dihed=>TorsionsPrm_()
    IF(.NOT. ASSOCIATED(Prm % dihed)) CALL Print_errors()
    Prm % imph=>ImpropersPrm_()
    IF(.NOT. ASSOCIATED(Prm % imph)) CALL Print_errors()

  END SUBROUTINE SystemPrm_
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE SystemPrm
