!!$/---------------------------------------------------------------------\
!!$   Copyright  © 2006-2007 Massimo Marchi <Massimo.Marchi at cea.fr>   |
!!$                                                                      |
!!$    This software is a computer program named oracDD whose            |
!!$    purpose is to simulate and model complex molecular systems.       |
!!$    The code is written in fortran 95 compliant with Technical        |
!!$    Report TR 15581, and uses MPI-1 routines for parallel             |
!!$    coding.                                                           |
!!$                                                                      |
!!$    This software is governed by the CeCILL license under             |
!!$    French law and abiding by the rules of distribution of            |
!!$    free software.  You can  use, modify and/ or redistribute         |
!!$    the software under the terms of the CeCILL icense as              |
!!$    circulated by CEA, CNRS and INRIA at the following URL            |
!!$    "http://www.cecill.info".                                         |
!!$                                                                      |
!!$    As a counterpart to the access to the source code and rights      |
!!$    to copy, modify and redistribute granted by the license,          |
!!$    users are provided only with a limited warranty and the           |
!!$    software's author, the holder of the economic rights, and         |
!!$    the successive licensors have only limited liability.             |
!!$                                                                      |
!!$    The fact that you are presently reading this means that you       |
!!$    have had knowledge of the CeCILL license and that you accept      |
!!$    its terms.                                                        |
!!$                                                                      |
!!$    You should have received a copy of the CeCILL license along       |
!!$    with this program; if not, you can collect copies on the URL's    |
!!$    "http://www.cecill.info/licences/Licence_CeCILL_V2-en.html"       |
!!$    "http://www.cecill.info/licences/Licence_CeCILL_V2-fr.html"       |
!!$                                                                      |
!!$----------------------------------------------------------------------/
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
  PUBLIC :: SystemPrm_, SystemPrm__Type, Prm, SystemPrm__Update&
       &, SystemPrm__Read, SystemPrm__Write
  TYPE :: SystemPrm__Type
     TYPE(SystemPrm__Chain), POINTER :: bonds(:)=>NULL(),angles(:)=>NULL(),&
          &dihed(:)=>NULL(),imph(:)=>NULL()
     TYPE(LennardJones__Type), POINTER :: LJ=>NULL()
     CHARACTER(len=max_atm), POINTER :: Types(:)=>NULL()
  END TYPE SystemPrm__Type

  TYPE(SystemPrm__Type), SAVE :: Prm

  REAL(8), DIMENSION(:), ALLOCATABLE  :: Params
  LOGICAL, DIMENSION(:), ALLOCATABLE  :: oks
CONTAINS
  SUBROUTINE SystemPrm_
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
  SUBROUTINE SystemPrm__Update
    Prm % bonds=> BondsPrm_()
    IF(.NOT. ASSOCIATED(Prm % bonds)) CALL Print_errors()
    Prm % angles=>AnglesPrm_()
    IF(.NOT. ASSOCIATED(Prm % angles)) CALL Print_errors()
    Prm % dihed=>TorsionsPrm_()
    IF(.NOT. ASSOCIATED(Prm % dihed)) CALL Print_errors()
    Prm % imph=>ImpropersPrm_()
    IF(.NOT. ASSOCIATED(Prm % imph)) CALL Print_errors()
  END SUBROUTINE SystemPrm__Update
  SUBROUTINE SystemPrm__Write
    CALL TypesPrm__Write
    CALL LennardJones__Write
    CALL BondsPrm__Write
    CALL AnglesPrm__Write
    CALL TorsionsPrm__Write
    CALL ImpropersPrm__Write
  END SUBROUTINE SystemPrm__Write
  SUBROUTINE SystemPrm__Read
    Prm % Types=>TypesPrm__Read()
    IF(.NOT. ASSOCIATED(Prm % Types)) CALL Print_errors()
    Prm % LJ=>LennardJones__Read() 
    IF(.NOT. ASSOCIATED(Prm % LJ)) CALL Print_errors()
    Prm % bonds=> BondsPrm__Read()
    IF(.NOT. ASSOCIATED(Prm % bonds)) CALL Print_errors()
    Prm % angles=>AnglesPrm__Read()
    IF(.NOT. ASSOCIATED(Prm % angles)) CALL Print_errors()
    Prm % dihed=>TorsionsPrm__Read()
    IF(.NOT. ASSOCIATED(Prm % dihed)) CALL Print_errors()
    Prm % imph=>ImpropersPrm__Read()
    IF(.NOT. ASSOCIATED(Prm % imph)) CALL Print_errors()
  END SUBROUTINE SystemPrm__Read
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE SystemPrm
