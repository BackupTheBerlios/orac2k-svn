!!$/---------------------------------------------------------------------\
!!$                                                                      |
!!$  Copyright (C) 2006-2007 Massimo Marchi <Massimo.Marchi@cea.fr>      |
!!$                                                                      |
!!$      This program is free software;  you  can  redistribute  it      |
!!$      and/or modify it under the terms of the GNU General Public      |
!!$      License version 2 as published  by  the  Free  Software         |
!!$      Foundation;                                                     |
!!$                                                                      |
!!$      This program is distributed in the hope that  it  will  be      |
!!$      useful, but WITHOUT ANY WARRANTY; without even the implied      |
!!$      warranty of MERCHANTABILITY or FITNESS  FOR  A  PARTICULAR      |
!!$      PURPOSE.   See  the  GNU  General  Public License for more      |
!!$      details.                                                        |
!!$                                                                      |
!!$      You should have received a copy of the GNU General  Public      |
!!$      License along with this program; if not, write to the Free      |
!!$      Software Foundation, Inc., 59  Temple  Place,  Suite  330,      |
!!$      Boston, MA  02111-1307  USA                                     |
!!$                                                                      |
!!$\---------------------------------------------------------------------/
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
  PUBLIC :: SystemPrm_, SystemPrm__Type,SystemPrm__Chain, Prm, SystemPrm__Update
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
    Prm % Types=>TypesPrm_()
    IF(.NOT. ASSOCIATED(Prm % Types)) CALL Print_errors()
    Prm % LJ=>LennardJones_() 
    IF(.NOT. ASSOCIATED(Prm % LJ)) CALL Print_errors()
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
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE SystemPrm
