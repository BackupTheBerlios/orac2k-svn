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
MODULE PDB

!!$***********************************************************************
!!$   Time-stamp: <2007-01-13 01:09:49 marchi>                           *
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

  USE IndSequence
  USE SystemTpg
  USE Constants, ONLY: max_pars,max_data,max_char, Used
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_w, errmsg_f
  USE Strings, ONLY: My_Fam, MyRead, MyPutNum
  USE Node
  USE MyParse
  USE STRPAK
  IMPLICIT none
  PRIVATE
  PUBLIC PDB_, AtomPdb, PDB__Write
  TYPE :: AtomPDB
     INTEGER :: Serial
     REAL(8) :: x,y,z
     CHARACTER(len=4) :: AtmName
  END TYPE AtomPDB
  TYPE ResiduePDB
     INTEGER :: No
     CHARACTER(len=3) :: ResName
     TYPE(AtomPdb), DIMENSION(:), ALLOCATABLE :: atm
  END TYPE ResiduePDB
  INTEGER, SAVE, POINTER :: Res_Atm(:,:)=>NULL()
  INTEGER, SAVE, POINTER :: Grp_Atm(:,:)=>NULL()
  TYPE :: Name_Exception
     SEQUENCE
     CHARACTER(len=max_char), ALLOCATABLE :: res(:)
     CHARACTER(len=max_char), ALLOCATABLE :: lab(:)
  END TYPE Name_Exception
  TYPE(Name_Exception), SAVE :: ex(4)
  TYPE(ResiduePDB), ALLOCATABLE :: ResPdb(:)
  INTEGER, SAVE, POINTER :: SltSlv(:,:)=>NULL()
  INTEGER, SAVE :: Res_Begins,Res_Ends
CONTAINS
  FUNCTION PDB_(Type, PDB_String, PDB__Coords) RESULT(out)
    LOGICAL :: out
    CHARACTER(len=*) :: Type
    INTEGER :: Res_Begins, Res_End
    CHARACTER(len=max_char)  :: PDB_string(:)
    TYPE(AtomPDB), POINTER :: PDB__Coords(:)

    out=.TRUE.
    CALL PDB__Init(Type, PDB__Coords)
    IF(.NOT. PDB__Read(Type, PDB_String)) THEN
       out=.FALSE.
       RETURN
    ELSE IF(.NOT. PDB__Validate(Type, PDB__Coords)) THEN
       out=.FALSE.
       RETURN
    END IF
  END FUNCTION PDB_
  SUBROUTINE PDB__Init(Type, PDB__Coords)
    CHARACTER(len=*) :: Type
    TYPE(AtomPDB), POINTER :: PDB__Coords(:)
    INTEGER :: nato
!!$
!!$--- Get the beginning and the end of each residue atom
!!$    
    Res_Atm=>IndSequence__Res()

!!$
!!$--- Get the beginning and the end of residue for each PDB Type
!!$    
    SltSlv=>IndSequence__sltslv_Res()
    IF(TRIM(Type) == 'Solute' .OR. TRIM(Type) == 'Template') THEN
       Res_Begins=SltSlv(1,1)
       Res_Ends=SltSlv(2,1)
    ELSE IF(TRIM(Type) == 'Solvent') THEN
       Res_Begins=SltSlv(1,2)
       Res_Ends=SltSlv(2,2)
    ELSE IF(TRIM(Type) == ' ') THEN
       Res_Begins=1
       Res_Ends=SIZE(Res_Atm,2)
    END IF

    nato=Res_Atm(2,Res_Ends)-Res_Atm(1,Res_Begins)+1
    IF(ASSOCIATED(PDB__Coords)) DEALLOCATE(PDB__Coords)
    ALLOCATE(PDB__Coords(nato))
    PDB__Coords(:) % x = 1.0D10
    PDB__Coords(:) % y = 1.0D10
    PDB__Coords(:) % z = 1.0D10
    PDB__Coords(:) % AtmName = 'h'    
    PDB__Coords(:) % Serial = 0

  END SUBROUTINE PDB__Init

  INCLUDE 'PDB__Read.f90'
  INCLUDE 'PDB__Write.f90'
  INCLUDE 'PDB__Validate.f90'
END MODULE PDB
