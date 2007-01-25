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
MODULE Atom
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Jan 25 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*

  USE Constants
  USE Errors, ONLY: Add_Errors=>Add, errmsg_f
  USE SystemTpg
  USE Cell
  USE AtomCnt
  USE AtomBox
  USE SimulationBox
  USE PDB
  IMPLICIT none
  PRIVATE
  PUBLIC Atom_, Atom__PDB, Atom__InitCoords
  TYPE :: Atom__
     CHARACTER(len=max_atm) :: Res, beta, betab
     INTEGER :: Res_No, Grp_No, Id_Res, Id_Type, Id_slv, Mol
     REAL(8) :: x,y,z      ! Coordinates orthogonal frame
     REAL(8) :: xa,ya,za   ! Coordinates reduced frame
     REAL(8) :: vx,vy,vz      ! Velocities orthogonal frame
     REAL(8) :: vxa,vya,vza   ! Velocities reduced frame
     REAL(8) :: chg, mass
  END TYPE Atom__
  TYPE(Atom__), ALLOCATABLE, SAVE :: Atoms(:)
CONTAINS
  FUNCTION Atom_() RESULT(out)
    LOGICAL :: out
    INTEGER :: n,m,s,nato

    out=.TRUE.
    IF(.NOT. ALLOCATED(AtomCnts)) THEN
       errmsg_f='The No. of Atoms of the system is unknown'
       CALL Add_Errors(-1,errmsg_f)
       out=.FALSE.
       RETURN
    END IF
    nato=SIZE(AtomCnts)
    ALLOCATE(Atoms(nato))
    DO n=1,nato
       Atoms(n) % Grp_No = AtomCnts(n) % Grp_No
       Atoms(n) % Res_No = AtomCnts(n) % Res_No
       Atoms(n) % Id_Type = AtomCnts(n) % Id_Type
       Atoms(n) % Id_Slv = AtomCnts(n) % Id_Slv
       Atoms(n) % Id_Res = AtomCnts(n) % Id_Res
       Atoms(n) % chg = AtomCnts(n) % chg
       Atoms(n) % mass = AtomCnts(n) % mass
       Atoms(n) % Res = AtomCnts(n) % Res
       Atoms(n) % beta = AtomCnts(n) % Beta
       Atoms(n) % betab = AtomCnts(n) % Betab
    END DO
    DO n=1,SIZE(Tpg % Mol_Atm)
       DO m=1,SIZE(Tpg % Mol_Atm(n) % g)
          s=Tpg % Mol_Atm(n) % g (m)
          Atoms(s) % Mol = n
       END DO
    END DO    
  END FUNCTION Atom_
  FUNCTION Atom__InitCoords() RESULT(out)
    LOGICAL :: out
    INTEGER :: n

    out=.TRUE.
    IF(.NOT. ALLOCATED(Atoms_InBox)) THEN
       errmsg_f='Cannot initialize atomic coordinates: Initial&
            & coordinates were expected, but none were either &
            &built or read in by SimulationBox Module.'
       CALL Add_Errors(-1,errmsg_f)
       out=.FALSE.
       RETURN
    END IF
    
    Atoms(:) % xa = Atoms_InBox(:) % x 
    Atoms(:) % ya = Atoms_InBox(:) % y 
    Atoms(:) % za = Atoms_InBox(:) % z 
    Atoms(:) % x = co(1,1)*Atoms(:) % xa+co(1,2)*Atoms(:) % ya+co(1,3)*Atoms(:) % za    
    Atoms(:) % y = co(2,1)*Atoms(:) % xa+co(2,2)*Atoms(:) % ya+co(2,3)*Atoms(:) % za    
    Atoms(:) % z = co(3,1)*Atoms(:) % xa+co(3,2)*Atoms(:) % ya+co(3,3)*Atoms(:) % za    

  END FUNCTION Atom__InitCoords
  SUBROUTINE Atom__PDB(unit)
    INTEGER :: unit
    TYPE(AtomPdb), ALLOCATABLE :: PDB__Coords(:)
    INTEGER :: n

    ALLOCATE(PDB__Coords(SIZE(Atoms)))
    DO n=1,SIZE(Atoms)
       PDB__Coords(n) % x = Atoms(n) % x
       PDB__Coords(n) % y = Atoms(n) % y
       PDB__Coords(n) % z = Atoms(n) % z
       PDB__Coords(n) % Serial =n 
    END DO
    CALL PDB__Write(unit,PDB__Coords)
    DEALLOCATE(PDB__Coords)
  END SUBROUTINE Atom__PDB
END MODULE Atom
