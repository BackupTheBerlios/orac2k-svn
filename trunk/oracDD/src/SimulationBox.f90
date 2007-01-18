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
MODULE SimulationBox

!!$***********************************************************************
!!$   Time-stamp: <2007-01-12 20:39:55 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Jan 12 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program oracDD ----*


!!$======================== DECLARATIONS ================================*


  USE Neighbors
  USE AddHydrogens_
  USE SystemTpg
  USE Cell
  USE Solvent
  USE Solute
  USE PDB
  USE IndSequence
  USE SecondarySeq
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f
  USE SystemPrm
  USE AtomCnt
  USE AtomBox
  IMPLICIT none
  PRIVATE
  PUBLIC :: SimulationBox_
  TYPE(AtomBox__), POINTER, SAVE :: Total(:)
  TYPE(AtomBox__), POINTER, SAVE :: Slv(:),Slt(:)
  REAL(8), PARAMETER :: Cube_length=10.0D0
  REAL(8), SAVE :: rcut
CONTAINS
  SUBROUTINE SimulationBox_
    INTEGER :: n,m,Begins, Ends,p, nccx,nccy,nccz
    
    TYPE(AtomPdb), POINTER :: PDB__Coords(:)=>NULL()

    IF(ALLOCATED(Secondary(1) % Line)) THEN
       IF(.NOT. ALLOCATED(PDB_Solute)) THEN
          errmsg_f='Solute is defined, but no file to read from '
          CALL Add_Errors(-1,errmsg_f)
          CALL Print_Errors()
       END IF
       IF(.NOT. PDB_('Solute', PDB_Solute, PDB__Coords)) CALL Print_Errors()
       IF(.NOT. AddHydrogens__(PDB__Coords)) CALL Print_Errors()

       
       IF(ASSOCIATED(PDB__Coords)) THEN
          CALL PDB__Write(PDB__Coords)
          IF(ALLOCATED(Secondary(2) % Line)) THEN
             CALL AtomBox_(PDB__Coords,Slt)
          END IF
          DEALLOCATE(PDB__Coords)
       END IF
    END IF

    IF(ALLOCATED(Secondary(2) % Line)) THEN
       IF(.NOT. ALLOCATED(PDB_Solvent)) THEN
          errmsg_f='Solvent is defined, but no file to read from '
          CALL Add_Errors(-1,errmsg_f)
          CALL Print_Errors()
       END IF
       IF(.NOT. PDB_('Solvent', PDB_Solvent, PDB__Coords)) CALL Print_Errors()
       IF(.NOT. AddHydrogens__(PDB__Coords)) CALL Print_Errors()
       IF(ASSOCIATED(PDB__Coords)) THEN
          CALL PDB__Write(PDB__Coords)
          CALL AtomBox_(PDB__Coords,Slv)
          DEALLOCATE(PDB__Coords)
       END IF

       IF(Solvent__Param % added /= 0) RETURN

       IF(.NOT. Solvent__Param % Build) RETURN
       IF(.NOT. AtomBox__BuildSlv(Slv)) CALL Print_Errors()

       IF(ALLOCATED(Secondary(1) % Line)) THEN
          n=SIZE(Slt)
          m=SIZE(Slv)
          ALLOCATE(Total(n+m))
          Total(1:n)=Slt
          Total(n+1:)=Slv
          DEALLOCATE(Slt,Slv)
          Slt=>Total(1:n)
          Slv=>Total(n+1:)
       ELSE
          m=SIZE(Slv)
          ALLOCATE(Total(m))
          DEALLOCATE(Slv)
          Slv=>Total
       END IF
       nccx=INT(a/Cube_Length)
       nccy=INT(b/Cube_Length)
       nccz=INT(c/Cube_Length)
       rcut=MAXVAL(Slv(:) % sigma)*2.0D0
       CALL AtomBox__ChgFrame(-1,Total)

       IF(.NOT. Neighbors_(rcut, nccx, nccy, nccz)) CALL Print_Errors()
       IF(.NOT. Neighbors__Atoms(Total(:) % x&
            &, Total(:) % y, Total(:) % z)) CALL Print_Errors()
       
       CALL AtomBox__ChgFrame(1,Total)
       DO n=1,SIZE(Total)
          WRITE(98,'(i6,2x,f9.3,3x,3f12.4)') Total(n) % Serial&
               &, Total(n) % Sigma, Total(n) % x, Total(n) % y, Total(n) % z
       END DO
    END IF


  END SUBROUTINE SimulationBox_
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE SimulationBox
