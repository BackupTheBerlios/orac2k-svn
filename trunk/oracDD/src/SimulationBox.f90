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

  
  USE SystemTpg
  USE Cell
  USE Solvent
  USE Solute
  USE PDB
  USE IndSequence
  USE SecondarySeq
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f
  IMPLICIT none

  PRIVATE
  PUBLIC :: SimulationBox_
  TYPE :: AtomsBox
     REAL(8) :: x,y,z
  END TYPE AtomsBox
  TYPE(AtomsBox), ALLOCATABLE, SAVE :: Slv(:),Slt(:)
CONTAINS
  SUBROUTINE SimulationBox_
    INTEGER :: n,m,Begins, Ends,p 
    TYPE(AtomPdb), POINTER :: PDB__Coords(:)

    IF(ALLOCATED(Secondary(1) % Line)) THEN
       IF(.NOT. ALLOCATED(PDB_Solute)) THEN
          errmsg_f='Solute is defined, but no file to read from '
          CALL Add_Errors(-1,errmsg_f)
          CALL Print_Errors()
       END IF
       IF(.NOT. PDB_('Solute', PDB_Solute, PDB__Coords)) THEN
          CALL Print_Errors()
       END IF
       
       
       IF(ASSOCIATED(PDB__Coords)) THEN
          DO n=1,SIZE(PDB__Coords)
             p=PDB__Coords(n) % Serial
             IF(p ==0) CYCLE
             WRITE(*,'(i5,2x,i5,4x,a5,3f10.4,2x,a5,2x,a5)') &
                  &  n,p&
                  &, TRIM(Tpg % atm(p) % a % Res) &
                  &, PDB__Coords(n) % x &
                  &, PDB__Coords(n) % y &
                  &, PDB__Coords(n) % z &
                  &, TRIM(Tpg % atm(p) % a % beta) & 
                  &, TRIM(Tpg % atm(p) % a % betab) 
          END DO
       END IF
    END IF
  END SUBROUTINE SimulationBox_
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE SimulationBox
