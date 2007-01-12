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

  USE Cell
  USE Solvent
  USE Solute
  USE PDB
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
    TYPE(ResiduePDB), POINTER :: PDB_Slv(:), PDB_Slt(:)
    INTEGER :: n,m

    IF(ALLOCATED(Secondary(1) % Line)) THEN
       IF(.NOT. ALLOCATED(PDB_Solute)) THEN
          errmsg_f='Solute is defined, but no file to read from '
          CALL Add_Errors(-1,errmsg_f)
          CALL Print_Errors()
       END IF
       IF(.NOT. PDB_('Solute ',PDB_Solute, PDB_Slt)) THEN
          CALL Print_Errors()
       ELSE
          IF(.NOT. PDB__Validate('Solute',PDB_Slt)) CALL Print_Errors()
       END IF
       
       
       DO n=1,SIZE(PDB_Slt)
          IF(ALLOCATED(PDB_Slt(n) % atm)) THEN
             WRITE(*,*) PDB_Slt(n) % No, PDB_Slt(n) % ResName
             DO m=1,SIZE(PDB_Slt(n) % atm)
                WRITE(*,*) PDB_Slt(n) % atm(m) % x&
                     &, PDB_Slt(n) % atm(m) % y, PDB_Slt(n) % atm(m) % z&
                     &, PDB_Slt(n) % atm(m) % Serial&
                     &, PDB_Slt(n) % atm(m) % AtmName
             END DO
          END IF
       END DO
    END IF
  END SUBROUTINE SimulationBox_
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE SimulationBox
