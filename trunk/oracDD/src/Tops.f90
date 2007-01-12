MODULE Tops

!!$***********************************************************************
!!$   Time-stamp: <2007-01-12 12:27:25 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Tue Nov 14 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*
  

  USE Solvent
  USE SecondarySeq
  USE Types
  USE Myparse
  USE Strings, ONLY: MY_Fam, My_Fxm,ST_Concat
  USE Hash_Tops
  USE Resid, ONLY: Resids, Residue=>Topology
  USE Node
  USE STRPAK
  USE Parameters_Globals
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f
  PRIVATE
  PUBLIC Tops_, Tops__Type, App_Char, Top__Store
  TYPE :: Tops__Type
     CHARACTER(len=max_char) :: FField=' ',Type=' ',Residue=' '
     CHARACTER(len=max_char) , DIMENSION(:,:), ALLOCATABLE :: bonds,imph&
          &,acc,acc_,don,don_,dele,ends
     CHARACTER(len=max_char), DIMENSION(:,:), ALLOCATABLE :: mass
     TYPE(list), DIMENSION(:), ALLOCATABLE :: group
  END TYPE Tops__Type
  TYPE(Tops__type), DIMENSION(:), ALLOCATABLE, SAVE, TARGET :: Res_char&
       &,Add_Char
  TYPE(Tops__type), DIMENSION(:), ALLOCATABLE, SAVE :: App_Char
  TYPE(Hash_Tops__Type), SAVE :: iKeys
  CHARACTER(len=max_char), DIMENSION(6), PARAMETER :: Oracs=(/'bonds',&
       &'imph ','don  ','acc  ','terma','backb'/)
  LOGICAL, ALLOCATABLE, SAVE :: ok_Residue(:)
  TYPE(Tops__Type), POINTER :: Resi

  CHARACTER(len=max_char), DIMENSION(:,:), ALLOCATABLE :: share
  TYPE(List), DIMENSION(:), ALLOCATABLE :: shareg
  TYPE(Resids), POINTER :: Ri
CONTAINS
  SUBROUTINE Tops_
    IMPLICIT NONE 
    INTEGER, SAVE :: i_L,n,m,p,o
    REAL(8) :: Time_begin, Time_end
    CALL Hash_Tops_

    IF(Solvent__Param % Build .AND. &
         &Solvent__Param % Added /=  0) &
         &CALL SecondarySeq__AddSlv(Solvent__Param % Added)

    ALLOCATE(Res_Char(SIZE(Residue)))
    ALLOCATE(ok_Residue(SIZE(Residue)))
    ok_residue=.FALSE.
    DO i_l=1,SIZE(Residue)
       Ri=>Residue(i_L)
       IF(TRIM(Ri % Type) == 'mass') THEN
          ok_Residue(i_L)=.TRUE.
          CYCLE
       END IF
       DO n=1,SIZE(Secondary)
          IF(.NOT. ALLOCATED(Secondary(n) % line)) CYCLE
          DO m=1,SIZE(Secondary (n) % line)
             IF(TRIM(Secondary (n) % line (m)) == TRIM(Ri% Type)) THEN
                ok_Residue(i_L)=.TRUE.
                EXIT
             END IF
          END DO
       END DO
       IF(.NOT. ALLOCATED(Patches)) CYCLE
       DO n=1,SIZE(Patches)
          IF(TRIM(patches(n) % pres) == TRIM(Ri% Type) &
               & .OR. TRIM(patches(n) % res) == TRIM(Ri% Type) ) THEN
             ok_Residue(i_L)=.TRUE.
             EXIT
          END IF
       END DO
    END DO
    CALL CPU_TIME(Time_Begin)
    DO i_l=1,SIZE(Residue)
       Ri=>Residue(i_L)
       IF(.NOT. ok_Residue(i_L)) CYCLE
       Res_Char(i_L) % FField=Ri % FField
       Res_Char(i_L) % Type=Ri % Type
       Res_Char(i_L) % Residue=Ri % Residue

       IF(Res_Char(i_L) % FField == 'CHARMM') THEN
          IF(Res_Char (i_L) % Type == 'mass') THEN
             CALL Mass_
          ELSE
             CALL Charmm_
          END IF
       ELSE IF(Residue(i_L) % FField == 'ORAC') THEN
          CALL Orac_
          CALL Charmm_
       ELSE
          STOP
       END IF
    END DO
    CALL Newresidues_
    CALL CPU_TIME(Time_End)
    WRITE(*,*) 'CPU Time = ',Time_End-Time_Begin,' in s '

  CONTAINS
    INCLUDE 'Tops_Orac.f90'
    INCLUDE 'Tops_Charmm.f90'
  END SUBROUTINE Tops_
  FUNCTION Tops__Store(Tops__, delete) RESULT(out)
    INTEGER, OPTIONAL :: delete
    TYPE(Tops__type), DIMENSION(:), ALLOCATABLE :: Tops__
    TYPE(Tops__type), DIMENSION(:), ALLOCATABLE, SAVE, TARGET :: store
    TYPE(Tops__type), DIMENSION(:), POINTER  :: out
    INTEGER :: o(2),n,m

    IF(PRESENT(delete)) THEN
       IF(ALLOCATED(Store)) DEALLOCATE(Store)
       out=>NULL()
       RETURN
    END IF
    IF(ALLOCATED(Store)) DEALLOCATE(Store)

    IF(ALLOCATED(Tops__)) THEN
       ALLOCATE(Store(SIZE(Tops__)))
       DO n=1,SIZE(Tops__)
          Store(n) % FField=Tops__(n) % FField
          Store(n) % Type=Tops__(n) % Type
          Store(n) % Residue=Tops__(n) % Residue
          IF(ALLOCATED(Tops__(n) % bonds)) THEN
             o=SHAPE(Tops__(n) % bonds)
             ALLOCATE(Store(n) % bonds(o(1), o(2)))
             Store(n) % Bonds=Tops__(n) % bonds
          END IF
          IF(ALLOCATED(Tops__(n) % imph)) THEN
             o=SHAPE(Tops__(n) % imph)
             ALLOCATE(Store(n) % imph(o(1), o(2)))
             Store(n) % imph=Tops__(n) % imph
          END IF
          IF(ALLOCATED(Tops__(n) % acc)) THEN
             o=SHAPE(Tops__(n) % acc)
             ALLOCATE(Store(n) % acc(o(1), o(2)))
             Store(n) % acc=Tops__(n) % acc
          END IF
          IF(ALLOCATED(Tops__(n) % don)) THEN
             o=SHAPE(Tops__(n) % don)
             ALLOCATE(Store(n) % don(o(1), o(2)))
             Store(n) % don=Tops__(n) % don
          END IF
          IF(ALLOCATED(Tops__(n) % acc_)) THEN
             o=SHAPE(Tops__(n) % acc_)
             ALLOCATE(Store(n) % acc_(o(1), o(2)))
             Store(n) % acc_=Tops__(n) % acc_
          END IF
          IF(ALLOCATED(Tops__(n) % don_)) THEN
             o=SHAPE(Tops__(n) % don_)
             ALLOCATE(Store(n) % don_(o(1), o(2)))
             Store(n) % don_=Tops__(n) % don_
          END IF
          IF(ALLOCATED(Tops__(n) % dele)) THEN
             o=SHAPE(Tops__(n) % dele)
             ALLOCATE(Store(n) % dele(o(1), o(2)))
             Store(n) % dele=Tops__(n) % dele
          END IF
          IF(ALLOCATED(Tops__(n) % ends)) THEN
             o=SHAPE(Tops__(n) % ends)
             ALLOCATE(Store(n) % ends(o(1), o(2)))
             Store(n) % ends=Tops__(n) % ends
          END IF
          IF(ALLOCATED(Tops__(n) % mass)) THEN
             o=SHAPE(Tops__(n) % mass)
             ALLOCATE(Store(n) % mass(o(1), o(2)))
             Store(n) % mass=Tops__(n) % mass
          END IF
          IF(ALLOCATED(Tops__(n) % group)) THEN
             ALLOCATE(Store(n) % group (SIZE(Tops__(n) % group)))
             DO m=1,SIZE(Tops__(n) % group)
                IF(ALLOCATED(Tops__(n) % group (m) % g)) THEN
                   ALLOCATE(Store(n) % group (m) % g (SIZE(Tops__(n) % group (m) % g)))
                   Store(n) % group(m) % g =Tops__(n) % group(m) % g
                END IF
             END DO
          END IF
       END DO
       out=>Store
    END IF
  END FUNCTION Tops__Store
  INCLUDE 'Newresidues_.f90'
END MODULE Tops
