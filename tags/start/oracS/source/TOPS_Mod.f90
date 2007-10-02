MODULE TOPS_Mod

!!$***********************************************************************
!!$   Time-stamp: <2006-12-15 13:57:56 marchi>                           *
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
  

  USE TYPES
  USE PARAMETERS_GLOBALS
  USE XERROR, ONLY: abort_now
  USE MYPARSE_Mod, ONLY: MY_parse=>parse
  USE STRINGS_MOD, ONLY: Resid, MY_Fam, ST_Concat

  CHARACTER(len=max_pars), DIMENSION(:), ALLOCATABLE, PRIVATE ::&
       & strngs
  PRIVATE Transform_CHARMM,Transform_ORAC, Transform_MASS
  TYPE(KHASH), DIMENSION(8), SAVE :: Hash_CHARMM,Hash_ORAC
  TYPE(Unit_Char), POINTER :: Resi

CONTAINS
  SUBROUTINE Transform(Residue)
    IMPLICIT NONE 
    TYPE(Resid), DIMENSION(:) :: Residue
    INTEGER :: i,n,p

    ALLOCATE(Res_Char(SIZE(Residue)))
    DO i=1,SIZE(Residue)
       Resi=>Res_Char(i)
       Resi%FField=Residue(i)%FField
       Resi%Type=Residue(i)%Type
       Resi%Residue=Residue(i)%Residue
       IF(Residue(i)%FField == 'CHARMM') THEN
          IF(Residue(i)%Type == 'mass') THEN
             CALL Transform_MASS(Residue(i))
          ELSE
             CALL Transform_CHARMM(Residue(i))
          END IF
       ELSE IF(Residue(i)%FField == 'ORAC') THEN
          CALL Transform_ORAC(Residue(i))
       END IF
    END DO
  END SUBROUTINE Transform
  FUNCTION IHash(types,hash)
    IMPLICIT NONE 
    INTEGER :: IHash
    CHARACTER(len=*) :: types
    TYPE(KHash), DIMENSION(:) :: hash
    CHARACTER(len=400) :: errmsg
    INTEGER :: n

    IHash=-1
    DO n=1,SIZE(hash)
       IF(TRIM(hash(n)%type) == TRIM(types)) THEN
          IHash=n
       END IF
    END DO
    IF(IHash == -1) THEN
       errmsg='ORAC internal: Type '&
            &//TRIM(types)//' not found in hash Tables'
       CALL abort_now(errmsg)
    END IF
  END FUNCTION IHash

  INCLUDE 'TOPS_Orac.f90'
  INCLUDE 'TOPS_Charmm.f90'
END MODULE TOPS_Mod
