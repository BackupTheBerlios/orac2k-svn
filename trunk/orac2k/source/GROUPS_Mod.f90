MODULE GROUPS_Mod

!!$***********************************************************************
!!$   Time-stamp: <2006-02-15 16:08:46 marchi2>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Wed Feb  8 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************
!!$---- This module is part of the program ORAC ----*

  USE Xerror_Mod
  USE INPUT_Mod, ONLY: Read_String, Parser, Parse_Numbers,err_open&
       &,err_end,err_unr,err_fnf,err_args
  LOGICAL, SAVE :: groups=.FALSE.
  INTEGER :: n_molecs
  PARAMETER(n_molecs=100)
  TYPE group
     INTEGER :: n
     CHARACTER(20) :: mode
     CHARACTER(20) :: type
     INTEGER, DIMENSION(:), POINTER :: index
  END TYPE group
  TYPE(group), DIMENSION(:), ALLOCATABLE, SAVE :: molecs

CONTAINS
  SUBROUTINE Init(ntap,nbun,grppt,mres)

!!$======================== DECLARATIONS ================================*

    IMPLICIT NONE 

!!$----------------------------- ARGUMENTS ------------------------------*

    INTEGER :: ntap,nbun,mres(:,:),grppt(:,:)
    CHARACTER(80) :: errmsg

!!$----------------------- VARIABLES IN COMMON --------------------------*
    
    INTEGER :: i,k,l,n,m,count
    INTEGER, DIMENSION(:), ALLOCATABLE :: index

!!$------------------------- LOCAL VARIABLES ----------------------------*

    ALLOCATE(index(ntap))
    i=1
    DO WHILE(molecs(i)%n /= 0)
       SELECT CASE(TRIM(molecs(i)%mode))
       CASE('atoms')
          DO n=1,molecs(i)%n
             IF(molecs(i)%index(n) < 1 .OR. molecs(i)%index(n) >&
                  & ntap) THEN
                errmsg='Group list contains atoms not in the '//&
                     &'allowed range.'
                CALL abort_now(errmsg)
             END IF
          END DO
       CASE('residue')
          count=0
          DO n=1,molecs(i)%n
             IF(molecs(i)%index(n) < 1 .OR. molecs(i)%index(n) >&
                  & nbun) THEN
                errmsg='Group list contains residues not in the '//&
                     &'allowed range.'
                CALL abort_now(errmsg)
             ELSE
                m=molecs(i)%index(n)
                DO k=mres(1,m),mres(2,m)
                   DO l=grppt(1,k),grppt(2,k)
                      count=count+1
                      index(count)=l
                   END DO
                END DO
             END IF
          END DO
          DEALLOCATE(molecs(i)%index)
          ALLOCATE(molecs(i)%index(count))
          molecs(i)%index=index
          molecs(i)%n=count
       END SELECT       
       i=i+1
    END DO
    DEALLOCATE(index)

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

  END SUBROUTINE Init
  INCLUDE 'GROUPS_Readit.f90'
END MODULE GROUPS_Mod
