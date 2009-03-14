!!$***********************************************************************
!!$   Time-stamp: <2005-03-03 21:03:05 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Feb 18 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program ORAC ----*
MODULE Class_Connect
  TYPE Connect
     INTEGER, POINTER :: m=>NULL()
     INTEGER, DIMENSION(:), POINTER :: cnt=>NULL()
  END TYPE Connect
  INTEGER, PRIVATE, SAVE :: nato, contacts
CONTAINS
  FUNCTION Init(concta,m,n) RESULT (out)
!!$======================== DECLARATIONS ================================*
    IMPLICIT none
!!$----------------------------- ARGUMENTS ------------------------------*
    INTEGER :: m,n
    INTEGER, TARGET :: concta(m,*)
    TYPE(Connect), DIMENSION(:), POINTER :: out
!!$------------------------- LOCAL VARIABLES ----------------------------*
    INTEGER :: i,o,p(1)
    LOGICAL, SAVE :: done=.FALSE.
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*
    nato = n
    p = MAXLOC(concta(1:n,1))
    contacts=p(1)
    ALLOCATE(out(n))
    DO i=1,n
       o=concta(i,1)
       out(i) % m => concta(i,1)
       out(i) % cnt => concta(i,2:o+1)
    END DO
    
  END FUNCTION Init
  SUBROUTINE Destroy(in)
!!$======================== DECLARATIONS ================================*
    IMPLICIT none
!!$----------------------------- ARGUMENTS ------------------------------*
    TYPE(Connect), DIMENSION(:), POINTER :: in
!!$------------------------- LOCAL VARIABLES ----------------------------*
    INTEGER :: i
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    DO i=1,nato
       NULLIFY(in(i) % m)
       NULLIFY(in(i) % cnt)
    END DO
    DEALLOCATE(in)
  END SUBROUTINE Destroy
  FUNCTION Print() RESULT(out)
    IMPLICIT none
    INTEGER :: out(2)
    out(1)=nato
    out(2)=contacts
  END FUNCTION Print
END MODULE Class_Connect
