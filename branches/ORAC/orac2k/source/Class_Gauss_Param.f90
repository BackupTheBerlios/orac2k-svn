!!$***********************************************************************
!!$   Time-stamp: <2005-02-19 18:13:07 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Mon Feb 14 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*
MODULE Class_Gauss_Param
  TYPE Gauss_Para
     LOGICAL :: initialized=.FALSE.
     CHARACTER(8) :: type
     REAL(8) :: n_plus,n_minus
     REAL(8) :: efact
     REAL(8) :: kt
  END TYPE Gauss_Para 
  TYPE(Gauss_Para), SAVE :: Param
CONTAINS
  FUNCTION Init(n_p,n_m,efact,t) RESULT(out)
!!$======================== DECLARATIONS ================================*
    IMPLICIT none
!!$----------------------------- ARGUMENTS ------------------------------*
    
    REAL(8), OPTIONAL :: n_p,n_m,t,efact
    LOGICAL :: out
    REAL(8), SAVE :: gascon = 8.3143D0

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    out=.FALSE.
    IF(.NOT. Param % Initialized) THEN
       IF(PRESENT(n_p) .AND. PRESENT(n_m) .AND. PRESENT(t) .AND.&
            & PRESENT(efact)) THEN 
          Param % initialized = .TRUE.
          Param % n_plus = n_p
          Param % n_minus = n_m
          Param % efact = efact
          Param % kt = gascon * t
          Param % Type = 'boltz'
          out=.TRUE.
       END IF
    END IF
  END FUNCTION Init

  FUNCTION Print() RESULT(out)
!!$======================== DECLARATIONS ================================*
    IMPLICIT none
!!$----------------------------- ARGUMENTS ------------------------------*
    TYPE(Gauss_Para) :: out
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    
    IF(Param % initialized) THEN
       out=Param
    ELSE
       out % initialized = .FALSE.
       out % n_plus=0.0
       out % n_minus=0.0
       out % kt=0.0
       out % Type='   '
    END IF
  END FUNCTION Print

  FUNCTION Found_Boltz(n,beta) RESULT(out)
!!$======================== DECLARATIONS ================================*
    IMPLICIT none
!!$----------------------------- ARGUMENTS ------------------------------*
    INTEGER :: n
    CHARACTER(7) :: beta(*)
    LOGICAL :: out
!!$------------------------- LOCAL VARIABLES ----------------------------*
    INTEGER :: i
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    out=.FALSE.
    DO i=1,n
       IF(beta(i) == 'boltz') out=.TRUE.
    END DO
  END FUNCTION Found_Boltz

END MODULE Class_Gauss_Param
