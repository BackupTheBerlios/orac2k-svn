MODULE Tet
  IMPLICIT NONE 
  PRIVATE
  PUBLIC Tet_, Ulo
  TYPE :: Ulo
     INTEGER, ALLOCATABLE :: p(:)
     REAL(8) :: g
  END TYPE Ulo
  TYPE(Ulo), DIMENSION(:), ALLOCATABLE, SAVE, TARGET :: a
  INTEGER, PARAMETER :: ntet=10
CONTAINS
  FUNCTION tet_(n) RESULT(out)
    INTEGER :: n
    TYPE(Ulo), POINTER :: out
    LOGICAL, SAVE :: first_Time=.TRUE.

    out=>NULL()
    IF(First_Time) THEN
       ALLOCATE(a(ntet))
       First_time=.FALSE.
    END IF
    IF(ALLOCATED(a(n) % p)) DEALLOCATE(a(n) % p)
    ALLOCATE(a(n) % p(4))
    a(n) % p = 3
    a(n) % g = 4.0D0
    out=>a(n)
  END FUNCTION tet_

END MODULE Tet
PROGRAM temp

!!$***********************************************************************
!!$   Time-stamp: <2002-09-27 10:49:30 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Tue Jan  9 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program  ----*


!!$======================== DECLARATIONS ================================*

  USE Tet
  IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

!!$----------------------- VARIABLES IN COMMON --------------------------*

  TYPE(Ulo)  :: a
  TYPE(Ulo), POINTER  :: b

!!$------------------------- LOCAL VARIABLES ----------------------------*


!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

  
  b=>Tet_(2)
  IF(ASSOCIATED(b)) ALLOCATE(a % p(SIZE(b % p)))
  a=b
  WRITE(*,*) SIZE(a % p)
  WRITE(*,*) a % p

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END PROGRAM temp
