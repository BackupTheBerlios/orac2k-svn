!!$***********************************************************************
!!$   Time-stamp: <2005-03-06 22:24:37 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Mar  3 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*
MODULE Module_Neighbors
  IMPLICIT none
  TYPE Neighbors
     INTEGER :: no
     INTEGER, DIMENSION (:), POINTER :: nb
  END TYPE Neighbors
  TYPE(Neighbors), DIMENSION (:), POINTER, SAVE :: neigha
  TYPE(Neighbors), DIMENSION (:), POINTER, SAVE :: neighb
  TYPE(Neighbors), DIMENSION (:), POINTER, SAVE :: neighc
CONTAINS
  SUBROUTINE Start(neigh,nmol)
    IMPLICIT none
    TYPE(Neighbors), DIMENSION (:), POINTER :: neigh
    INTEGER :: nmol

    INTEGER :: i

    ALLOCATE(neigh(nmol))
    neigh(:) % no = 0
  END SUBROUTINE Start
  SUBROUTINE Delete(neigh)
    IMPLICIT none
    INTEGER :: i
    TYPE(Neighbors), DIMENSION (:), POINTER :: neigh

    IF(ASSOCIATED(neigh)) THEN
       IF(SIZE(neigh) > 0) THEN
          DO i=1,SIZE(neigh)
             IF(ASSOCIATED(neigh(i) % nb)) THEN
                DEALLOCATE(neigh(i) % nb)
             END IF
          END DO
          DEALLOCATE(neigh)
       END IF
    END IF
  END SUBROUTINE Delete
END MODULE Module_Neighbors
