!!$***********************************************************************
!!$   Time-stamp: <2005-12-08 17:24:47 marchi>                           *
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
    DO i=1,nmol
       neigh(i) % no = 0
       NULLIFY(neigh(i) % nb)
    END DO
  END SUBROUTINE Start
  SUBROUTINE Delete(neigh)
    IMPLICIT none
    INTEGER :: i,ierr
    TYPE(Neighbors), DIMENSION (:), POINTER :: neigh
    
    IF(ASSOCIATED(neigh)) THEN
       DO i=1,SIZE(neigh)
          IF(ASSOCIATED(neigh(i) % nb)) THEN
             DEALLOCATE(neigh(i) % nb, STAT=ierr)
             IF(ierr /= 0) THEN
                WRITE(*,*) 'ierr not zero',ierr
             END IF
             NULLIFY(neigh(i) % nb)
          END IF
       END DO

       DEALLOCATE(neigh, STAT=ierr)
       NULLIFY(neigh)
    END IF
  END SUBROUTINE Delete
END MODULE Module_Neighbors

