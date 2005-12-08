MODULE Neighbors_Mod

!!$***********************************************************************
!!$   Time-stamp: <2005-12-08 14:03:32 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Feb  4 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

  TYPE Neighbors
     INTEGER :: m
     INTEGER, DIMENSION(:), POINTER :: neigh
  END TYPE Neighbors
  TYPE(Neighbors), DIMENSION (:), TARGET, ALLOCATABLE, SAVE :: nn

!!$---- This subroutine is part of the program ORAC ----*

CONTAINS
  SUBROUTINE Neighbors_init(nnlpp,ngrp,index)

!!$---- This subroutine is part of the program ORAC ----*

    IMPLICIT none

    INTEGER, TARGET :: nnlpp(*)
    INTEGER :: ngrp
    INTEGER, OPTIONAL :: index(*)

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: i, m, n, count
    LOGICAL :: indx

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*
    
    indx=.FALSE.
    IF(PRESENT(index)) indx=.TRUE.

    ALLOCATE(nn(ngrp))

    count=0
    n=0
    DO i=1,ngrp
       IF(indx) THEN
          IF(index(i) == 2) THEN
             n=n+1
             m=nnlpp(1+count)
             nn(n)%m=m
             nn(n)%Neigh=>nnlpp(1+count+1:1+count+m)
             count=count+m+1
          ENDIF
       ELSE
          n=n+1
          m=nnlpp(1+count)
          nn(n)%m=m
          nn(n)%Neigh=>nnlpp(1+count+1:1+count+m)
          count=count+m+1
       END IF
    END DO
  END SUBROUTINE Neighbors_init
END MODULE Neighbors_Mod
