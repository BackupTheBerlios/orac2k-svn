!!$***********************************************************************
!!$   Time-stamp: <2005-03-03 12:46:25 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Tue Mar  1 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*
MODULE Module_Node
  IMPLICIT none

  PUBLIC :: Start, Add, Extract, Cleanup
  INTEGER, PARAMETER, PUBLIC :: max_char = 7

  TYPE, PRIVATE :: NODE_PAR
     CHARACTER(len=max_char), DIMENSION(:), POINTER :: labels
     REAL(8), DIMENSION(:), POINTER :: pars
     TYPE(NODE_PAR), POINTER :: next
  END TYPE NODE_PAR
  TYPE(NODE_PAR), POINTER, PRIVATE, SAVE :: root
  TYPE(NODE_PAR), POINTER, PRIVATE :: current
CONTAINS
  SUBROUTINE Start()
    IMPLICIT none
    IF(ASSOCIATED(root)) THEN
       CALL Cleanup()
    END IF
    ALLOCATE(current)
    NULLIFY(current%next)
    root=>current
  END SUBROUTINE Start
  SUBROUTINE Add(labels,n_l,pars,n_p)
    IMPLICIT none
    REAL(8), DIMENSION(:), INTENT(IN) :: pars
    CHARACTER(len=max_char), DIMENSION (:), INTENT(IN) :: labels
    INTEGER, INTENT(IN) :: n_l,n_p

    ALLOCATE(current % labels (n_l))
    ALLOCATE(current % pars (n_p))
    current % labels = labels(1:n_l)
    current % pars = pars(1:n_p)
    ALLOCATE(current % next)
    current=>current % next
    NULLIFY(current % next)
  END SUBROUTINE Add
  SUBROUTINE Extract(labels,pars,No_Param)
    IMPLICIT none
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: pars
    CHARACTER(len=max_char), DIMENSION (:,:), ALLOCATABLE :: labels
    INTEGER, INTENT(OUT) :: No_Param
    INTEGER :: pars_max=-1, labels_max=-1
    INTEGER :: count,n_l,n_p

    current=>root
    count=0
    DO WHILE(ASSOCIATED(current % next))
       IF(pars_max < SIZE(current % pars)) pars_max=SIZE(current % pars)
       IF(labels_max < SIZE(current % labels)) labels_max=SIZE(current % labels)
       count=count+1
       current=>current % next
    END DO
    ALLOCATE(pars(pars_max,count))
    ALLOCATE(labels(labels_max,count))
    No_Param=count
    current=>root
    count=0
    DO WHILE(ASSOCIATED(current % next))
       count=count+1
       n_p=SIZE(current % pars)
       n_l=SIZE(current % labels)
       labels(1:n_l,count)=current % labels(1:n_l)
       pars(1:n_p,count)=current % pars(1:n_p)
       current=>current % next
    END DO
  END SUBROUTINE Extract
  SUBROUTINE Cleanup()
    IMPLICIT none
    TYPE(NODE_PAR), POINTER :: dummy
    INTEGER :: count
    
    current=>root
    DO WHILE(ASSOCIATED(current % next))
       dummy=>current % next
       DEALLOCATE(current)
       current=>dummy
    END DO
    NULLIFY(dummy, root)
  END SUBROUTINE Cleanup
END MODULE Module_Node
PROGRAM main
  USE Module_Node
  IMPLICIT none
  CHARACTER(len=max_char), DIMENSION (3) :: labels
  REAL(8), DIMENSION (4) :: pars
  INTEGER :: n_l,n_p
  CHARACTER(len=max_char), DIMENSION (:,:), ALLOCATABLE :: lbend
  REAL(8), DIMENSION (:,:), ALLOCATABLE :: pbend
  INTEGER :: n=100,i,nbend

  CALL Start()
  n_l=3
  n_p=2
  DO i=1,n
     IF(MOD(i,2) == 0) THEN
        labels(1)='ca'
        labels(2)='o'
        labels(3)='cb'
        pars(1)=12.34D0
        pars(2)=2.4D0
     ELSE
        labels(1)='cb'
        labels(2)='o2'
        labels(3)='p2'
        pars(1)=10.00D0
        pars(2)=1.89D0
     END IF
     CALL Add(labels,n_l,pars,n_p)
  END DO
  WRITE(*,*) n_l
  CALL Extract(lbend,pbend,nbend)
  WRITE(*,*) nbend,lbend(1:3,1:3)
  CALL Cleanup()
  DEALLOCATE(lbend,pbend)
END PROGRAM main

