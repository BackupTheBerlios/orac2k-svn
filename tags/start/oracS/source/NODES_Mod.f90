!!$***********************************************************************
!!$   Time-stamp: <2006-11-24 12:52:44 marchi>                           *
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
MODULE NODES_Mod
  USE CONSTANTS, ONLY: max_char
  PUBLIC :: Start, Add, Extract, Cleanup
  TYPE, PRIVATE :: NODE_PAR
     CHARACTER(len=max_char)  :: labels
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
  SUBROUTINE Add(labels)
    IMPLICIT none
    CHARACTER(len=*), INTENT(IN) :: labels

    current % labels = labels
    ALLOCATE(current % next)
    current=>current % next
    NULLIFY(current % next)
  END SUBROUTINE Add
  SUBROUTINE Extract(labs,No_Param)
    IMPLICIT none
    CHARACTER(len=*), DIMENSION (:), ALLOCATABLE :: labs
    INTEGER, INTENT(OUT) :: No_Param
    INTEGER :: labels_max=-1
    INTEGER :: count

    current=>root
    count=0
    DO WHILE(ASSOCIATED(current % next))
       count=count+1
       current=>current % next
    END DO
    ALLOCATE(labs(count))
    No_Param=count
    current=>root
    count=0
    DO WHILE(ASSOCIATED(current % next))
       count=count+1
       labs(count)=current % labels
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
END MODULE NODES_Mod
