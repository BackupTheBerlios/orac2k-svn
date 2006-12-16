!!$***********************************************************************
!!$   Time-stamp: <2006-12-13 17:34:37 marchi>                           *
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
MODULE Linked_Int4_D
  USE CONSTANTS, ONLY: max_char
  PRIVATE
  PUBLIC :: Start, Add, Extract, Cleanup, NODE_INT4_D
  TYPE :: NODE_Int4_D
     INTEGER, DIMENSION(:), POINTER  :: labels=>NULL()
     TYPE(NODE_Int4_D), POINTER :: next=>NULL()
  END TYPE NODE_INT4_D
  TYPE(NODE_Int4_D), POINTER :: root
  TYPE(NODE_Int4_D), POINTER :: current
  INTEGER, SAVE :: count_a=0
  INTERFACE Extract
     MODULE PROCEDURE Extract_ST
     MODULE PROCEDURE Extract_DYN
  END INTERFACE
CONTAINS
  SUBROUTINE Start()
    IMPLICIT none
    count_a=0
    IF(ASSOCIATED(root)) THEN
       CALL Cleanup()
    END IF
    ALLOCATE(root)
    NULLIFY(root%next)
    current=>root
  END SUBROUTINE Start
  SUBROUTINE Add(labels, count_out)
    IMPLICIT none
    INTEGER, DIMENSION(:), INTENT(IN) :: labels
    INTEGER, INTENT(OUT) :: count_out
    INTEGER :: kk
    TYPE(NODE_Int4_D), POINTER :: new_node

    
    count_a=count_a+1

    ALLOCATE(new_node)
    ALLOCATE(new_node % labels(SIZE(labels)))
    new_node % labels = labels
    NULLIFY(new_node % next)

    current % next => new_node
    current => current % next

    count_out=count_a
  END SUBROUTINE Add

  SUBROUTINE Extract_ST(labs,end_of_list,init)
    IMPLICIT none
    INTEGER, INTENT(IN), OPTIONAL :: init
    INTEGER, DIMENSION(:), INTENT(OUT) :: labs
    LOGICAL, INTENT(OUT)  :: end_of_list
    IF(PRESENT(init)) THEN
       IF(init == 0) current=>root
    END IF
    end_of_list=.NOT. ASSOCIATED(current % next )
    IF(end_of_list) RETURN
    current=>current % next
    labs=current % labels
  END SUBROUTINE Extract_ST
  SUBROUTINE Extract_DYN(labs,end_of_list,Dyn,init)
    IMPLICIT none
    INTEGER, INTENT(IN), OPTIONAL :: init
    LOGICAL, INTENT(IN)  :: Dyn
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: labs
    LOGICAL, INTENT(OUT)  :: end_of_list
    IF(PRESENT(init)) THEN
       IF(init == 0) current=>root
    END IF
    end_of_list=.NOT. ASSOCIATED(current % next )
    IF(end_of_list) RETURN
    current=>current % next
    IF(ALLOCATED(labs)) DEALLOCATE(labs)
    ALLOCATE(labs(SIZE(current % labels)))
    labs=current % labels
  END SUBROUTINE Extract_DYN

  SUBROUTINE Cleanup()
    CALL Remove(root)
  END SUBROUTINE Cleanup
  RECURSIVE SUBROUTINE Remove(old_node)
    TYPE(NODE_Int4_D), POINTER :: old_node
    INTEGER, SAVE :: uppa=0
    
    IF(.NOT. ASSOCIATED(old_node) ) RETURN
    CALL Remove(old_node % next)
    IF(ASSOCIATED(old_node % labels)) DEALLOCATE(old_node % labels)
    DEALLOCATE(old_node)
  END SUBROUTINE Remove
END MODULE Linked_Int4_D
MODULE Linked_Real8_D
  USE CONSTANTS, ONLY: max_char
  PRIVATE
  PUBLIC :: Start, Add, Extract, Cleanup, NODE_INT4_D
  TYPE :: NODE_Real8_D
     REAL(8), DIMENSION(:), POINTER  :: labels=>NULL()
     TYPE(NODE_Real8_D), POINTER :: next=>NULL()
  END TYPE NODE_REAL8_D
  TYPE(NODE_Real8_D), POINTER :: root
  TYPE(NODE_Real8_D), POINTER :: current
  INTEGER, SAVE :: count_a=0
  INTERFACE Extract
     MODULE PROCEDURE Extract_ST
     MODULE PROCEDURE Extract_DYN
  END INTERFACE
CONTAINS
  SUBROUTINE Start()
    IMPLICIT none
    count_a=0
    IF(ASSOCIATED(root)) THEN
       CALL Cleanup()
    END IF
    ALLOCATE(root)
    NULLIFY(root%next)
    current=>root
  END SUBROUTINE Start
  SUBROUTINE Add(labels, count_out)
    IMPLICIT none
    REAL(8), DIMENSION(:), INTENT(IN) :: labels
    INTEGER, INTENT(OUT) :: count_out
    INTEGER :: kk
    TYPE(NODE_Real8_D), POINTER :: new_node

    
    count_a=count_a+1

    ALLOCATE(new_node)
    ALLOCATE(new_node % labels(SIZE(labels)))
    new_node % labels = labels
    NULLIFY(new_node % next)

    current % next => new_node
    current => current % next

    count_out=count_a
  END SUBROUTINE Add

  SUBROUTINE Extract_ST(labs,end_of_list,init)
    IMPLICIT none
    INTEGER, INTENT(IN), OPTIONAL :: init
    REAL(8), DIMENSION(:), INTENT(OUT) :: labs
    LOGICAL, INTENT(OUT)  :: end_of_list
    IF(PRESENT(init)) THEN
       IF(init == 0) current=>root
    END IF
    end_of_list=.NOT. ASSOCIATED(current % next )
    IF(end_of_list) RETURN
    current=>current % next
    labs=current % labels
  END SUBROUTINE Extract_ST
  SUBROUTINE Extract_DYN(labs,end_of_list,Dyn,init)
    IMPLICIT none
    INTEGER, INTENT(IN), OPTIONAL :: init
    LOGICAL, INTENT(IN)  :: Dyn
    REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: labs
    LOGICAL, INTENT(OUT)  :: end_of_list
    IF(PRESENT(init)) THEN
       IF(init == 0) current=>root
    END IF
    end_of_list=.NOT. ASSOCIATED(current % next )
    IF(end_of_list) RETURN
    current=>current % next
    IF(ALLOCATED(labs)) DEALLOCATE(labs)
    ALLOCATE(labs(SIZE(current % labels)))
    labs=current % labels
  END SUBROUTINE Extract_DYN

  SUBROUTINE Cleanup()
    CALL Remove(root)
  END SUBROUTINE Cleanup
  RECURSIVE SUBROUTINE Remove(old_node)
    TYPE(NODE_Real8_D), POINTER :: old_node
    INTEGER, SAVE :: uppa=0
    
    IF(.NOT. ASSOCIATED(old_node) ) RETURN
    CALL Remove(old_node % next)
    IF(ASSOCIATED(old_node % labels)) DEALLOCATE(old_node % labels)
    DEALLOCATE(old_node)
  END SUBROUTINE Remove
END MODULE Linked_Real8_D
