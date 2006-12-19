MODULE TYPES_Utils

!!$***********************************************************************
!!$   Time-stamp: <2006-12-19 09:35:05 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Mon Dec  4 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*

  USE TYPES
  IMPLICIT none
  INTERFACE Move
     MODULE PROCEDURE Move_Real_1
     MODULE PROCEDURE Move_Real_2
  END INTERFACE
  INTERFACE copy
     MODULE PROCEDURE Copy_UnitChar0
     MODULE PROCEDURE Copy_UnitChar1
     MODULE PROCEDURE Copy_UnitChar2
  END INTERFACE
CONTAINS
  SUBROUTINE Move_Real_1(a,b)
    REAL(8), DIMENSION(:), ALLOCATABLE :: a,b
    INTEGER :: m_Size
    m_Size=SIZE(a)
    ALLOCATE(b(m_Size))
    b=a
    DEALLOCATE(a)
    
  END SUBROUTINE Move_Real_1
  SUBROUTINE Move_Real_2(a,b)
    REAL(8), DIMENSION(:,:), ALLOCATABLE  :: a,b
    INTEGER :: m_Size,n_Size
    m_Size=SIZE(a,1)
    n_Size=SIZE(a,2)
    ALLOCATE(b(m_Size,n_size))
    b=a
    DEALLOCATE(a)    
  END SUBROUTINE Move_Real_2
  SUBROUTINE Copy_UnitChar1(a,b,i,o)
    TYPE(Unit_Char), DIMENSION (:), ALLOCATABLE  :: a,b
    INTEGER :: i,o
    INTEGER :: m_size,j,n,m

    CALL Alloc(a (i) % bonds, b (o) % bonds)
    CALL Alloc(a (i) % imph , b (o) % imph )
    CALL Alloc(a (i) % acc, b (o) % acc)
    CALL Alloc(a (i) % acc_, b (o) % acc_)
    CALL Alloc(a (i) % don , b (o) % don )
    CALL Alloc(a (i) % don_, b (o) % don_)
    CALL Alloc(a (i) % dele, b (o) % dele)
    CALL Alloc(a (i) % ends, b (o) % ends)
    IF(ALLOCATED(a (i) % group)) THEN
       IF(ALLOCATED(b (o) % group)) DEALLOCATE(b(o) % group)
       ALLOCATE(b (o ) % group (SIZE(a (i) % group)))
       DO m=1,SIZE(a (i) % group)
          ALLOCATE(b (o) % group (m) % g (SIZE(a (i) % group (m) % g)))
       END DO
    END IF
    IF(ALLOCATED(a (i) % mass)) THEN
       IF(ALLOCATED(b (o) % mass)) DEALLOCATE( b(o) % mass)
       ALLOCATE(b (o) % mass (2,SIZE(a(i) % mass)))
    END IF
    b(o)=a(i)
  CONTAINS
    SUBROUTINE Alloc(tpg_o,tpg_n)
      CHARACTER(len=max_char), DIMENSION (:,:),  ALLOCATABLE :: tpg_o, tpg_n
      INTEGER :: p1,p2
      
      IF(ALLOCATED(Tpg_o)) THEN
         p1=SIZE(Tpg_o,1)
         p2=SIZE(Tpg_o,2)
         IF(ALLOCATED(Tpg_n)) DEALLOCATE(Tpg_n)
         ALLOCATE(Tpg_n(p1,p2))
      END IF
    END SUBROUTINE Alloc
  END SUBROUTINE Copy_UnitChar1
  SUBROUTINE Copy_UnitChar2(a,b,i)
    TYPE(Unit_Char), DIMENSION (:) :: a
    TYPE(Unit_Char) :: b
    INTEGER :: i
    INTEGER :: m_size,j,n,m
    CALL Alloc(a (i) % bonds, b % bonds)
    CALL Alloc(a (i) % imph , b % imph )
    CALL Alloc(a (i) % acc,  b % acc)
    CALL Alloc(a (i) % acc_, b % acc_)
    CALL Alloc(a (i) % don , b % don )
    CALL Alloc(a (i) % don_, b % don_)
    CALL Alloc(a (i) % dele, b % dele)
    CALL Alloc(a (i) % ends, b % ends)
    IF(ALLOCATED(a (i) % group)) THEN
       ALLOCATE(b % group (SIZE(a (i) % group)))
       DO m=1,SIZE(a (i) % group)
          ALLOCATE(b % group (m) % g (SIZE(a (i) % group (m) % g)))
       END DO
    END IF
    IF(ALLOCATED(a (i) % mass)) THEN
       ALLOCATE(b % mass (2,SIZE(a(i) % mass)))
    END IF
    b=a(i)
  CONTAINS
    SUBROUTINE Alloc(tpg_o,tpg_n)
      CHARACTER(len=max_char), DIMENSION (:,:),  ALLOCATABLE :: tpg_o, tpg_n
      INTEGER :: p1,p2
      
      IF(ALLOCATED(Tpg_o)) THEN
         p1=SIZE(Tpg_o,1)
         p2=SIZE(Tpg_o,2)
         IF(ALLOCATED(Tpg_n)) DEALLOCATE(Tpg_n)
         ALLOCATE(Tpg_n(p1,p2))
      END IF
    END SUBROUTINE Alloc
  END SUBROUTINE Copy_UnitChar2
  SUBROUTINE Copy_UnitChar0(a,b)
    TYPE(Unit_Char), DIMENSION (:), ALLOCATABLE  :: a,b
    INTEGER :: i,o
    INTEGER :: m_size,j,n,m

    ALLOCATE(b(SIZE(a)))
    DO i=1,SIZE(a)
       CALL Alloc(a (i) % bonds, b (i) % bonds)
       CALL Alloc(a (i) % imph , b (i) % imph )
       CALL Alloc(a (i) % acc, b (i) % acc)
       CALL Alloc(a (i) % acc_, b (i) % acc_)
       CALL Alloc(a (i) % don , b (i) % don )
       CALL Alloc(a (i) % don_, b (i) % don_)
       CALL Alloc(a (i) % dele, b (i) % dele)
       CALL Alloc(a (i) % ends, b (i) % ends)
       IF(ALLOCATED(a (i) % group)) THEN
          ALLOCATE(b (i) % group (SIZE(a (i) % group)))
          DO m=1,SIZE(a (i) % group)
             ALLOCATE(b (i) % group (m) % g (SIZE(a (i) % group (m) % g)))
          END DO
       END IF
       IF(ALLOCATED(a (i) % mass)) THEN
          ALLOCATE(b (i) % mass (2,SIZE(a(i) % mass,2)))
       END IF
    END DO
    b=a
  CONTAINS
    SUBROUTINE Alloc(tpg_o,tpg_n)
      CHARACTER(len=max_char), DIMENSION (:,:),  ALLOCATABLE :: tpg_o, tpg_n
      INTEGER :: p1,p2
      
      IF(ALLOCATED(Tpg_o)) THEN
         p1=SIZE(Tpg_o,1)
         p2=SIZE(Tpg_o,2)
         IF(ALLOCATED(Tpg_n)) DEALLOCATE(Tpg_n)
         ALLOCATE(Tpg_n(p1,p2))
      END IF
    END SUBROUTINE Alloc
  END SUBROUTINE Copy_UnitChar0
END MODULE TYPES_Utils
