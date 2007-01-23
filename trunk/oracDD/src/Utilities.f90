MODULE Utilities

!!$***********************************************************************
!!$   Time-stamp: <2007-01-23 11:12:11 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Tue Jan 23 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program oracDD ----*

  IMPLICIT none
  PRIVATE
  PUBLIC MyShape

  INTERFACE MYSHAPE
     MODULE PROCEDURE Shapes_I4_1
     MODULE PROCEDURE Shapes_I4_2
     MODULE PROCEDURE Shapes_I4_3
     MODULE PROCEDURE Shapes_R8_1
     MODULE PROCEDURE Shapes_R8_2
     MODULE PROCEDURE Shapes_R8_3
  END INTERFACE
  INTERFACE MYSIZE
     MODULE PROCEDURE MySize_I4_1
     MODULE PROCEDURE MySize_I4_2
     MODULE PROCEDURE MySize_I4_3
     MODULE PROCEDURE MySize_R8_1
     MODULE PROCEDURE MySize_R8_2
     MODULE PROCEDURE MySize_R8_3
  END INTERFACE
CONTAINS
  FUNCTION Shapes_I4_1(arr) RESULT(out)
    INTEGER, ALLOCATABLE :: arr(:)
    INTEGER :: out(1)
    
    IF(ALLOCATED(arr)) THEN
       out=SHAPE(arr)
    ELSE
       out(1)=0
    END IF
  END FUNCTION Shapes_I4_1
  FUNCTION Shapes_I4_2(arr) RESULT(out)
    INTEGER, ALLOCATABLE :: arr(:,:)
    INTEGER :: out(2)
    
    IF(ALLOCATED(arr)) THEN
       out=SHAPE(arr)
    ELSE
       out(:)=0
    END IF
  END FUNCTION Shapes_I4_2
  FUNCTION Shapes_I4_3(arr) RESULT(out)
    INTEGER, ALLOCATABLE :: arr(:,:,:)
    INTEGER :: out(3)
    
    IF(ALLOCATED(arr)) THEN
       out=SHAPE(arr)
    ELSE
       out(:)=0
    END IF
  END FUNCTION Shapes_I4_3
  FUNCTION Shapes_R8_1(arr) RESULT(out)
    REAL(8), ALLOCATABLE :: arr(:)
    INTEGER :: out(1)
    
    IF(ALLOCATED(arr)) THEN
       out=SHAPE(arr)
    ELSE
       out(1)=0
    END IF
  END FUNCTION Shapes_R8_1
  FUNCTION Shapes_R8_2(arr) RESULT(out)
    REAL(8), ALLOCATABLE :: arr(:,:)
    INTEGER :: out(2)
    
    IF(ALLOCATED(arr)) THEN
       out=SHAPE(arr)
    ELSE
       out(:)=0
    END IF
  END FUNCTION Shapes_R8_2
  FUNCTION Shapes_R8_3(arr) RESULT(out)
    REAL(8), ALLOCATABLE :: arr(:,:,:)
    INTEGER :: out(3)
    
    IF(ALLOCATED(arr)) THEN
       out=SHAPE(arr)
    ELSE
       out(:)=0
    END IF
  END FUNCTION Shapes_R8_3
  FUNCTION MySize_I4_1(arr) RESULT(out)
    INTEGER, ALLOCATABLE :: arr(:)
    INTEGER :: out
    IF(ALLOCATED(arr)) THEN
       out=SIZE(arr)
    ELSE
       out=0
    END IF
  END FUNCTION MySize_I4_1
  FUNCTION MySize_I4_2(arr,ndim) RESULT(out)
    INTEGER, ALLOCATABLE :: arr(:,:)
    INTEGER, OPTIONAL :: ndim
    INTEGER :: out
    INTEGER :: n
    n=0
    IF(PRESENT(ndim)) n=ndim
    IF(ALLOCATED(arr)) THEN
       IF(n /= 0) THEN
          out=SIZE(arr, n)
       ELSE
          out=SIZE(arr)
       END IF
    ELSE
       out=0
    END IF
  END FUNCTION MySize_I4_2
  FUNCTION MySize_I4_3(arr,ndim) RESULT(out)
    INTEGER, ALLOCATABLE :: arr(:,:,:)
    INTEGER, OPTIONAL :: ndim
    INTEGER :: out
    INTEGER :: n
    n=0
    IF(PRESENT(ndim)) n=ndim
    IF(ALLOCATED(arr)) THEN
       IF(n /= 0) THEN
          out=SIZE(arr, n)
       ELSE
          out=SIZE(arr)
       END IF
    ELSE
       out=0
    END IF
  END FUNCTION MySize_I4_3
  FUNCTION MySize_R8_1(arr) RESULT(out)
    REAL(8), ALLOCATABLE :: arr(:)
    INTEGER :: out
    IF(ALLOCATED(arr)) THEN
       out=SIZE(arr)
    ELSE
       out=0
    END IF
  END FUNCTION MySize_R8_1
  FUNCTION MySize_R8_2(arr,ndim) RESULT(out)
    REAL(8), ALLOCATABLE :: arr(:,:)
    INTEGER, OPTIONAL :: ndim
    INTEGER :: out
    INTEGER :: n

    n=0
    IF(PRESENT(ndim)) n=ndim
    IF(ALLOCATED(arr)) THEN
       IF(n /= 0) THEN
          out=SIZE(arr, n)
       ELSE
          out=SIZE(arr)
       END IF
    ELSE
       out=0
    END IF
  END FUNCTION MySize_R8_2
  FUNCTION MySize_R8_3(arr,ndim) RESULT(out)
    REAL(8), ALLOCATABLE :: arr(:,:,:)
    INTEGER, OPTIONAL :: ndim
    INTEGER :: out
    INTEGER :: n
    n=0
    IF(PRESENT(ndim)) n=ndim
    IF(ALLOCATED(arr)) THEN
       IF(n /= 0) THEN
          out=SIZE(arr, n)
       ELSE
          out=SIZE(arr)
       END IF
    ELSE
       out=0
    END IF
  END FUNCTION MySize_R8_3
END MODULE Utilities
