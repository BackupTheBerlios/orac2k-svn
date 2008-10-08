!!$/---------------------------------------------------------------------\
!!$   Copyright  © 2006-2007 Massimo Marchi <Massimo.Marchi at cea.fr>   |
!!$                                                                      |
!!$    This software is a computer program named oracDD whose            |
!!$    purpose is to simulate and model complex molecular systems.       |
!!$    The code is written in fortran 95 compliant with Technical        |
!!$    Report TR 15581, and uses MPI-1 routines for parallel             |
!!$    coding.                                                           |
!!$                                                                      |
!!$    This software is governed by the CeCILL license under             |
!!$    French law and abiding by the rules of distribution of            |
!!$    free software.  You can  use, modify and/ or redistribute         |
!!$    the software under the terms of the CeCILL icense as              |
!!$    circulated by CEA, CNRS and INRIA at the following URL            |
!!$    "http://www.cecill.info".                                         |
!!$                                                                      |
!!$    As a counterpart to the access to the source code and rights      |
!!$    to copy, modify and redistribute granted by the license,          |
!!$    users are provided only with a limited warranty and the           |
!!$    software's author, the holder of the economic rights, and         |
!!$    the successive licensors have only limited liability.             |
!!$                                                                      |
!!$    The fact that you are presently reading this means that you       |
!!$    have had knowledge of the CeCILL license and that you accept      |
!!$    its terms.                                                        |
!!$                                                                      |
!!$    You should have received a copy of the CeCILL license along       |
!!$    with this program; if not, you can collect copies on the URL's    |
!!$    "http://www.cecill.info/licences/Licence_CeCILL_V2-en.html"       |
!!$    "http://www.cecill.info/licences/Licence_CeCILL_V2-fr.html"       |
!!$                                                                      |
!!$----------------------------------------------------------------------/
MODULE Node
  USE CONSTANTS
  PRIVATE
  PUBLIC :: Node_,Node__Delete,Node__Push,Node__Pop&
       &,Node__Size
  TYPE :: NODES
     CHARACTER(len=1), DIMENSION(:), POINTER  :: vector=>NULL()
     TYPE(NODES), POINTER :: next=>NULL()
  END TYPE NODES

  TYPE(NODES), POINTER :: root=>NULL()
  TYPE(NODES), POINTER :: current=>NULL()
  INTEGER, SAVE :: counter=0
  INTERFACE Node__Push
     MODULE PROCEDURE Node_I__Push
     MODULE PROCEDURE Node_I_1__Push
     MODULE PROCEDURE Node_R8__Push
     MODULE PROCEDURE Node_R8_1__Push
     MODULE PROCEDURE Node_CL__Push
     MODULE PROCEDURE Node_CL_1__Push
  END INTERFACE
  INTERFACE Node__Pop
     MODULE PROCEDURE Node_I_1__Pop
     MODULE PROCEDURE Node_I__Pop
     MODULE PROCEDURE Node_R8__Pop
     MODULE PROCEDURE Node_R8_1__Pop
     MODULE PROCEDURE Node_CL__Pop
     MODULE PROCEDURE Node_CL_1__Pop
  END INTERFACE
CONTAINS
  FUNCTION Node_() RESULT(out)
    IMPLICIT none
    LOGICAL :: out
    IF(counter /= 0) THEN
       counter=0
       CALL Node__Delete()
    END IF
    ALLOCATE(root)
    NULLIFY(root%next)
    current=>root
    out=.TRUE.
  END FUNCTION Node_

  FUNCTION Node__Size() RESULT(out)
    INTEGER :: out
    out=counter
  END FUNCTION Node__Size

  SUBROUTINE Node__delete()
    CALL Remove(root)
    NULLIFY(root)
    NULLIFY(current)
  END SUBROUTINE Node__delete
  RECURSIVE SUBROUTINE Remove(old_node)
    TYPE(NODES), POINTER :: old_node
    IF(.NOT. ASSOCIATED(old_node) ) RETURN
    CALL Remove(old_node % next)
    IF(ASSOCIATED(old_node % vector)) DEALLOCATE(old_node % vector)
    DEALLOCATE(old_node)
  END SUBROUTINE Remove

!!$
!!$--- Overloading 
!!$

  SUBROUTINE Node_C__Push(vect)
    IMPLICIT none
    CHARACTER(len=1), DIMENSION(:), INTENT(IN) :: vect

    TYPE(NODES), POINTER :: new_node
    counter=counter+1
    ALLOCATE(new_node)
    ALLOCATE(new_node % vector(SIZE(vect)))
    new_node % vector = vect
    NULLIFY(new_node % next)

    current % next => new_node
    current => current % next

  END SUBROUTINE Node_C__Push

  FUNCTION Node_C__Pop(vect) RESULT(out)
    IMPLICIT none
    LOGICAL :: out
    CHARACTER(len=1), DIMENSION(:), POINTER :: vect
    INTEGER, SAVE :: count0=0

    IF(count0 == 0) THEN
       current=>root
    END IF
    
    count0=count0+1
    out=ASSOCIATED(current % next )
    IF(.NOT. out) THEN
       count0=0
       RETURN
    END IF
 
    current=>current % next
    IF(ASSOCIATED(vect)) DEALLOCATE(Vect)
    ALLOCATE(Vect(SIZE(current % vector)))
    vect=current % vector

  END FUNCTION Node_C__Pop

!!$
!!$--- Hack to get around bug on the TRANSFER() function of the pgi compiler!
!!$--- version 6.2-5
!!$

  FUNCTION MyTransfer_Str2Char(vect) RESULT(out)
    CHARACTER(len=Max_Char), DIMENSION(:) :: vect
    CHARACTER(len=1), DIMENSION(:), POINTER :: out
    CHARACTER(len=1), DIMENSION(:), ALLOCATABLE, TARGET, SAVE :: cvect0
    INTEGER :: n,m,p

    n=Max_char*SIZE(vect)
    IF(ALLOCATED(cvect0)) DEALLOCATE(cvect0)
    ALLOCATE(cvect0(n))
    p=0
    DO n=1,SIZE(vect)
       DO m=1,LEN(vect(n))
          p=p+1
          cvect0(p)=vect(n)(m:m)
       END DO
    END DO
    out=>cvect0
  END FUNCTION MyTransfer_Str2Char

  FUNCTION MyTransfer_Char2Str(cvect) RESULT(out)
    CHARACTER(len=1), DIMENSION(:) :: cvect

    CHARACTER(len=Max_Char), DIMENSION(:), POINTER :: out
    CHARACTER(len=Max_Char), DIMENSION(:), ALLOCATABLE, TARGET, SAVE :: vect0
    INTEGER :: n,m,p

    n=SIZE(cvect)/Max_Char
    IF(ALLOCATED(vect0)) DEALLOCATE(vect0)
    ALLOCATE(vect0(n))

    p=0
    DO n=1,SIZE(vect0)
       DO m=1,LEN(vect0(n))
          p=p+1
          vect0(n)(m:m)=cvect(p)
       END DO
    END DO
    out=>vect0
  END FUNCTION MyTransfer_Char2Str
!!$
!!$--- Hack to get around bug in the pgi compiler end!
!!$

!!$
!!$--- Overloaded routines
!!$

  SUBROUTINE Node_CL__Push(vect)
    CHARACTER(len=max_Char), DIMENSION(:), INTENT(IN) :: vect
    CHARACTER(len=1), DIMENSION(:), POINTER :: cvect

    cvect=>MYTransfer_Str2Char(vect)

    CALL Node_C__Push(cvect)

  END SUBROUTINE Node_CL__Push

  SUBROUTINE Node_CL_1__Push(vect)
    CHARACTER(len=max_Char), INTENT(IN) :: vect
    CHARACTER(len=1), DIMENSION(:), POINTER :: cvect
    CHARACTER(len=max_Char) :: vect0(1)

    vect0(1)=vect
    cvect=>MYTransfer_Str2Char(vect0)

    CALL Node_C__Push(cvect)

  END SUBROUTINE Node_CL_1__Push

  FUNCTION Node_CL__Pop(vect) RESULT(out)
    IMPLICIT none
    LOGICAL :: out
    CHARACTER(len=max_char), DIMENSION(:), POINTER :: vect
    CHARACTER(len=max_char), DIMENSION(:), POINTER :: aux
    CHARACTER(len=1), DIMENSION(:), POINTER :: cvect
    
    out=Node_C__Pop(cvect)

    IF(.NOT. out) RETURN

    IF(ASSOCIATED(vect)) DEALLOCATE(vect)
    aux=>MyTransfer_Char2Str(cvect)
    ALLOCATE(vect(SIZE(aux)))
    vect=aux
    DEALLOCATE(cvect)
  END FUNCTION Node_CL__Pop
  FUNCTION Node_CL_1__Pop(vect) RESULT(out)
    IMPLICIT none
    LOGICAL :: out
    CHARACTER(len=max_char), INTENT(OUT) :: vect
    CHARACTER(len=max_char), DIMENSION(:), POINTER :: aux
    CHARACTER(len=1), DIMENSION(:), POINTER :: cvect
    
    out=Node_C__Pop(cvect)
    IF(.NOT. out) RETURN

    aux=>MyTransfer_Char2Str(cvect)
    vect=aux(1)
    DEALLOCATE(cvect)
  END FUNCTION Node_CL_1__Pop

  SUBROUTINE Node_R8__Push(vect)
    REAL(8), DIMENSION(:), INTENT(IN) :: vect
    CHARACTER(len=1), DIMENSION(:), ALLOCATABLE :: cvect1
    CHARACTER(len=1), DIMENSION(1) :: cvect
    INTEGER :: length

    length=SIZE(TRANSFER(vect,cvect))

    ALLOCATE(cvect1(length))
    cvect1=TRANSFER(vect,cvect1)
    CALL Node_C__Push(cvect1)
  END SUBROUTINE Node_R8__Push

  SUBROUTINE Node_R8_1__Push(vect)
    REAL(8), INTENT(IN) :: vect
    CHARACTER(len=1), DIMENSION(:), ALLOCATABLE :: cvect1
    CHARACTER(len=1), DIMENSION(1) :: cvect
    INTEGER :: length
    length=SIZE(TRANSFER(vect,cvect))
    ALLOCATE(cvect1(length))
    cvect1=TRANSFER(vect,cvect1)
    CALL Node_C__Push(cvect1)
  END SUBROUTINE Node_R8_1__Push
  FUNCTION Node_R8__Pop(vect) RESULT(out)
    IMPLICIT none
    LOGICAL :: out
    REAL(8), DIMENSION(:), POINTER :: vect
    REAL(8), DIMENSION(:), POINTER  :: aux
    CHARACTER(len=1), DIMENSION(:), POINTER :: cvect
    
    out=Node_C__Pop(cvect)
    IF(.NOT. out) RETURN
    IF(ASSOCIATED(vect)) DEALLOCATE(vect)
    ALLOCATE(vect(SIZE(TRANSFER(cvect,aux))))
    vect=TRANSFER(cvect,vect)
    DEALLOCATE(cvect)
  END FUNCTION Node_R8__Pop
  FUNCTION Node_R8_1__Pop(vect) RESULT(out)
    IMPLICIT none
    LOGICAL :: out
    REAL(8), INTENT(OUT) :: vect
    CHARACTER(len=1), DIMENSION(:), POINTER :: cvect
    
    out=Node_C__Pop(cvect)
    IF(.NOT. out) RETURN
    vect=TRANSFER(cvect,vect)
    DEALLOCATE(cvect)
  END FUNCTION Node_R8_1__Pop
  SUBROUTINE Node_I__Push(vect)
    INTEGER, DIMENSION(:), INTENT(IN) :: vect
    CHARACTER(len=1), DIMENSION(:), ALLOCATABLE :: cvect1
    CHARACTER(len=1), DIMENSION(1) :: cvect
    INTEGER :: length
    length=SIZE(TRANSFER(vect,cvect))
    ALLOCATE(cvect1(length))
    cvect1=TRANSFER(vect,cvect1)
    CALL Node_C__Push(cvect1)
  END SUBROUTINE Node_I__Push
  SUBROUTINE Node_I_1__Push(vect)
    INTEGER, INTENT(IN) :: vect
    CHARACTER(len=1), DIMENSION(:), ALLOCATABLE :: cvect1
    CHARACTER(len=1), DIMENSION(1) :: cvect
    INTEGER :: length
    length=SIZE(TRANSFER(vect,cvect))
    ALLOCATE(cvect1(length))
    cvect1=TRANSFER(vect,cvect1)
    CALL Node_C__Push(cvect1)
  END SUBROUTINE Node_I_1__Push
  FUNCTION Node_I__Pop(vect) RESULT(out)
    IMPLICIT none
    LOGICAL :: out
    INTEGER, DIMENSION(:), POINTER :: vect
    INTEGER, DIMENSION(:), POINTER  :: aux
    CHARACTER(len=1), DIMENSION(:), POINTER :: cvect
    
    out=Node_C__Pop(cvect)
    IF(.NOT. out) RETURN
    IF(ASSOCIATED(vect)) DEALLOCATE(vect)

    ALLOCATE(vect(SIZE(TRANSFER(cvect,aux))))
    vect=TRANSFER(cvect,vect)
    DEALLOCATE(cvect)
  END FUNCTION Node_I__Pop
  FUNCTION Node_I_1__Pop(vect) RESULT(out)
    IMPLICIT none
    LOGICAL :: out
    INTEGER, INTENT(OUT) :: vect
    CHARACTER(len=1), DIMENSION(:), POINTER :: cvect
    
    out=Node_C__Pop(cvect)
    IF(.NOT. out) RETURN
    vect=TRANSFER(cvect,vect)
    DEALLOCATE(cvect)
  END FUNCTION Node_I_1__Pop
END MODULE Node
