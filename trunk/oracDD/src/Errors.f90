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
MODULE Node_Error
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
END MODULE Node_Error
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
MODULE ERROR_List
  USE CONSTANTS, ONLY: max_char,max_data
  USE Print_Defs
  IMPLICIT NONE 
  INTEGER, SAVE :: counter=0, count_f=0, count_w=0
  INTEGER, PARAMETER :: max_err=4*max_char, max_err_long=12000
  INTEGER, PARAMETER :: Max_errors=50
  TYPE Orac_Errors
     INTEGER :: tag
     INTEGER :: count
     CHARACTER(len=max_err_long) :: err_text
     TYPE(Orac_Errors), POINTER :: next
  END TYPE Orac_Errors
  TYPE(Orac_Errors),  POINTER, SAVE :: root, current
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
  SUBROUTINE Add(No,error)
    IMPLICIT none
    CHARACTER(len=*), INTENT(IN) :: error
    INTEGER, INTENT(IN) :: No
    IF(counter == 0) THEN
       CALL Start()
    END IF
    IF(counter > Max_errors) RETURN
    counter=counter+1
    IF(No < 0) THEN
       count_f=count_f+1
       current % count = count_f
    ELSE
       count_w=count_w+1
       current % count = count_w
    END IF
    current % tag = No
    current % err_text = error
    ALLOCATE(current % next)
    current=>current % next
    NULLIFY(current % next)
  END SUBROUTINE Add
  SUBROUTINE Cleanup()
    IMPLICIT none
    TYPE(Orac_Errors), POINTER :: dummy
    
    current=>root
    DO WHILE(ASSOCIATED(current % next))
       dummy=>current % next
       DEALLOCATE(current)
       current=>dummy
    END DO
    NULLIFY(dummy, root)
    count_w=0
  END SUBROUTINE Cleanup
END MODULE ERROR_List
MODULE Errors
  USE Node_Error
  USE ERROR_List
  USE TYPES, ONLY: list
  IMPLICIT NONE 
  PRIVATE
  PUBLIC abort_now, abort_later, warning, add, print_errors,&
       & error_args,error_unr ,error_file,error_other&
       &,Setup_Errors,errmsg_f,errmsg_w, Print_warnings

  CHARACTER(LEN=45), SAVE :: err_arg_no='.. Arguments to command must be at least '
  TYPE fatal_err
     CHARACTER(74), DIMENSION (3) :: top= &
          &(/ '**************************************************************************'&
          &,'*                                                                        *'&
          &,'*-- The following FATAL errors were found:         ----------------------*'/)
     CHARACTER(74), DIMENSION (:), ALLOCATABLE :: body
     CHARACTER(74), DIMENSION(1) :: intrabodies = &
          (/ '*                                                                        *'/)
     CHARACTER(74), DIMENSION(2) :: bottom= &
          (/ '*                                                                        *' &
          ,'**************************************************************************' /)
  END TYPE fatal_err

  TYPE warning_err
     CHARACTER(74), DIMENSION (3) :: top= &
          (/ '**************************************************************************' &
          &,'*                                                                        *'&
          &,'*-- The following WARNING errors were found:         --------------------*'/)
     CHARACTER(74), DIMENSION(1) :: intrabodies = &
          (/ '*                                                                        *'/)
     CHARACTER(74), DIMENSION (:), ALLOCATABLE :: body
     CHARACTER(74), DIMENSION(2) :: bottom= &
          (/ '*                                                                        *' &
          ,'**************************************************************************' /)
  END TYPE warning_err

  TYPE(List), SAVE :: error_args,error_unr,error_file,error_other
  CHARACTER(len=max_data) :: errmsg_f,errmsg_w
CONTAINS

  SUBROUTINE abort_now(msg)
    !*******************************************************************
    !                                                                  *
    !  Abort NOW.                                                      *
    !  msg1: Main Error message                                        *
    !                                                                  *
    !*******************************************************************
    IMPLICIT NONE
    CHARACTER(*)  :: msg

    TYPE(fatal_err) :: msg2
    INTEGER :: n,m,i,mn,i_nn,nn,nn_old,nlast
    INTEGER, PARAMETER :: nl=66
    CHARACTER(len=max_char) :: vector0
    CHARACTER(len=max_char), POINTER :: vector(:)=>NULL()

    n=LEN_TRIM(msg)+1 ! Add a space to stop DO WHILE at the end of the string

!!$-- Make the error box

    nn=0
    nlast=0
    IF(.NOT. Node_()) STOP
    DO WHILE(nlast /= n)
       nn_old=nn          
       nlast=MIN(nn+nl,n)
       nn=SCAN(msg(1:nlast),' ',BACK=.TRUE.)
       i_nn=nn-nn_old-1
       vector0=' '
       vector0(1:1)='*'
       vector0(74:74)='*'
       vector0(5:i_nn+5-1)=msg(nn_old+1:nn-1)
       CALL Node__Push(vector0)
    END DO
        
    m=Node__Size()
    
    ALLOCATE(msg2 % body(m))
    i=0
    DO WHILE(Node__Pop(vector))
       i=i+1
       msg2%body(i)=vector(1)
    END DO
    
    WRITE(kprint,'(2x,a)') msg2%top
    WRITE(kprint,'(2x,a)') msg2%body(1:m)
    WRITE(kprint,'(2x,a)') msg2%bottom
    IF(ASSOCIATED(vector)) DEALLOCATE(vector)
    IF(ALLOCATED(msg2 % body)) DEALLOCATE(msg2 % body)
    
    STOP

  END SUBROUTINE abort_now
  SUBROUTINE abort_later(msg)
!!$!*******************************************************************
!!$!                                                                  *
!!$!  Abort.                                                          *
!!$!  msg1: Main Error message                                        *
!!$!                                                                  *
!!$!*******************************************************************
    IMPLICIT NONE
    CHARACTER(*), OPTIONAL  :: msg

    TYPE(fatal_err) :: msg2
    INTEGER :: n,m,i,mn,i_nn,nn,nn_old,nlast
    INTEGER, PARAMETER :: nl=66
    INTEGER, SAVE :: counter=0
    CHARACTER(len=max_char) :: vector0
    CHARACTER(len=max_char), POINTER :: vector(:)=>NULL()
    
    IF(PRESENT(msg)) THEN
       n=LEN_TRIM(msg)+1 ! Add a space to stop DO WHILE at the end of the string

!!$-- Make the error box

       nn=0
       nlast=0
       IF(.NOT. Node_()) STOP
       DO WHILE(nlast /= n)
          nn_old=nn          
          nlast=MIN(nn+nl,n)
          nn=SCAN(msg(1:nlast),' ',BACK=.TRUE.)
          i_nn=nn-nn_old-1
          vector0=' '
          vector0(1:1)='*'
          vector0(74:74)='*'
          vector0(5:i_nn+5-1)=msg(nn_old+1:nn-1)
          CALL Node__Push(vector0)
       END DO
       
       m=Node__Size()
       ALLOCATE(msg2 % body(m))
       i=0
       DO WHILE(Node__Pop(vector))
          i=i+1
          msg2%body(i)=vector(1)
       END DO

       IF(counter == 0) WRITE(kprint,'(2x,a)') msg2%top
       WRITE(kprint,'(2x,a)') msg2%intrabodies
       WRITE(kprint,'(2x,a)') msg2%body(1:m)
    ELSE
       WRITE(kprint,'(2x,a)') msg2%bottom
    END IF
    counter=counter+1
    IF(ASSOCIATED(vector)) DEALLOCATE(vector)
    IF(ALLOCATED(msg2 % body)) DEALLOCATE(msg2 % body)
  END SUBROUTINE abort_later
  SUBROUTINE Warning(msg, reset)
!!$!*******************************************************************
!!$!                                                                  *
!!$!  Warn.                                                           *
!!$!  msg1: Main Error message                                        *
!!$!                                                                  *
!!$!*******************************************************************
    IMPLICIT NONE
    CHARACTER(*), OPTIONAL  :: msg
    INTEGER, OPTIONAL  :: reset

    TYPE(warning_err) :: msg2
    INTEGER :: n,m,i,mn,i_nn,nn,nn_old,nlast
    INTEGER, PARAMETER :: nl=66
    INTEGER, SAVE :: counter=0
    CHARACTER(len=max_char) :: vector0
    CHARACTER(len=max_char), POINTER :: vector(:)=>NULL()


    IF(PRESENT(reset)) THEN
       counter=0    
       RETURN
    END IF
    IF(PRESENT(msg)) THEN
       n=LEN_TRIM(msg)+1 ! Add a space to stop DO WHILE at the end of the string
       
!!$-- Make the error box


       nn=0
       nlast=0
       IF(.NOT. Node_()) STOP
       DO WHILE(nlast /= n)
          nn_old=nn          
          nlast=MIN(nn+nl,n)
          nn=SCAN(msg(1:nlast),' ',BACK=.TRUE.)
          i_nn=nn-nn_old-1
          vector0=' '
          vector0(1:1)='*'
          vector0(74:74)='*'
          vector0(5:i_nn+5-1)=msg(nn_old+1:nn-1)
          CALL Node__Push(vector0)
       END DO
       
       m=Node__Size()
       ALLOCATE(msg2 % body(m))
       i=0
       DO WHILE(Node__Pop(vector))
          i=i+1
          msg2%body(i)=vector(1)
       END DO

       IF(counter == 0) WRITE(kprint,'(2x,a)') msg2%top
       WRITE(kprint,'(2x,a)') msg2%intrabodies
       WRITE(kprint,'(2x,a)') msg2%body(1:m)
    ELSE
       WRITE(kprint,'(2x,a)') msg2%bottom
    END IF
    IF(ASSOCIATED(vector)) DEALLOCATE(vector)
    IF(ALLOCATED(msg2 % body)) DEALLOCATE(msg2 % body)
    counter=counter+1
  END SUBROUTINE Warning
  SUBROUTINE Print_Warnings()
    LOGICAL :: warnings_found
    CHARACTER(len=max_err_long) :: error
    CHARACTER(len=10) :: dummy
    INTEGER :: No,No_Warn,count

    error=' '
    No_Warn=0
    warnings_found=.FALSE.
    current=>root
    DO WHILE(ASSOCIATED(current % next))
       No=current % tag
       count=current % count
       WRITE(dummy,'(i3,'') --'')') count
       Error=TRIM(dummy)//' '//current % err_text
       IF(No > 0) THEN
          current % tag = 0
          No_Warn=No_Warn+1
          IF(No_Warn == 1) CALL Warning(' ',0)
          CALL Warning(Error)
          warnings_found=.TRUE.
       END IF
       current=>current % next
    END DO
    IF(warnings_found) THEN
       CALL Warning(); count_w=0
    END IF
  END SUBROUTINE Print_Warnings
  SUBROUTINE Print_Errors()
    LOGICAL :: stop_run
    CHARACTER(len=max_err_long) :: error
    CHARACTER(len=10) :: dummy
    INTEGER :: No,count

    
    error=' '
    IF(.NOT. ASSOCIATED(root)) RETURN
    stop_run=.FALSE.
    current=>root    
    DO WHILE(ASSOCIATED(current % next))
       No=current % tag
       count=current % count
       WRITE(dummy,'(i3,'') --'')') count
       Error=TRIM(dummy)//' '//current % err_text
       IF(No < 0) THEN
          stop_run=.TRUE.
          CALL Abort_Later(Error)
       END IF
       current=>current % next
    END DO
    IF(stop_run) THEN
       CALL Abort_Later
       STOP
    END IF
  END SUBROUTINE Print_Errors
  SUBROUTINE Setup_Errors
    ALLOCATE(error_args % g (4),error_unr % g (4), error_file % g (1)&
         &, error_other % g (1)) 

    error_args % g (1)= 'OPEN keyword not found '
    error_args % g (2)= 'Number of arguments must be at least '
    error_args % g (3)= 'Number of arguments must not exceed '
    error_args % g (4)= 'Number of arguments must be only '
    error_unr % g (1)='Unrecognized command  ---> '
    error_unr % g (2)='Unrecognized subcommand -> '
    error_unr % g (3)='Unrecognized keyword ----> '
    error_unr % g (4)='UNSUPPORTED  COMMAND ----> '
    error_file % g (1)='...or missing &END'
    error_other % g (1)=' file not found '

  END SUBROUTINE Setup_Errors
END MODULE Errors
