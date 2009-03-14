MODULE NeiCONSTANTS
!!$***********************************************************************
!!$   Time-stamp: <2009-02-24 18:16:55 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Nov 24 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program ORAC ----*
  INTEGER, PARAMETER :: max_pars=200,max_char_tree = 80,&
       & max_char_long = 12000, max_data=12000,max_char=120,max_atm=10
  CHARACTER(len=1), DIMENSION(2), PARAMETER :: Comms=(/'!','#'/)
  CHARACTER(len=max_char), DIMENSION(9), PARAMETER :: Used=(/'ATOM  ','H&
       &ETATM','CONECT','SSBOND','CRYST1','SEQRES','HELIX ','SHEET ','CI&
       &SPEP'/)
  REAL(8), PARAMETER :: Huge=1.0D10,Tiny=1.0D-10
END MODULE NEICONSTANTS
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
  USE neiCONSTANTS
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
MODULE NeiERROR_List
  USE neiCONSTANTS, ONLY: max_char,max_data
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
  INTEGER, SAVE :: kprint=6
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
END MODULE NeiERROR_List
MODULE NeiErrors
  USE Node_Error
  USE NeiERROR_List
  IMPLICIT NONE 
  PRIVATE
  PUBLIC abort_now, abort_later, warning, add, print_errors,&
       & error_args,error_unr ,error_file,error_other&
       &,Setup_Errors,errmsg_f,errmsg_w, Print_warnings,List

  TYPE List
     CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE :: g
  END TYPE List
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
END MODULE NeiErrors
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
MODULE NeiNode
  USE NeiCONSTANTS
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
END MODULE NeiNode
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
MODULE Neighbors

!!$***********************************************************************
!!$   Time-stamp: <2007-01-19 18:57:55 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Jan 18 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program oracDD ----*

  USE NeiErrors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f
  USE NeiNode
  IMPLICIT none
  PRIVATE
  PUBLIC Neighbors_, Neighbors__Particles, ind_xyz, Neighbors__Ind&
       &, Neighbors__Chain, chain_xyz, Head_xyz

  TYPE :: Neighbors__Ind
     INTEGER :: i,j,k
  END TYPE Neighbors__Ind     
  TYPE :: Neighbors__Chain
     INTEGER :: i,j,k
     INTEGER :: p
  END TYPE Neighbors__Chain
  TYPE(Neighbors__Ind), ALLOCATABLE, SAVE :: Ind_xyz(:)
  TYPE(Neighbors__Chain), ALLOCATABLE, SAVE :: Chain_xyz(:)
  INTEGER, ALLOCATABLE, SAVE :: Head_xyz(:)
  INTEGER, SAVE :: ncx,ncy,ncz
CONTAINS
!!$
!!$--- Constructor for Ind_xyz
!!$
  FUNCTION Neighbors_(rcut,nx,ny,nz,co) RESULT(out)
    Real(8) :: co(3,3)
    LOGICAL :: out
    INTEGER :: nx,ny,nz
    REAL(8) :: rcut

    INTEGER :: vect0(3)
    INTEGER, POINTER :: vect(:)
    REAL(8) :: sqcut,dx,dy,dz,rmin
    INTEGER :: imax,jmax,kmax,i,j,k,istart,jstart,kstart,warnx,warny, warnz&
         &,nxmax, nymax,nzmax,nind,jend,kend

    IF(ALLOCATED(ind_xyz)) DEALLOCATE(ind_xyz)

    out=.TRUE.
    ncx=nx; ncy=ny; ncz=nz

    sqcut=rcut**2

    dx=2.d0/ncx
    dy=2.d0/ncy
    dz=2.d0/ncz
    imax=0
    jmax=0
    kmax=0

    vect0=(/0, 0, 0/) 
    IF(.NOT. Node_()) STOP

    CALL Node__Push(vect0)   

    istart=0
    DO i=istart,ncx-1
       jend=ncy-1
       IF(i == 0) jend=0
       DO j=1-ncy,jend
          kend=ncz-1
          IF(i == 0 .AND. j == 0) kend=0
          DO k=1-ncz,kend
             rmin=dist_ijk(i,j,k,dx,dy,dz)
             IF(rmin < sqcut) then
                vect0=(/i, j, k/) 
                IF(imax < abs(i)) imax=abs(i)
                IF(jmax < abs(j)) jmax=abs(j)
                IF(kmax < abs(k)) kmax=abs(k)
                IF(.NOT. (i == 0 .AND. j == 0 .AND. k == 0)) THEN
                   CALL Node__Push(vect0)
                END IF
             END IF
          END DO
       END DO
    END DO
    nind=Node__Size()
    ALLOCATE(ind_xyz(nind))
    nind=0
    DO WHILE(Node__Pop(vect))
       nind=nind+1
       ind_xyz(nind) % i=vect(1)
       ind_xyz(nind) % j=vect(2)
       ind_xyz(nind) % k=vect(3)
    END DO
    nxmax=(ncx+1)/2
    nymax=(ncy+1)/2
    nzmax=(ncz+1)/2
    warnx=0
    warny=0
    warnz=0
    IF(imax.ge.nxmax) warnx=1
    IF(jmax.ge.nymax) warny=1
    IF(kmax.ge.nzmax) warnz=1
    IF(warnx == 1 .OR. warny == 1 .OR. warnz == 1) THEN
       errmsg_f='Neighbor cells might be counted twice: Lower the&
            & cutoff or increase the No. of cell '
       IF(warnx == 1) THEN
          errmsg_f=TRIM(errmsg_f)//' along x '
       END IF
       IF(warny == 1) THEN
          errmsg_f=TRIM(errmsg_f)//' along y '
       END IF
       IF(warnz == 1) THEN
          errmsg_f=TRIM(errmsg_f)//' along z '
       END IF
       CALL Add_Errors(-1,errmsg_f)
       out=.TRUE. 
    END IF
  CONTAINS
    FUNCTION dist_ijk(ni,nj,nk,dx,dy,dz) RESULT(out)
      REAL(8) :: out
      REAL(8) ::  dx,dy,dz
      INTEGER ::  ni,nj,nk

      REAL(8) ::  d,dmin,dt
      REAL(8) ::  lx,ly,lz
      REAL(8) ::  mx,my,mz
      REAL(8) ::  dmx,dmy,dmz
      REAL(8) ::  msq,dmsq,lambda,s

      INTEGER, PARAMETER :: nv(8,3)=RESHAPE((/&
           & 0, 1, 0, 0, 1, 1, 0, 1 &
           &,0, 0, 1, 0, 1, 0, 1, 1 &
           &,0, 0, 0, 1, 0, 1, 1, 1 &
           &/),(/8, 3/))
      INTEGER ::  i,j,imin,jmin
      INTEGER, SAVE ::  ndtmax=0

!!$
!!$--- Minimum distance between corners of the cells (0, 0, 0) and 
!!$--- (ni, nj, nk)
!!$

      dmin=1.0D8
      do i=1,8
         do j=1,8
            lx=(ni+nv(j,1)-nv(i,1))*dx
            ly=(nj+nv(j,2)-nv(i,2))*dy
            lz=(nk+nv(j,3)-nv(i,3))*dz

            mx=co(1,1)*lx+co(1,2)*ly+co(1,3)*lz
            my=co(2,1)*lx+co(2,2)*ly+co(2,3)*lz
            mz=co(3,1)*lx+co(3,2)*ly+co(3,3)*lz

            d=mx*mx+my*my+mz*mz
            if(d.lt.dmin) then
               dmin=d
               imin=i
               jmin=j
            endif

         end do
      end do
!!$
!!$--- Check if the minimal distance is not on the edge of the cube 
!!$

      lx=(ni+nv(jmin,1)-nv(imin,1))*dx
      ly=(nj+nv(jmin,2)-nv(imin,2))*dy
      lz=(nk+nv(jmin,3)-nv(imin,3))*dz

      mx=co(1,1)*lx+co(1,2)*ly+co(1,3)*lz
      my=co(2,1)*lx+co(2,2)*ly+co(2,3)*lz
      mz=co(3,1)*lx+co(3,2)*ly+co(3,3)*lz

      msq=mx*mx+my*my+mz*mz
!!$
!!$--- Loop on the three edges
!!$

      do i=1,3
         dt=sign(1.,0.5-nv(imin,i))*dx
         
         dmx=co(1,i)*dt
         dmy=co(2,i)*dt
         dmz=co(3,i)*dt
         
         s=mx*dmx+my*dmy+mz*dmz
         dmsq=dmx*dmx+dmy*dmy+dmz*dmz
         lambda=-s/dmsq
         
         if((lambda.gt.0.) .and. (lambda.lt.1.)) then
            d=msq-s*s/dmsq
            if(d.lt.dmin) dmin=d
         endif
      enddo

      do i=1,3
         dt=sign(1.,0.5-nv(jmin,i))*dx
         
         dmx=co(1,i)*dt
         dmy=co(2,i)*dt
         dmz=co(3,i)*dt
         
         s=mx*dmx+my*dmy+mz*dmz
         dmsq=dmx*dmx+dmy*dmy+dmz*dmz
         lambda=-s/dmsq
         
         if((lambda.gt.0.) .and. (lambda.lt.1.)) then
            d=msq-s*s/dmsq
            if(d.lt.dmin) dmin=d
         endif
      enddo

      out=dmin
    END FUNCTION dist_ijk
  END FUNCTION Neighbors_
!!$
!!$--- Constructor for Chain_xyz and Head_xyz
!!$
  FUNCTION Neighbors__Particles(x,y,z,gr) RESULT(out)
    LOGICAL :: out
    INTEGER, OPTIONAL :: gr(:)
    REAL(8) :: x(:),y(:),z(:)
    REAL(8) :: x1,y1,z1,dx,dy,dz
    INTEGER :: n,nx,ny,nz,natp,numcell,l
    
    out=.TRUE.
    natp=SIZE(x)
    IF(.NOT. Neighbors__Valid()) THEN
       out=.FALSE.
       errmsg_f='Must construct cell indeces before atomic indeces can be obtained'
       CALL Add_Errors(-1,errmsg_f)
       RETURN
    END IF
    IF(.NOT. ALLOCATED(Head_xyz)) THEN
       ALLOCATE(Head_xyz(ncx*ncy*ncz))
       ALLOCATE(Chain_xyz(natp))
    ELSE
       DEALLOCATE(Head_xyz,Chain_xyz)
       ALLOCATE(Head_xyz(ncx*ncy*ncz))
       ALLOCATE(Chain_xyz(natp))
    END IF
    Head_xyz=0
    Chain_xyz (:) % p = 0

!!$=======================================================================
!!$     Compute chain list for the system
!!$=======================================================================

    dx=2.d0/ncx
    dy=2.d0/ncy
    dz=2.d0/ncz
    
    IF(PRESENT(gr))THEN
       DO n=1,natp
          Chain_xyz (n) % p=0
          IF(Gr(n) == 0) CYCLE
          x1=x(n)/dx
          y1=y(n)/dy
          z1=z(n)/dz
          nx=INT(x1)+(SIGN(1.D0,x1-INT(x1))-1.)/2
          ny=INT(y1)+(SIGN(1.D0,y1-INT(y1))-1.)/2
          nz=INT(z1)+(sign(1.d0,z1-int(z1))-1.)/2
          nx=MOD(MOD(nx,ncx)+ncx,ncx)
          ny=MOD(MOD(ny,ncy)+ncy,ncy)
          nz=MOD(MOD(nz,ncz)+ncz,ncz)
          Chain_xyz (n) % i=nx
          Chain_xyz (n) % j=ny
          Chain_xyz (n) % k=nz
          numcell=nz+ncz*(ny+ncy*nx)+1
          Chain_xyz (n) % p=Head_xyz(numcell)
          Head_xyz(numcell)=n
       END DO
    ELSE
       DO n=1,natp
          x1=x(n)/dx
          y1=y(n)/dy
          z1=z(n)/dz
          nx=INT(x1)+(SIGN(1.D0,x1-INT(x1))-1.)/2
          ny=INT(y1)+(SIGN(1.D0,y1-INT(y1))-1.)/2
          nz=INT(z1)+(sign(1.d0,z1-int(z1))-1.)/2
          nx=MOD(MOD(nx,ncx)+ncx,ncx)
          ny=MOD(MOD(ny,ncy)+ncy,ncy)
          nz=MOD(MOD(nz,ncz)+ncz,ncz)
          Chain_xyz (n) % i=nx
          Chain_xyz (n) % j=ny
          Chain_xyz (n) % k=nz
          numcell=nz+ncz*(ny+ncy*nx)+1
          Chain_xyz (n) % p=Head_xyz(numcell)
          Head_xyz(numcell)=n       
       END DO
    END IF
  END FUNCTION Neighbors__Particles
  SUBROUTINE Neighbors__Delete
    DEALLOCATE(Ind_xyz,Chain_xyz,Head_xyz)
  END SUBROUTINE Neighbors__Delete
  FUNCTION Neighbors__Valid() RESULT(out)
    LOGICAL :: out
    out=ALLOCATED(ind_xyz)
  END FUNCTION Neighbors__Valid
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE Neighbors
