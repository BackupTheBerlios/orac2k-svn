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
  USE Node
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
    nn_old=0
    IF(.NOT. Node_()) STOP
    
    DO WHILE(nn_old+nlast < n)
       nn_old=nn
       nlast=MIN(nn_old+nl,n)
       nn=SCAN(msg(1:nlast),' ',BACK=.TRUE.)
       i_nn=nn-1-nn_old
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
    
    WRITE(*,'(2x,a)') msg2%top
    WRITE(*,'(2x,a)') msg2%body(1:m)
    WRITE(*,'(2x,a)') msg2%bottom
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
       nn_old=0
       IF(.NOT. Node_()) STOP

       DO WHILE(nn_old+nlast < n)
          nn_old=nn
          nlast=MIN(nn_old+nl,n)
          nn=SCAN(msg(1:nlast),' ',BACK=.TRUE.)
          i_nn=nn-1-nn_old
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

       IF(counter == 0) WRITE(*,'(2x,a)') msg2%top
       WRITE(*,'(2x,a)') msg2%intrabodies
       WRITE(*,'(2x,a)') msg2%body(1:m)
    ELSE
       WRITE(*,'(2x,a)') msg2%bottom
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
       nn_old=0
       IF(.NOT. Node_()) STOP

       DO WHILE(nn_old+nlast < n)
          nn_old=nn
          nlast=MIN(nn_old+nl,n)
          nn=SCAN(msg(1:nlast),' ',BACK=.TRUE.)
          i_nn=nn-1-nn_old
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

       IF(counter == 0) WRITE(*,'(2x,a)') msg2%top
       WRITE(*,'(2x,a)') msg2%intrabodies
       WRITE(*,'(2x,a)') msg2%body(1:m)
    ELSE
       WRITE(*,'(2x,a)') msg2%bottom
    END IF
    IF(ASSOCIATED(vector)) DEALLOCATE(vector)
    IF(ALLOCATED(msg2 % body)) DEALLOCATE(msg2 % body)
    counter=counter+1
  END SUBROUTINE Warning
  SUBROUTINE Print_Warnings()
    LOGICAL :: warnings_found
    CHARACTER(len=max_err_long) :: error
    CHARACTER(len=10) :: dummy
    INTEGER :: No,count,No_Warn

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
          No_Warn=No_Warn+1
          IF(No_Warn == 1) CALL Warning(' ',0)
          CALL Warning(Error)
          warnings_found=.TRUE.
       END IF
       current=>current % next
    END DO
    IF(warnings_found) THEN
       CALL Warning(); CALL CleanUp()
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
