MODULE ERROR_List
  USE CONSTANTS, ONLY: max_char
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
  TYPE(Orac_Errors),  POINTER :: root, current
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
    INTEGER :: count
    
    current=>root
    DO WHILE(ASSOCIATED(current % next))
       dummy=>current % next
       DEALLOCATE(current)
       current=>dummy
    END DO
    NULLIFY(dummy, root)
  END SUBROUTINE Cleanup
END MODULE ERROR_List
MODULE ERROR_Mod
  USE ERROR_List
  USE TYPES, ONLY: list
  IMPLICIT NONE 
  PRIVATE
  PUBLIC abort_now, abort_later, warning, add, print_errors,&
       & error_args,error_unr ,error_file,error_other,Setup_Errors

  CHARACTER(LEN=45), SAVE :: err_arg_no='.. Arguments to command must be at least '
  TYPE fatal_err
     CHARACTER(74), DIMENSION (3) :: top= &
          &(/ '**************************************************************************'&
          &,'*                                                                        *'&
          &,'*-- The following FATAL errors were found:         ----------------------*'/)
     CHARACTER(74), DIMENSION (:), POINTER :: body
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
          ,'*-- The following WARNING errors were found:         --------------------*'/)
     CHARACTER(74), DIMENSION(1) :: intrabodies = &
          (/ '*                                                                        *'/)
     CHARACTER(74), DIMENSION (:), POINTER :: body
     CHARACTER(74), DIMENSION(2) :: bottom= &
          (/ '*                                                                        *' &
          ,'**************************************************************************' /)
  END TYPE warning_err

  TYPE(List), SAVE :: error_args,error_unr,error_file,error_other
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
    INTEGER :: n,m,i,mn
    INTEGER, PARAMETER :: nl=66

    n=LEN_TRIM(msg)

    m=INT(CEILING(FLOAT(n)/FLOAT(nl)))
    ! Make error box

    ALLOCATE(msg2%body(m))

    msg2%body(:)=' '
    msg2%body(:)(1:1)='*'
    msg2%body(:)(74:74)='*'

    DO i=1,m-1
       msg2%body(i)(5:nl+5-1)=msg((i-1)*nl+1:i*nl)
    END DO
    mn=nl
    IF(MOD(n,nl) /= 0) mn=MOD(n,nl)
    
    msg2%body(m)(5:nl+5-1)=' '
    msg2%body(m)(5:5+mn-1)=msg((m-1)*nl+1:n)

    WRITE(*,'(2x,a)') msg2%top
    WRITE(*,'(2x,a)') msg2%body(1:m)
    WRITE(*,'(2x,a)') msg2%bottom
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
    INTEGER :: n,m,i,mn
    INTEGER, PARAMETER :: nl=66
    INTEGER, SAVE :: counter=0
    
    IF(PRESENT(msg)) THEN
       n=LEN_TRIM(msg)

       m=INT(CEILING(FLOAT(n)/FLOAT(nl)))
       ! Make error box
       
       ALLOCATE(msg2%body(m))
       
       msg2%body(:)=' '
       msg2%body(:)(1:1)='*'
       msg2%body(:)(74:74)='*'

       DO i=1,m-1
          msg2%body(i)(5:nl+5-1)=msg((i-1)*nl+1:i*nl)
       END DO
       mn=nl
       IF(MOD(n,nl) /= 0) mn=MOD(n,nl)
       
       msg2%body(m)(5:nl+5-1)=' '
       msg2%body(m)(5:5+mn-1)=msg((m-1)*nl+1:n)

       IF(counter == 0) WRITE(*,'(2x,a)') msg2%top
       WRITE(*,'(2x,a)') msg2%intrabodies
       WRITE(*,'(2x,a)') msg2%body(1:m)
    ELSE
       WRITE(*,'(2x,a)') msg2%bottom
    END IF
    counter=counter+1
  END SUBROUTINE abort_later
  SUBROUTINE Warning(msg)
!!$!*******************************************************************
!!$!                                                                  *
!!$!  Warn.                                                           *
!!$!  msg1: Main Error message                                        *
!!$!                                                                  *
!!$!*******************************************************************
    IMPLICIT NONE
    CHARACTER(*), OPTIONAL  :: msg

    TYPE(warning_err) :: msg2
    INTEGER :: n,m,i,mn
    INTEGER, PARAMETER :: nl=66
    INTEGER, SAVE :: counter=0

    IF(PRESENT(msg)) THEN
       n=LEN_TRIM(msg)
       
       m=INT(CEILING(FLOAT(n)/FLOAT(nl)))

!!$-- Make the error box
       
       ALLOCATE(msg2%body(m))
       
       msg2%body(:)=' '
       msg2%body(:)(1:1)='*'
       msg2%body(:)(74:74)='*'
       
       DO i=1,m-1
          msg2%body(i)(5:nl+5-1)=msg((i-1)*nl+1:i*nl)
       END DO
       mn=nl
       IF(MOD(n,nl) /= 0) mn=MOD(n,nl)
       
       msg2%body(m)(5:nl+5-1)=' '
       msg2%body(m)(5:5+mn-1)=msg((m-1)*nl+1:n)
       
       IF(counter == 0) WRITE(*,'(2x,a)') msg2%top
       WRITE(*,'(2x,a)') msg2%intrabodies
       WRITE(*,'(2x,a)') msg2%body(1:m)
    ELSE
       WRITE(*,'(2x,a)') msg2%bottom
    END IF
    counter=counter+1
  END SUBROUTINE Warning
  SUBROUTINE Print_Errors()
    IMPLICIT none
    LOGICAL :: stop_run,warnings_found
    CHARACTER(len=max_err_long) :: error
    CHARACTER(len=10) :: dummy
    INTEGER :: No,count
    
    error=' '
    IF(.NOT. ASSOCIATED(root)) RETURN
    stop_run=.FALSE.
    warnings_found=.FALSE.

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
    current=>root    
    DO WHILE(ASSOCIATED(current % next))
       No=current % tag
       count=current % count
       WRITE(dummy,'(i3,'') --'')') count
       Error=TRIM(dummy)//' '//current % err_text
       IF(No > 0) THEN
          CALL Warning(Error)
          warnings_found=.TRUE.
       END IF
       current=>current % next
    END DO
    IF(stop_run) THEN
       CALL Abort_Later
       STOP
    END IF
    IF(warnings_found) CALL Warning()
    CALL Cleanup()
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
END MODULE ERROR_Mod
