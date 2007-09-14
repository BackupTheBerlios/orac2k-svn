Module xerror_mod

  TYPE fatal_err
     CHARACTER(74), DIMENSION (6) :: top= &
          (/ '**************************************************************************' &
          ,'*                                                                        *' & 
          ,'*-- Fatal ----E R R O R ----- E R R E U R -----E R R O R E !!!----Fatal--*' &
          ,'*                                                                        *' &
          ,'*-- The following errors were found:          ---------------------------*' &
          ,'*                                                                        *' /)
     CHARACTER(74), DIMENSION (:), POINTER :: body
     CHARACTER(74), DIMENSION(4) :: bottom= &
          (/ '*                                                                        *' &
          ,'*-- Fatal ----E R R O R ----- E R R E U R -----E R R O R E ***----Fatal--*' &
          ,'*                                                                        *' &
          ,'**************************************************************************' /)
  END TYPE fatal_err

  TYPE warning_err
     CHARACTER(74), DIMENSION (6) :: top= &
          (/ '**************************************************************************' &
          ,'*                                                                        *' & 
          ,'*-- Warning---E R R O R ----- E R R E U R -----E R R O R E !!!---Warning-*' &
          ,'*                                                                        *' &
          ,'*-- The following errors were found:          ---------------------------*' &
          ,'*                                                                        *' /)
     CHARACTER(74), DIMENSION (:), POINTER :: body
     CHARACTER(74), DIMENSION(4) :: bottom= &
          (/ '*                                                                        *' &
          ,'*-- Warning---E R R O R ----- E R R E U R -----E R R O R E ***---Warning-*' &
          ,'*                                                                        *' &
          ,'**************************************************************************' /)
  END TYPE warning_err

  CHARACTER(LEN=45), SAVE :: err_arg_no='.. Arguments to command must be at least '



  !***********************************************************************
  !   Time-stamp: <1999-12-06 10:51:27 marchi>                           *
  !                                                                      *
  !  This module handles error conditions and abort on request           *
  !                                                                      *
  !======================================================================*
  !                                                                      *
  !              Author:  Massimo Marchi                                 *
  !              CEA/Centre d'Etudes Saclay, FRANCE                      *
  !                                                                      *
  !              - Sun Dec  5 1999 -                                     *
  !                                                                      *
  !***********************************************************************
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
  SUBROUTINE write_error(msg,write_msgs,type)
    IMPLICIT NONE
    CHARACTER(*) :: msg
    CHARACTER(*), OPTIONAL :: type
    LOGICAL, OPTIONAL :: write_msgs
  END SUBROUTINE write_error
END MODULE xerror_mod
