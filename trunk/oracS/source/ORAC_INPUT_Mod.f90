MODULE ORAC_INPUT_Mod
  
!!$***********************************************************************
!!$   Time-stamp: <2006-11-24 12:56:53 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Tue Nov 21 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This MODULE is part of the program ORAC ----*

  USE XERROR, ONLY: abort_now
  USE STRPAK
  USE TYPES
  USE STRINGS_Mod, ONLY: MY_Fam, ST_Clean_Line
  IMPLICIT none
  PRIVATE
  PUBLIC Store,Modify,Single_String,kinput,file,Options,Input_String
  INTEGER, SAVE :: kinput
  CHARACTER(len=max_char), SAVE :: file
  CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE, SAVE :: data
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: ind
  CHARACTER(len=max_data), SAVE :: Input_String
  CHARACTER(len=120) :: errmsg
  INTEGER, EXTERNAL :: iargc
  LOGICAL, SAVE :: prints=.FALSE.
CONTAINS
  SUBROUTINE Store
    CHARACTER(len=max_char) :: line,ui
  
    LOGICAL :: end_of_file=.FALSE.
    INTEGER :: io,c,iopt,iflag
    
    CALL CHANNEL(io)
    kinput=io
    OPEN(unit=kinput, file=file, form='FORMATTED')
    c=0
    DO WHILE(.NOT. end_of_file)
       READ(kinput,'(a)',IOSTAT=iopt) line
       IF(iopt /= 0) EXIT
       IF(ST_Clean_Line(line)) CYCLE
       iflag=0
       DO WHILE(iflag == 0) 
          CALL STR_RAS(CHAR(09),'   ',line,iflag)
       END DO
       IF(prints) WRITE(*,*) TRIM(line)
       c=c+1
    END DO
    ALLOCATE(data(c),ind(c))

    REWIND(kinput)
    c=0
    DO WHILE(.NOT. end_of_file)
       READ(kinput,'(a)',IOSTAT=iopt) line
       IF(ST_Clean_Line(line)) CYCLE
       IF(iopt /= 0) EXIT
       c=c+1
       iflag=0
       DO WHILE(iflag == 0) 
          CALL STR_RAS(CHAR(09),'   ',line,iflag)
       END DO
       data(c)=ADJUSTL(TRIM(line))
       ind(c)=LEN_TRIM(line)-LEN_TRIM(ADJUSTL(line))
    END DO
  END SUBROUTINE Store
  SUBROUTINE Modify
    INTEGER :: n,n_d,n_trim,m,o
    LOGICAL :: ok
    CHARACTER(len=max_char) :: dummy
    
    n_d=SIZE(data)

    DO n=1,n_d
       IF(data(n)(1:1) == '&') THEN
          IF(data(n)(1:4) == '&END') THEN
             data(n)(5:5)='}'
          ELSE
             CALL STR_SHFT(1,data(n))
             data(n)(1:1)='{'
          END IF
       END IF
    END DO

    n=0
    DO WHILE(n < n_d)
       n=n+1
       IF(MY_Fam('&',data(n))) CYCLE       
       IF(data(n)(1:3) == 'END') THEN
          data(n)(4:4)='}'
          DO o=n-1,1,-1
             IF(MY_Fam('&',data(o))) CYCLE
             dummy=data(o)
             CALL TRANUC(dummy)
             IF(dummy == data(o) .AND. (ind(n) == ind(o))) THEN
                m=o
                EXIT
             END IF
          END DO
          CALL STR_SHFT(1,data(m))
          data(m)(1:1)='{'
          DO o=m+1,n-1
             CALL STR_SHFT(1,data(o))
             data(o)(1:1)='{'
             n_trim=LEN_TRIM(data(o))
             data(o)(n_trim+1:n_trim+1)='}'
          END DO
       END IF

    END DO

    n=0
    DO WHILE(n < n_d)
       n=n+1
       IF(MY_Fam('&',data(n)) .OR. MY_Fam('{',data(n)) &
            &.OR. MY_Fam('}',data(n))) CYCLE
       CALL STR_SHFT(1,data(n))
       data(n)(1:1)='{'
       n_trim=LEN_TRIM(data(n))
       data(n)(n_trim+1:n_trim+1)='}'
    END DO

  END SUBROUTINE Modify
  SUBROUTINE Single_String
    INTEGER :: n
    Input_String=' '
    DO n=1,SIZE(data)
       Input_String=TRIM(Input_String)//' '//TRIM(data(n))
    END DO
    DEALLOCATE(data,ind)

  END SUBROUTINE Single_String
  SUBROUTINE Options
    IMPLICIT NONE 
    INTEGER :: n,ntot
    CHARACTER(len=max_char) :: line, linea
    LOGICAL :: ok

    n=0
    ntot=IARGC()
    IF(ntot == 0) THEN
       CALL GETARG(0,line)
       errmsg='Usage: '//TRIM(line)//' [-V] input > output '
       CALL abort_now(errmsg)
    END IF
    ok=.FALSE.
    DO WHILE(n < ntot)
       n=n+1
       CALL GETARG(n,line)
       linea=TRIM(line)
       IF(linea(1:1) == '-') THEN
          line=linea(2:)
          IF(line == 'V') THEN
             prints=.TRUE.
          ELSE
             errmsg='Unknown input option. Found '//TRIM(linea)
             CALL abort_now(errmsg)
          END IF
       ELSE
          file=linea
          ok=.TRUE.
          EXIT
       END IF
    END DO
    IF(.NOT. ok) THEN
       CALL GETARG(0,line)
       errmsg='Usage: '//TRIM(line)//' [-V] input > output '
       CALL abort_now(errmsg)
    END IF
  END SUBROUTINE Options

END MODULE ORAC_INPUT_Mod
