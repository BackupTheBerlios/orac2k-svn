MODULE Strings

!!$***********************************************************************
!!$   Time-stamp: <2007-01-08 12:41:56 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Mon Nov 13 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

  USE Errors, ONLY: Add_Errors=> Add, errmsg_f, Print_Errors
  USE CONSTANTS, ONLY: max_pars, max_char, Comms
  USE TYPES
  USE STRPAK
CONTAINS
  FUNCTION My_FAM(str1,str2)
    IMPLICIT NONE 
    CHARACTER(len=*) :: str1,str2
    LOGICAL :: My_Fam
    INTEGER :: j1,j2,iflag

    My_FAM=.FALSE.
    CALL STR_FAM(str1,str2,j1,j2,iflag)
    IF(iflag == 0) MY_Fam=.TRUE.
  END FUNCTION My_FAM
  FUNCTION My_FXM(str1,str2)
    IMPLICIT NONE 
    CHARACTER(len=*) :: str1,str2
    LOGICAL :: My_Fxm
    INTEGER :: j1,j2,iflag

    My_FXM=.FALSE.
    CALL STR_FXM(str1,str2,j1,j2,iflag)
    IF(iflag == 0) MY_Fxm=.TRUE.
  END FUNCTION My_FXM
  FUNCTION ST_Concat(i,strngs)
    IMPLICIT NONE 
    CHARACTER(len=max_char) :: ST_Concat
    INTEGER :: i
    CHARACTER(len=max_pars), DIMENSION(:) :: strngs

    INTEGER :: n,m,pad
    CHARACTER(len=max_char) :: line,aux 

    CALL STR_FILL(' ',line)
    
    pad=1
    m=0
    DO n=i,SIZE(strngs)
       m=m+LEN_TRIM(strngs(n))+pad
    END DO
    IF(m > LEN(line)) THEN
       WRITE(aux,'(i3)') LEN(line)
       errmsg_f='Cannot contatenate strings. Final string length larger&
            & than hard bound '//TRIM(aux)
       CALL Add_Errors(-1,errmsg_f)
       CALL Print_Errors()
    END IF
    DO n=i,SIZE(strngs)
       line=TRIM(line)//' '//TRIM(strngs(n))
    END DO
    ST_Concat=line
  END FUNCTION ST_Concat
  FUNCTION ST_Clean_Line(line)
    IMPLICIT NONE 
    LOGICAL :: ST_Clean_Line
    CHARACTER(len=max_char) :: line
    INTEGER :: i

    ST_Clean_Line=.FALSE.
    DO i=1,SIZE(Comms)
       CALL STR_TRIM(Comms(i),line)
    END DO
    IF(LEN_TRIM(line) == 0) ST_Clean_Line=.TRUE.
  END FUNCTION ST_Clean_Line
END MODULE Strings
