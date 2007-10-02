MODULE STRINGS_Mod

!!$***********************************************************************
!!$   Time-stamp: <2006-11-24 12:43:16 marchi>                           *
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

    INTEGER :: n
    CHARACTER(len=max_char) :: line
    CALL STR_FILL(' ',line)
    WRITE(line,'(10a7)') (ADJUSTL(strngs(n)),n=i,SIZE(strngs))
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
END MODULE STRINGS_Mod
