MODULE INPUT_Mod

!!$***********************************************************************
!!$   Time-stamp: <2005-10-17 18:57:00 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Sun Oct 16 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*

  USE Xerror_Mod

  CHARACTER(22), SAVE :: err_open='OPEN keyword not found'
  CHARACTER(37), DIMENSION(3), SAVE :: err_args=(/ &
       & 'Number of arguments must be at least '   &
       &,'Number of arguments must not exceed  '   &
       &,'Number of arguments must be only     '/)
  CHARACTER(20), SAVE :: err_end='...or missing &END'
  CHARACTER(27), DIMENSION(4), SAVE :: err_unr=&
       &(/'Unrecognized command  ---> '&
       &, 'Unrecognized subcommand -> '&
       &, 'Unrecognized keyword ----> '&
       &, 'UNSUPPORTED  COMMAND ----> '/)
  CHARACTER(15), SAVE :: err_fnf=' file not found'
  CHARACTER(80), SAVE :: line,strngs(40)
  CHARACTER(1),  DIMENSION(2), SAVE :: sep=(/' ',','/)&
       &,comm=(/'(',')'/)
  CHARACTER(8) :: fmt
  INTERFACE Read_String
     MODULE PROCEDURE Read_string_i4
     MODULE PROCEDURE Read_string_r8
  END INTERFACE
  CONTAINS
    SUBROUTINE Read_String_i4(string,out)
      IMPLICIT NONE 
      CHARACTER(80) :: string
      INTEGER :: ivalue,out
      CALL fndfmt(1,string,fmt)
      READ(string,fmt,err=20) ivalue
      out=ivalue
      RETURN

20    CALL abort_now('internal reading error: wrong format?? TAB character??')

    END SUBROUTINE Read_String_i4
    SUBROUTINE Read_String_r8(string,out)
      IMPLICIT NONE 
      CHARACTER(80) :: string
      REAL(8) :: rvalue,out
      CALL fndfmt(2,string,fmt)
      READ(string,fmt,err=20) rvalue
      out=rvalue
      RETURN

20    CALL abort_now('internal reading error: wrong format?? TAB character??')

    END SUBROUTINE Read_String_r8
    SUBROUTINE parser(line,strings,nword)
      IMPLICIT NONE 
      INTEGER :: nword_max
      PARAMETER (nword_max=40)
      CHARACTER(80) :: strings(nword_Max),line
      INTEGER :: nword
      
      INTEGER :: iret
      CHARACTER(80) :: errmsg
      
      CALL parse(line,sep,2,comm,strings,nword_max,nword,iret,errmsg)
      IF(iret /=0) CALL abort_now(errmsg)
      RETURN
    END SUBROUTINE parser
END MODULE INPUT_Mod
