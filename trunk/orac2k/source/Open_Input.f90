SUBROUTINE Open_Input(kread)

!!$***********************************************************************
!!$   Time-stamp: <2006-02-08 14:31:26 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Sun Feb  5 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program  ----*


!!$======================== DECLARATIONS ================================*

  IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

  INTEGER :: kread

!!$------------------------- LOCAL VARIABLES ----------------------------*

  INTEGER :: IARGC,nargs,nchars
  CHARACTER(80) :: buffer,fileinput,errmsg,command
  LOGICAL, SAVE :: ok=.FALSE.

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*


  nargs=IARGC()
  IF(nargs >= 1) THEN
     ok=.TRUE.
     CALL GETARG(1, buffer)
     fileinput=TRIM(buffer)
     OPEN(unit=kread,file=fileinput,form='FORMATTED')
  END IF

  IF(.NOT. ok) THEN
     CALL GETARG(0, buffer)
     command = TRIM(buffer)
     nchars=LEN_TRIM(buffer)
     errmsg='Usage: '// command(1:nchars) //' input [ > output ]'
     CALL xerror(errmsg,80,1,2)
  END IF

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE Open_Input
