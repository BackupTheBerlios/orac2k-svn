      SUBROUTINE find_length_run(stop_anl)

************************************************************************
*   Time-stamp: <2008-03-10 17:04:24 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed Nov 26 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER stop_anl

*----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'unit.h'

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER rec,unit
      LOGICAL ok
      REAL*4 faux
      LOGICAL near0

*----------------------- EXECUTABLE STATEMENTS ------------------------*

*=======================================================================
*----- Read box and timestep -------------------------------------------
*=======================================================================

      unit=kwrite_dump_i(1)
      ok=.TRUE.
      rec=0
100   CONTINUE
         rec=rec+1
         READ(unit=unit,rec=rec,ERR=200) faux
         GOTO 100
200   CONTINUE 
      stop_anl=rec-1

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
