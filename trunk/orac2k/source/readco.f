      SUBROUTINE readco(restart_in,restart_out,restart_read
     &     ,restart_write,co,oc,volume)

************************************************************************
*   Time-stamp: <99/10/11 18:31:57 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Jul  2 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8  co(3,3),oc(3,3),volume
      CHARACTER*80 restart_in,restart_out
      LOGICAL restart_read,restart_write

*----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'unit.h'

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,kdump
      REAL*8  oc1(3,3)
      CHARACTER*80 file

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      kdump=kdump_in
      file=restart_in
      IF((.NOT. restart_read) .AND. restart_write) THEN
         kdump=kdump_out
         file=restart_out
      END IF
      OPEN(unit=kdump,file=file,form='UNFORMATTED',status='OLD')      
      REWIND kdump
      READ(kdump) ((co(i,j),j=1,3),i=1,3),((oc(i,j),j=1,3),i=1,3)
      CALL matinv(3,3,co,oc1,volume)
      volume=volume*boxl**3
      CLOSE(kdump)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
