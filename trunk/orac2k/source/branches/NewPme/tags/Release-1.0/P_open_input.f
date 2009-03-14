      SUBROUTINE P_open_input(myid,numprocs,ncube,nbyte)

************************************************************************
*   Time-stamp: <2006-02-05 16:27:23 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed Jul  8 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER myid,numprocs,ncube,nbyte

*----------------------- VARIABLES IN COMMON --------------------------*

#ifdef MPI      
      INCLUDE 'mpif.h'
      include 'mpi_size.h'

      INTEGER status(MPI_STATUS_SIZE)

*------------------------- LOCAL VARIABLES ----------------------------*

      CHARACTER*80 line,errmsg
      INTEGER iret,i
      INTEGER ierr

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      iret=1
      IF(myid .EQ. 0) THEN
         OPEN (unit=99,file="./.ORAC.INPUT")
1000     READ(5,'(a80)',END=2000) line
         WRITE(*,'(a80)') line
         WRITE(99,'(a80)') line
         GOTO 1000
2000     CONTINUE
         iret=0
         CLOSE(99)
      END IF
      CALL P_Barrier

      CALL MPI_BCAST(iret,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
      IF(iret .NE. 0) THEN
         call MPI_FINALIZE(ierr)
         errmsg=
     &        ' Badly sent or received msg. Abort. '
         CALL xerror(errmsg,80,1,2)
      END IF
#endif

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
