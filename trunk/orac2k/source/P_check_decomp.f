      SUBROUTINE P_check_decomp(nstart,nend,node,nbyte,nprocs,tag,iret
     &     ,errmsg)

************************************************************************
*   Time-stamp: <2005-01-28 19:36:32 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Feb 21 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nstart,nend,node,nbyte,nprocs,iret
      CHARACTER*80 errmsg,tag
      INCLUDE 'mpif.h'
      INCLUDE 'mpi_size.h'

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER temp,i
      INTEGER ierr
      INTEGER starts(512),ends(512)
      LOGICAL ok

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      ok=.TRUE.
      CALL MPI_ALLGATHER(nstart,1,MPI_INTEGER4,starts,1,MPI_INTEGER4
     &     ,MPI_COMM_WORLD,ierr)
      CALL MPI_ALLGATHER(nend,1,MPI_INTEGER4,ends,1,MPI_INTEGER4
     &     ,MPI_COMM_WORLD,ierr)
      IF(nprocs .NE. 1) THEN
         DO i=2,nprocs
            temp=ends(i-1)+1
            IF(temp .NE. starts(i)) ok=.FALSE.
         END DO
      END IF
      IF(.NOT. ok) THEN
         iret=1
         WRITE(errmsg,'(''Atom decomposition check failed:''
     &        ,'' Decomposition tagged '',a12,'' overlaps.'')') tag
         WRITE(*,'(''Atom decomposition check failed:''
     &        ,'' Decomposition tagged '',a12,'' overlaps.'')') tag
      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
