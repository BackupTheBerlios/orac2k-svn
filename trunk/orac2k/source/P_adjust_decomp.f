      SUBROUTINE P_adjust_decomp(nstart,nend,nlocal,node,nbyte,nprocs)

************************************************************************
*   Time-stamp: <2005-01-28 16:36:27 marchi>                             *
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

      INTEGER nstart,nend,nlocal,node,nbyte,nprocs
      INCLUDE 'mpif.h'
      include 'mpi_size.h'


*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER temp,ierr,i
      INTEGER starts(512),ends(512)
      LOGICAL ok

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      CALL MPI_ALLGATHER(nstart,1,MPI_INTEGER4,starts,1,MPI_INTEGER4
     &     ,MPI_COMM_WORLD,ierr)
      CALL MPI_ALLGATHER(nend,1,MPI_INTEGER4,ends,1,MPI_INTEGER4
     &     ,MPI_COMM_WORLD,ierr)

      IF(nprocs .NE. 1) THEN
         DO i=2,nprocs
            temp=ends(i-1)+1
            IF(temp .NE. starts(i)) starts(i)=temp
         END DO
         nstart=starts(node+1)
         nlocal=nend-nstart+1
      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
