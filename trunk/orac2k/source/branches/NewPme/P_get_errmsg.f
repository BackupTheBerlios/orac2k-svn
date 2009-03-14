      SUBROUTINE P_get_errmsg(iret,errmsg,nerr,node,nprocs,ncube,nbyte)

************************************************************************
*   Time-stamp: <2005-01-28 16:33:39 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Fri Feb 19 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER iret,nerr,node,nprocs,ncube,nbyte
      CHARACTER*1 errmsg(80)
      INCLUDE 'mpif.h'
      include 'mpi_size.h'

*------------------------- LOCAL VARIABLES ----------------------------*


      INTEGER ierr,nchar,i,irets(512),ia
      CHARACTER*1 errs(80,512)
      LOGICAL ok

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      CALL MPI_ALLGATHER(iret,1,MPI_INTEGER4,irets,1,MPI_INTEGER4
     &     ,MPI_COMM_WORLD,ierr)

      ok=.TRUE.
      DO i=1,nprocs
         IF(irets(i) .NE. 0) THEN
            ia=i-1
            ok=.FALSE.
         END IF
      END DO
      IF(.NOT. ok) THEN
         CALL MPI_BCAST(errmsg,80,MPI_CHARACTER, ia, MPI_COMM_WORLD,
     &        ierr)
         CALL xerror(errmsg,80,1,2)
      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
