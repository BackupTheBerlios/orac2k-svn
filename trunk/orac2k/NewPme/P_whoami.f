      SUBROUTINE P_whoami(node,nprocs,ncube,nbyte)

************************************************************************
*   Time-stamp: <2009-03-12 10:16:03 marchi>                             *
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

      INTEGER node,nprocs,ncube,nbyte

*----------------------- VARIABLES IN COMMON --------------------------*

#ifdef MPI
      INCLUDE 'mpif.h'
      INTEGER ierr,status(MPI_STATUS_SIZE)

*------------------------- LOCAL VARIABLES ----------------------------*

      CHARACTER*80 errmsg
      INTEGER nodea,ierra,nprocsa,comm,iret

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      ierra=0
      call MPI_INIT( ierra )
      comm=MPI_COMM_WORLD
      call MPI_COMM_RANK( comm, nodea, ierra )
      call MPI_COMM_SIZE( comm, nprocsa, ierra )
      
      node=nodea
      nprocs=nprocsa
      ncube = 0
10    IF (2**ncube .LT. nprocs) THEN
         ncube = ncube + 1
         GOTO 10
      END IF
#endif
#ifdef SHMEM
      IF(nprocs .NE. 2**ncube) THEN
         iret=1
         errmsg=' Number of processors *MUST* be a power of two '
     &        / /'in SHMEM mode.'
         CALL P_get_errmsg(iret,errmsg,80,node,nprocs,ncube,nbyte)
      END IF
#endif

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
