      SUBROUTINE P_get_iret(iret,node,nprocs,ncube,nbyte)

************************************************************************
*   Time-stamp: <2005-01-28 16:36:18 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Jul 11 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER iret,node,nprocs,ncube,nbyte
      INCLUDE 'mpif.h'
      include 'mpi_size.h'


*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,ierr
      INTEGER irets(512)
      
*----------------------- EXECUTABLE STATEMENTS ------------------------*


      CALL MPI_ALLGATHER(iret,1,MPI_INTEGER4,irets,1,MPI_INTEGER4
     &     ,MPI_COMM_WORLD,ierr)
      iret=0
      DO i=1,nprocs
         IF(irets(i) .NE. 0) iret=irets(i)
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
