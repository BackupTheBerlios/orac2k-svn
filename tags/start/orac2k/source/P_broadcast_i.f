      SUBROUTINE P_broadcast_i(data,n)

************************************************************************
*   Time-stamp: <2005-01-28 19:36:21 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Jul  9 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER  data(*)
      INTEGER n
      INCLUDE 'mpif.h'
      INCLUDE 'mpi_size.h'

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER ierr

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      CALL MPI_BCAST(data,n,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
