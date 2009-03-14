      SUBROUTINE P_barrier

************************************************************************
*   Time-stamp: <04/12/18 14:40:04 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Mon Feb 15 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------- VARIABLES IN COMMON --------------------------*

#ifdef MPI
      INCLUDE 'mpif.h'

*------------------------- LOCAL VARIABLES ----------------------------*

#if defined T3E
      INTEGER*8 comm,ierr
#else
      INTEGER comm,ierr
#endif

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      comm=MPI_COMM_WORLD
      CALL MPI_BARRIER(comm,ierr)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*
#endif
      RETURN
      END
