SUBROUTINE P_isend_i(buf,count,dest)

!!$***********************************************************************
!!$   Time-stamp: <04/11/11 16:02:25 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Nov 11 2004 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*

  IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

  INTEGER :: buf(*),count,dest

!!$----------------------- VARIABLES IN COMMON --------------------------*

  INCLUDE 'mpif.h'

!!$------------------------- LOCAL VARIABLES ----------------------------*


!!$----------------------- EXECUTABLE STATEMENTS ------------------------*


  CALL MPI_ISEND(buf, count, MPI_INTEGER4, dest, 0, MPI_COMM_WORLD, ierr)

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE P_isend_i
