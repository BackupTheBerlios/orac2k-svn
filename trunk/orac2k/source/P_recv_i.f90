SUBROUTINE P_recv_i(buf,count,source,tag)

!!$***********************************************************************
!!$   Time-stamp: <04/11/11 16:09:21 marchi>                           *
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

 INTEGER :: buf(*),count,source,tag


!!$----------------------- VARIABLES IN COMMON --------------------------*

  INCLUDE 'mpif.h'

!!$------------------------- LOCAL VARIABLES ----------------------------*

  INTEGER :: ierr

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*


  CALL MPI_RECV(buf, count, MPI_INTEGER4, source, tag, MPI_COMM_WORLD, ierr)


!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE P_recv_i
