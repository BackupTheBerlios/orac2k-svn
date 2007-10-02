SUBROUTINE P_expand_i(o,nstart,nend,nlocal,node,nprocs)

!!$***********************************************************************
!!$   Time-stamp: <2005-01-28 16:32:31 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Sat Mar 10 2001 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program orac ----*


!!$======================== DECLARATIONS ================================*

  IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

  INTEGER :: nstart,nend,nlocal,node,nprocs,ncube
  INTEGER :: o(*)
  include 'mpif.h'
  include 'mpi_size.h'
  
!!$------------------------- LOCAL VARIABLES ----------------------------*
  
  INTEGER :: n,nstart_b,nend_b,nlocal_b,i,ierr,nato
  INTEGER :: locals(256),displ(256)


!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

!!$------------------------------------------------------------------------
!!$--- Create a one dimensional array
!!$------------------------------------------------------------------------

  CALL MPI_ALLGATHER(nlocal,1,MPI_INTEGER4,locals,1 &
       ,MPI_INTEGER4,MPI_COMM_WORLD,ierr)
  displ(1)=0
  DO i=2,nprocs
     displ(i)=displ(i-1)+locals(i-1)
  END DO
  CALL MPI_ALLGATHERV(o(nstart),nlocal,MPI_INTEGER4,o,locals,displ&
       &,MPI_INTEGER4,MPI_COMM_WORLD,ierr) 

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE P_expand_i
