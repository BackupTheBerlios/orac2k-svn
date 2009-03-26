SUBROUTINE P_expand_r8x3(x,y,z,nstart,nend,nlocal,node,nprocs,ptr)

!!$***********************************************************************
!!$   Time-stamp: <2009-03-25 22:27:53 marchi>                           *
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

  INTEGER :: nstart,nend,nlocal,node,nprocs,ncube,ptr
  REAL(8) :: x(*),y(*),z(*)
  include 'mpif.h'
  include 'mpi_size.h'
!!$------------------------- LOCAL VARIABLES ----------------------------*
  
  INTEGER :: n,nstart_b,nend_b,nlocal_b,i,nato
  INTEGER :: locals(256),displ(256)
  INTEGER :: ierr

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

!!$------------------------------------------------------------------------
!!$--- Create a one dimensional array
!!$------------------------------------------------------------------------

  If(Nprocs == 1) Return
  CALL MPI_ALLGATHER(nlocal,1,MPI_INTEGER4,locals,1 &
       ,MPI_INTEGER4,MPI_COMM_WORLD,ierr)
  displ(1)=ptr-1
  DO i=2,nprocs
     displ(i)=displ(i-1)+locals(i-1)
  END DO
  CALL MPI_ALLGATHERV(x(nstart),nlocal,MPI_REAL8,x,locals,displ&
       &,MPI_REAL8,MPI_COMM_WORLD,ierr) 
  CALL MPI_ALLGATHERV(y(nstart),nlocal,MPI_REAL8,y,locals,displ&
       &,MPI_REAL8,MPI_COMM_WORLD,ierr) 
  CALL MPI_ALLGATHERV(z(nstart),nlocal,MPI_REAL8,z,locals,displ&
       &,MPI_REAL8,MPI_COMM_WORLD,ierr) 

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE P_expand_r8x3
