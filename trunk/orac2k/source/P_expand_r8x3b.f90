SUBROUTINE P_expand_r8x3b(x,y,z,ptr,nstart,nend,nlocal,node,nprocs)

!!$***********************************************************************
!!$   Time-stamp: <04/11/10 14:41:15 marchi>                           *
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
  
!!$------------------------- LOCAL VARIABLES ----------------------------*
  
  REAL(8), DIMENSION (:), ALLOCATABLE :: fp,fp_out
  INTEGER :: n,nstart_b,nend_b,nlocal_b,i,ierr,nato
  INTEGER :: locals(256),displ(256)


!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

!!$------------------------------------------------------------------------
!!$--- Create a one dimensional array
!!$------------------------------------------------------------------------

  nstart_b=nstart+ptr-1
  nend_b=nend+ptr-1

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

END SUBROUTINE P_expand_r8x3b
