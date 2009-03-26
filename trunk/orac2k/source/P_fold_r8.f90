SUBROUTINE P_fold_r8(n,x,nstart,nend,nlocal,node,nprocs)

!!$************************************************************************
!!$*   Time-stamp: <2009-03-25 22:34:39 marchi>                             *
!!$*                                                                      *
!!$*                                                                      *
!!$*                                                                      *
!!$*======================================================================*
!!$*                                                                      *
!!$*              Author:  Massimo Marchi                                 *
!!$*              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$*                                                                      *
!!$*              - Mon Feb 22 1999 -                                     *
!!$*                                                                      *
!!$************************************************************************

!!$*---- This subroutine is part of the program ORAC ----*


!!$*======================== DECLARATIONS ================================*

  IMPLICIT none
  
!!$*----------------------------- ARGUMENTS ------------------------------*
      
  INTEGER :: n,nstart,nend,nlocal,node,nprocs
  REAL(8) :: x(*)
  include 'mpif.h'
  include 'mpi_size.h'
      
!!$*------------------------- LOCAL VARIABLES ----------------------------*

  REAL(8), Allocatable  :: fp(:)
  INTEGER :: nstart_b,nend_b,nlocal_b,i,ierr

!!$*----------------------- EXECUTABLE STATEMENTS ------------------------*

!!$------------------------------------------------------------------------
!!$--- Create a one dimensional array
!!$------------------------------------------------------------------------

  If(nprocs == 1) Return
  Allocate(fp(n))
  CALL MPI_ALLREDUCE(x,fp,n,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

!!$------------------------------------------------------------------------
!!$--- Copy the vector in the original ones
!!$------------------------------------------------------------------------

  DO i=1,n
     x(i)=fp(i)
  END DO

!!$*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE P_fold_r8
