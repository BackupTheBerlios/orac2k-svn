SUBROUTINE P_merge_vecr8(x,n)

!!$************************************************************************
!!$*   Time-stamp: <2005-01-28 16:35:54 marchi>                             *
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
      
  REAL(8) :: x(*)
  INTEGER :: n
  include 'mpif.h'
  include 'mpi_size.h'

      
!!$*------------------------- LOCAL VARIABLES ----------------------------*

  REAL(8), DIMENSION (:), POINTER :: fp
  INTEGER :: ierr

!!$*----------------------- EXECUTABLE STATEMENTS ------------------------*

!!$------------------------------------------------------------------------
!!$--- Reduce 
!!$------------------------------------------------------------------------


  ALLOCATE (fp(n))
  CALL MPI_ALLREDUCE(x,fp,n,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

!!$------------------------------------------------------------------------
!!$--- Copy the vector in the original ones
!!$------------------------------------------------------------------------
  
  x(1:n)=fp

  DEALLOCATE(fp)

!!$*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE P_merge_vecr8
