SUBROUTINE P_merge_i(o)

!!$************************************************************************
!!$*   Time-stamp: <2005-01-28 16:35:34 marchi>                             *
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
      
  INTEGER :: o
  include 'mpif.h'
  include 'mpi_size.h'
      
!!$*------------------------- LOCAL VARIABLES ----------------------------*

  INTEGER :: fp
  INTEGER :: n
  INTEGER :: ierr

!!$*----------------------- EXECUTABLE STATEMENTS ------------------------*

!!$------------------------------------------------------------------------
!!$--- Reduce 
!!$------------------------------------------------------------------------

  n=1
  CALL MPI_ALLREDUCE(o,fp,n,MPI_INTEGER4,MPI_SUM,MPI_COMM_WORLD,ierr)

!!$------------------------------------------------------------------------
!!$--- Copy the vector in the original ones
!!$------------------------------------------------------------------------
  
  o=fp

!!$*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE P_merge_i
