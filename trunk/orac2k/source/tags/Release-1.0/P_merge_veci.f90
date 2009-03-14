SUBROUTINE P_merge_veci(o,n)

!!$************************************************************************
!!$*   Time-stamp: <2005-01-28 16:36:05 marchi>                             *
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
      
  INTEGER :: o(*)
  INTEGER :: n
  include 'mpif.h'
  include 'mpi_size.h'

      
!!$*------------------------- LOCAL VARIABLES ----------------------------*

  INTEGER, DIMENSION (:), POINTER :: fp
  INTEGER :: ierr

!!$*----------------------- EXECUTABLE STATEMENTS ------------------------*

!!$------------------------------------------------------------------------
!!$--- Reduce 
!!$------------------------------------------------------------------------


  ALLOCATE (fp(n))
  CALL MPI_ALLREDUCE(o,fp,n,MPI_INTEGER4,MPI_SUM,MPI_COMM_WORLD,ierr)

!!$------------------------------------------------------------------------
!!$--- Copy the vector in the original ones
!!$------------------------------------------------------------------------
  
  o(1:n)=fp
  DEALLOCATE(fp)

!!$*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE P_merge_veci
