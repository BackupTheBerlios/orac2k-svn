SUBROUTINE P_fold_i(n,o,nstart,nend,nlocal,node,nprocs)

!!$************************************************************************
!!$*   Time-stamp: <2005-01-28 16:33:23 marchi>                             *
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
  INTEGER :: o(*)
  include 'mpif.h'
  include 'mpi_size.h'
      
!!$*------------------------- LOCAL VARIABLES ----------------------------*

  INTEGER, DIMENSION (:), POINTER :: fp
  INTEGER :: nstart_b,nend_b,nlocal_b,i,ierr

!!$*----------------------- EXECUTABLE STATEMENTS ------------------------*

!!$------------------------------------------------------------------------
!!$--- Create a one dimensional array
!!$------------------------------------------------------------------------

  ALLOCATE (fp(n))
  CALL MPI_ALLREDUCE(o,fp,n,MPI_INTEGER4,MPI_SUM,MPI_COMM_WORLD,ierr)

!!$------------------------------------------------------------------------
!!$--- Copy the vector in the original ones
!!$------------------------------------------------------------------------

  DO i=1,n
     o(i)=fp(i)
  END DO

  DEALLOCATE(fp)

!!$*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE P_fold_i
