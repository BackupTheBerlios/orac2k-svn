SUBROUTINE P_fold_r8(n,x,nstart,nend,nlocal,node,nprocs)

!!$************************************************************************
!!$*   Time-stamp: <2005-01-28 16:33:05 marchi>                             *
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

  INTEGER, PARAMETER :: m_tot=94000
  REAL(8)  :: fp(m_tot)
  INTEGER :: nstart_b,nend_b,nlocal_b,i,ierr

!!$*----------------------- EXECUTABLE STATEMENTS ------------------------*

!!$------------------------------------------------------------------------
!!$--- Create a one dimensional array
!!$------------------------------------------------------------------------

  CALL MPI_ALLREDUCE(x,fp,n,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

!!$------------------------------------------------------------------------
!!$--- Copy the vector in the original ones
!!$------------------------------------------------------------------------

  DO i=1,n
     x(i)=fp(i)
  END DO

!!$*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE P_fold_r8
