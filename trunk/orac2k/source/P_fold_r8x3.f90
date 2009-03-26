SUBROUTINE P_fold_r8x3(n,x,y,z,nstart,nend,nlocal,node,nprocs)

!!$************************************************************************
!!$*   Time-stamp: <2009-03-25 22:34:57 marchi>                             *
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
  REAL(8) :: x(*),y(*),z(*)
  include 'mpif.h'
  include 'mpi_size.h'
      
!!$*------------------------- LOCAL VARIABLES ----------------------------*

  INTEGER, PARAMETER :: m_tot=94000
  REAL(8), DIMENSION(:), ALLOCATABLE  :: fp,fp_out
  INTEGER :: nstart_b,nend_b,nlocal_b,ntot,i
  INTEGER :: ierr

!!$*----------------------- EXECUTABLE STATEMENTS ------------------------*

!!$------------------------------------------------------------------------
!!$--- Create a one dimensional array
!!$------------------------------------------------------------------------

  If(nprocs == 1) Return
  ntot=n*3
  ALLOCATE(fp(ntot),fp_out(ntot))
  DO i=1,n
     fp((i-1)*3+1)=x(i)
     fp((i-1)*3+2)=y(i)
     fp((i-1)*3+3)=z(i)
  END DO

  CALL MPI_ALLREDUCE(fp,fp_out,ntot,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

!!$------------------------------------------------------------------------
!!$--- Copy the vector in the original ones
!!$------------------------------------------------------------------------

  DO i=1,n
     x(i)=fp_out((i-1)*3+1)
     y(i)=fp_out((i-1)*3+2)
     z(i)=fp_out((i-1)*3+3)
  END DO

  DEALLOCATE(fp,fp_out)
!!$*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE P_fold_r8x3
