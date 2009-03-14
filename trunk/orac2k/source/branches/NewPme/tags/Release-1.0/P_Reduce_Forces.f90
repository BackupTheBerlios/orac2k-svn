SUBROUTINE P_Reduce_Forces(x,y,z,nstart_a,nend_a,nlocal_a,node_a,nprocs_a)

!!$************************************************************************
!!$*   Time-stamp: <2005-01-28 16:34:00 marchi>                             *
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
      
  REAL(8) :: x(*),y(*),z(*)
  INTEGER, OPTIONAL :: nstart_a,nend_a,nlocal_a,node_a,nprocs_a
      
!!$*------------------------- LOCAL VARIABLES ----------------------------*

  include 'mpif.h'
  include 'mpi_size.h'
  INTEGER :: ierr

  REAL(8), DIMENSION(:), ALLOCATABLE, SAVE :: fp,fp_out
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: starts,ends,locals
  INTEGER, SAVE :: node,nprocs,nstart,nend,nlocal
  INTEGER :: n,ntot,i,k,ii,fp_length
  
!!$*----------------------- EXECUTABLE STATEMENTS ------------------------*

!!$------------------------------------------------------------------------
!!$--- Create a one dimensional array
!!$------------------------------------------------------------------------

  IF(PRESENT(nstart_a)) THEN
     node=node_a
     nprocs=nprocs_a
     nstart=nstart_a
     nend=nend_a
     nlocal=nlocal_a

     ALLOCATE(starts(nprocs),ends(nprocs),locals(nprocs))
     CALL MPI_ALLGATHER(nstart,1,MPI_INTEGER4,starts,1,MPI_INTEGER4&
          &,MPI_COMM_WORLD,ierr)
     CALL MPI_ALLGATHER(nend,1,MPI_INTEGER4,ends,1,MPI_INTEGER4&
          &,MPI_COMM_WORLD,ierr)
     CALL MPI_ALLGATHER(nlocal,1,MPI_INTEGER4,locals,1,MPI_INTEGER4&
          &,MPI_COMM_WORLD,ierr)
     n=0
     DO i=1,nprocs
        IF(n < locals(i)) n=locals(i)
     END DO
     ntot=n*3
     ALLOCATE(fp(ntot),fp_out(ntot))
     RETURN
  END IF
  DO k=0,nprocs-1
     DO ii=1,locals(k+1)
        i=starts(k+1)-1+ii
        fp((ii-1)*3+1)=x(i)
        fp((ii-1)*3+2)=y(i)
        fp((ii-1)*3+3)=z(i)
     END DO
     fp_length=locals(k+1)*3
     CALL MPI_REDUCE(fp, fp_out, fp_length, MPI_REAL8, MPI_SUM, k&
          &, MPI_COMM_WORLD, ierr)

!!$------------------------------------------------------------------------
!!$--- Copy the vector in the original ones
!!$------------------------------------------------------------------------

     IF(k == node) THEN
        DO ii=1,nlocal
           i=nstart-1+ii
           x(i)=fp_out((ii-1)*3+1)
           y(i)=fp_out((ii-1)*3+2)
           z(i)=fp_out((ii-1)*3+3)
        END DO
     END IF
  END DO

!!$*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE P_Reduce_Forces
