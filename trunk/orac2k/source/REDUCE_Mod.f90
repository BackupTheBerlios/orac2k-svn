MODULE REDUCE
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors,errmsg_f
  PRIVATE
  PUBLIC Init, Reduce_Forces, Reduce_Phi, No_of_Calls
  REAL(8), DIMENSION(:), ALLOCATABLE, SAVE :: fp,fp_out
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: starts,ends,locals
  INTEGER, SAVE :: node,nprocs,nstart,nend,nlocal
  INTEGER, SAVE :: No_of_Calls=0
CONTAINS
  SUBROUTINE Init(nstart_a,nend_a,nlocal_a,node_a,nprocs_a)
    INTEGER :: nstart_a,nend_a,nlocal_a,node_a,nprocs_a
#ifdef PARALLEL    
    include 'mpif.h'
    include 'mpi_size.h'
#endif
    No_of_Calls=No_of_Calls+1
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
  END SUBROUTINE Init
  SUBROUTINE Reduce_Forces(x,y,z,nstart_a,nend_a,nlocal_a,node_a,nprocs_a)

!!$************************************************************************
!!$*   Time-stamp: <2007-11-23 10:41:47 marchi>                             *
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

#ifdef PARALLEL
    
    include 'mpif.h'
    include 'mpi_size.h'
#endif

    INTEGER :: ierr    
    INTEGER :: n,ntot,i,k,ii,fp_length
    
!!$*----------------------- EXECUTABLE STATEMENTS ------------------------*
    
!!$------------------------------------------------------------------------
!!$--- Create a one dimensional array
!!$------------------------------------------------------------------------
    
    IF(PRESENT(nstart_a)) THEN
       CALL Init(nstart_a,nend_a,nlocal_a,node_a,nprocs_a)
       RETURN
    END IF
    IF(No_of_Calls < 1) THEN
       errmsg_f='Need Initialisation to reduce three vectors.'
       CALL Add_Errors(-1,errmsg_f)
       CALL Print_Errors()
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
    
  END SUBROUTINE Reduce_Forces
  SUBROUTINE Reduce_Phi(x)
    REAL(8) :: x(*)
    include 'mpif.h'
    include 'mpi_size.h'
    INTEGER :: ierr
    
    INTEGER :: n,ntot,i,k,ii,fp_length


    IF(No_of_Calls < 1) THEN
       errmsg_f='Need Initialisation to reduce vector.'
       CALL Add_Errors(-1,errmsg_f)
       CALL Print_Errors()
    END IF
    DO k=0,nprocs-1
       DO ii=1,locals(k+1)
          i=starts(k+1)-1+ii
          fp(ii)=x(i)
       END DO
       fp_length=locals(k+1)
       CALL MPI_REDUCE(fp, fp_out, fp_length, MPI_REAL8, MPI_SUM, k&
            &, MPI_COMM_WORLD, ierr)
       
!!$------------------------------------------------------------------------
!!$--- Copy the vector in the original ones
!!$------------------------------------------------------------------------
       
       IF(k == node) THEN
          DO ii=1,nlocal
             i=nstart-1+ii
             x(i)=fp_out(ii)
          END DO
       END IF
    END DO
    
  END SUBROUTINE Reduce_Phi
END MODULE REDUCE
