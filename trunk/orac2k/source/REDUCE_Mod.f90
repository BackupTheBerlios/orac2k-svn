MODULE REDUCE
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors,errmsg_f
  PRIVATE
  PUBLIC Init, Reduce_Forces, Reduce_Phi, No_of_Calls
  REAL(8), DIMENSION(:), ALLOCATABLE, SAVE :: fp,fp_out
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: starts,ends,locals
  INTEGER, SAVE :: node,nprocs,nstart,nend,nlocal
  INTEGER, SAVE :: No_of_Calls=0,natom
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
    
#ifdef PARALLEL    
    CALL MPI_ALLREDUCE(nlocal,natom,1,MPI_INTEGER4,MPI_SUM,MPI_COMM_WORLD,ierr)
    Write(*,*) natom
    n=natom
    ntot=n*3
    ALLOCATE(fp(ntot),fp_out(ntot))
#endif
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

#ifndef PARALLEL
    RETURN
#endif     
#ifdef PARALLEL
    IF(PRESENT(nstart_a)) THEN
       CALL Init(nstart_a,nend_a,nlocal_a,node_a,nprocs_a)
       RETURN
    END IF
    IF(No_of_Calls < 1) THEN
       errmsg_f='Need Initialisation to reduce three vectors.'
       CALL Add_Errors(-1,errmsg_f)
       CALL Print_Errors()
    END IF
    Do k=1,natom
       fp((k-1)*3+1)=x(k)
       fp((k-1)*3+2)=y(k)
       fp((k-1)*3+3)=z(k)
    End Do
    fp_length=natom*3
    CALL MPI_ALLREDUCE(fp, fp_out, fp_length, MPI_REAL8, MPI_SUM,&
         & MPI_COMM_WORLD, ierr) 
       
!!$------------------------------------------------------------------------
!!$--- Copy the vector in the original ones
!!$------------------------------------------------------------------------
    Do k=1,natom
       x(k)=fp_out((k-1)*3+1)
       y(k)=fp_out((k-1)*3+2)
       z(k)=fp_out((k-1)*3+3)
    End Do
#endif
    
    
!!$*----------------- END OF EXECUTABLE STATEMENTS -----------------------*
    
  END SUBROUTINE Reduce_Forces
  SUBROUTINE Reduce_Phi(x)
    REAL(8) :: x(*)
#ifdef PARALLEL
    include 'mpif.h'
    include 'mpi_size.h'
#endif
    INTEGER :: ierr
    
    INTEGER :: n,ntot,i,k,ii,fp_length


#ifdef PARALLEL
    IF(No_of_Calls < 1) THEN
       errmsg_f='Need Initialisation to reduce vector.'
       CALL Add_Errors(-1,errmsg_f)
       CALL Print_Errors()
    END IF

    Do k=1,natom
       fp(k)=x(k)
    End Do
    fp_length=natom
    CALL MPI_ALLREDUCE(fp, fp_out, fp_length, MPI_REAL8, MPI_SUM,&
         & MPI_COMM_WORLD, ierr)
       
!!$------------------------------------------------------------------------
!!$--- Copy the vector in the original ones
!!$------------------------------------------------------------------------
    Do k=1,natom
       x(k)=fp_out(k)
    End Do
#endif
  END SUBROUTINE Reduce_Phi
END MODULE REDUCE
