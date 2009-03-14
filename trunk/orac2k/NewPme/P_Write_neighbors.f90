SUBROUTINE P_Write_neighbors(rshell,npp,kprint,node,nprocs)

!!$***********************************************************************
!!$   Time-stamp: <2005-01-28 19:29:59 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Tue Jan  4 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*
  
  IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*
  INCLUDE 'mpif.h'
  INCLUDE 'mpi_size.h'
  INTEGER :: npp,kprint,node,nprocs
  CHARACTER(1) :: rshell
  
!!$------------------------- LOCAL VARIABLES ----------------------------*

  INTEGER, DIMENSION (:), ALLOCATABLE ::  npps
  INTEGER :: ierr,i

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

  ALLOCATE(npps(nprocs))
  CALL MPI_ALLGATHER(npp,1,MPI_INTEGER4,npps,1,MPI_INTEGER4&
         &,MPI_COMM_WORLD,ierr)
  IF(node .EQ. 0) THEN
     WRITE(kprint,100) 
     WRITE(kprint,200) rshell
     WRITE(kprint,300) (npps(i),i=1,nprocs)
  END IF
  DEALLOCATE(npps)
  
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*
100 FORMAT(/'********************************** Neighbors *************',& 
         &'********************')
200 FORMAT(20x,' Shell = ',a1)
300 FORMAT(10x,10i8)
END SUBROUTINE P_Write_neighbors
