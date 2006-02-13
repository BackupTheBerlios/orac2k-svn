SUBROUTINE Pme_init(node_first,nprocs,nfft1,nfft2,nfft3,nfft3_start&
     &,nfft3_local,nfft2_start,nfft2_local,iret,errmsg)

!!$***********************************************************************
!!$   Time-stamp: <2006-02-13 15:34:14 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Dec  9 2004 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*


  USE rfft3d
  IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

  INTEGER ::  iret,nfft1,nfft2,nfft3,nfft3_start,nfft3_local&
       &,nfft2_start,nfft2_local,node_first,nprocs

  CHARACTER(80) ::  errmsg

!!$----------------------- VARIABLES IN COMMON --------------------------*

#ifdef PARALLEL 
  INCLUDE 'mpif.h'
#endif

!!$------------------------- LOCAL VARIABLES ----------------------------*

  REAL(8) ::  scale,Q,work
  REAL(8), DIMENSION (:,:,:), ALLOCATABLE :: dummy
    
  INTEGER :: i,k1,k2,k3

  INTEGER :: ierr,node,comm1d
  INTEGER :: dims(3),ndim
  LOGICAL :: periods(3),reorder

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*


  
  IF(node_first == 0) WRITE(*,100) 

  comm1d=MPI_COMM_WORLD

  CALL do_rfft3d(0,dummy,nfft1,nfft2,nfft3,nfft3_start,nfft3_local&
       &,nfft2_start,nfft2_local,k1,k2,k3,comm1d)

100 FORMAT(/22x,'Finding optimal parameters for FFTWs.'/&
     &     22x,'     This will take a while...'/ /) 

END SUBROUTINE Pme_init
