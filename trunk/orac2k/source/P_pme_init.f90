SUBROUTINE P_pme_init(node,nprocs,nodex,nodey,nodez,npy,npz,ictxt&
     &,descQ,fftable,nfft1,nfft2,nfft3,nfft3_start,nfft3_local&
     &,nfft2_start,nfft2_local,iret,errmsg)

!!$***********************************************************************
!!$   Time-stamp: <2005-01-10 13:00:47 marchi>                           *
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

  INTEGER ::  iret,node,nprocs,nodex,nodey,nodez,npy,npz,ictxt,nfft1,&
       &nfft2,nfft3,nfft3_start,nfft3_local,descQ(*),nfft2_start,nfft2_local
  REAL(8) :: fftable(*)
  CHARACTER(80) ::  errmsg

!!$----------------------- VARIABLES IN COMMON --------------------------*

!!$------------------------- LOCAL VARIABLES ----------------------------*

  REAL(8) ::  scale,Q,work
  REAL(8), DIMENSION (:,:,:), ALLOCATABLE :: dummy
    
  INTEGER :: i,k1,k2,k3
  INTEGER ::info,iceil,nax,nay,descQa(12),nfft1a,nfft2a,nfft3a&
       &,npya,npza,ictxta,isys(4)
  INTEGER, SAVE :: zero=0

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*


  IF(node ==0) WRITE(*,100)
  CALL do_rfft3d(0,dummy,nfft1,nfft2,nfft3,nfft3_start,nfft3_local&
       &,nfft2_start,nfft2_local,k1,k2,k3)

100 FORMAT(/22x,'Finding optimal parameters for FFTWs.'/&
     &     22x,'     This will take a while...'/ /) 

END SUBROUTINE P_pme_init
