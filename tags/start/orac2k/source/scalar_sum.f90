SUBROUTINE scalar_sum(node,nodey,nodez,nay,naz,Q,ewaldcof,volume&
     &,recip,bsp_mod1,bsp_mod2,bsp_mod3,nfft1,nfft2,nfft3,nfftdim1&
     &,nfftdim2,nfftdim3,eer,vir,rkcut)

!!$***********************************************************************
!!$   Time-stamp: <2004-12-09 19:50:13 marchi>                           *
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

  IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

  INTEGER :: nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,node,nodey&
       &,nodez,nay,naz
  REAL(8) ::  Q(nfftdim1,nfftdim2,nfftdim3)
  REAL(8) ::  bsp_mod1(nfft1),bsp_mod2(nfft2),bsp_mod3(nfft3)&
       &,ewaldcof,volume,rkcut
  REAL(8) :: eer,vir(3,3)
  REAL(8) :: recip(3,3)


!!$----------------------- VARIABLES IN COMMON --------------------------*

!!$------------------------- LOCAL VARIABLES ----------------------------*


!!$----------------------- EXECUTABLE STATEMENTS ------------------------*





!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE scalar_sum
