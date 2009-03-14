SUBROUTINE Setup_fake_decompa(nstart_h,nend_h,nlocal_h,nstart_ah&
     &,nend_ah,nlocal_ah,nstart_uh,nend_uh,nlocal_uh,ntap,ngrp,nbun)

!!$***********************************************************************
!!$   Time-stamp: <2005-01-29 11:34:41 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Sat Jan 29 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*

  IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*
  
  INTEGER nstart_h,nend_h,nlocal_h,nstart_ah,nend_ah,nlocal_ah&
       &,nstart_uh,nend_uh,nlocal_uh,ntap,ngrp,nbun

!!$----------------------- VARIABLES IN COMMON --------------------------*

!!$------------------------- LOCAL VARIABLES ----------------------------*


!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

  nstart_h=1
  nend_h=ngrp
  nlocal_h=ngrp
  nstart_ah=1
  nend_ah=ntap
  nlocal_ah=ntap
  nstart_uh=1
  nend_uh=nbun
  nlocal_uh=nbun

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE Setup_fake_decompa
