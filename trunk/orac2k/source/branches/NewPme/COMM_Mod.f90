MODULE COMM_Mod

!!$***********************************************************************
!!$   Time-stamp: <2006-04-07 16:21:36 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Apr  7 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*

  INTEGER, SAVE :: node,nprocs,ncube

!!$======================== DECLARATIONS ================================*

CONTAINS
  SUBROUTINE Init(nodea,nprocsa,ncubea)
    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

    INTEGER :: nodea,nprocsa,ncubea

!!$----------------------- VARIABLES IN COMMON --------------------------*

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    node=nodea
    nprocsa=nprocs
    ncube=ncubea

  END SUBROUTINE Init

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE COMM_Mod
