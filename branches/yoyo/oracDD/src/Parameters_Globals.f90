MODULE Parameters_globals

!!$***********************************************************************
!!$   Time-stamp: <2007-01-05 15:21:42 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Wed Nov 22 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program  ORAC ----*

  USE TYPES
  INTEGER, SAVE :: ktpg_read=0,kpar_read=0,kbin=0
  CHARACTER(len=max_char), SAVE :: ftpg_read,fpar_read,fbin
  CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE :: input_data

  TYPE(keys), SAVE :: in_str
  TYPE(Param), DIMENSION(2) :: Secondary_Seq
  TYPE(PATCH), DIMENSION(:), ALLOCATABLE :: patches  
END MODULE Parameters_globals
