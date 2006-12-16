MODULE PARAMETERS_GLOBALS

!!$***********************************************************************
!!$   Time-stamp: <2006-12-04 13:56:49 marchi>                           *
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
  TYPE(Resid), DIMENSION(:), ALLOCATABLE, SAVE  :: Topology,Paras
  CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE :: input_data
  TYPE(keys), SAVE :: in_str
  TYPE(Param), DIMENSION(2) :: Secondary_Seq
  TYPE(PATCH), DIMENSION(:), ALLOCATABLE :: patches  
  TYPE(Unit_Char), DIMENSION(:), ALLOCATABLE, SAVE, TARGET :: Res_char&
       &,Add_Char
  TYPE(Unit_Char), DIMENSION(:), ALLOCATABLE, SAVE :: App_Char
END MODULE PARAMETERS_GLOBALS
