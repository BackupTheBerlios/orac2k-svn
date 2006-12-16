MODULE CONSTANTS
!!$***********************************************************************
!!$   Time-stamp: <2006-11-24 12:43:24 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Nov 24 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program ORAC ----*
  INTEGER, PARAMETER :: max_pars=200,max_char_tree = 80,&
       & max_char_long = 12000, max_data=12000,max_char=120
  CHARACTER(len=1), DIMENSION(2), PARAMETER :: Comms=(/'!','#'/)

END MODULE CONSTANTS
