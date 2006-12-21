MODULE CONSTANTS
!!$***********************************************************************
!!$   Time-stamp: <2006-12-21 12:29:39 marchi>                           *
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
  CHARACTER(len=max_char), DIMENSION(9), PARAMETER :: Used=(/'ATOM  ','H&
       &ETATM','CONECT','SSBOND','CRYST1','SEQRES','HELIX ','SHEET ','CI&
       &SPEP'/)
  REAL(8), PARAMETER :: Huge=1.0D10,Tiny=1.0D-10
END MODULE CONSTANTS
