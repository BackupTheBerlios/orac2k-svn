!!$/---------------------------------------------------------------------\
!!$                                                                      |
!!$  Copyright (C) 2006-2007 Massimo Marchi <Massimo.Marchi@cea.fr>      |
!!$                                                                      |
!!$      This program is free software;  you  can  redistribute  it      |
!!$      and/or modify it under the terms of the GNU General Public      |
!!$      License version 2 as published  by  the  Free  Software         |
!!$      Foundation;                                                     |
!!$                                                                      |
!!$      This program is distributed in the hope that  it  will  be      |
!!$      useful, but WITHOUT ANY WARRANTY; without even the implied      |
!!$      warranty of MERCHANTABILITY or FITNESS  FOR  A  PARTICULAR      |
!!$      PURPOSE.   See  the  GNU  General  Public License for more      |
!!$      details.                                                        |
!!$                                                                      |
!!$      You should have received a copy of the GNU General  Public      |
!!$      License along with this program; if not, write to the Free      |
!!$      Software Foundation, Inc., 59  Temple  Place,  Suite  330,      |
!!$      Boston, MA  02111-1307  USA                                     |
!!$                                                                      |
!!$\---------------------------------------------------------------------/
MODULE CONSTANTS
!!$***********************************************************************
!!$   Time-stamp: <2007-01-12 20:28:23 marchi>                           *
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
       & max_char_long = 12000, max_data=12000,max_char=120,max_atm=10
  CHARACTER(len=1), DIMENSION(2), PARAMETER :: Comms=(/'!','#'/)
  CHARACTER(len=max_char), DIMENSION(9), PARAMETER :: Used=(/'ATOM  ','H&
       &ETATM','CONECT','SSBOND','CRYST1','SEQRES','HELIX ','SHEET ','CI&
       &SPEP'/)
  REAL(8), PARAMETER :: Huge=1.0D10,Tiny=1.0D-10
END MODULE CONSTANTS
