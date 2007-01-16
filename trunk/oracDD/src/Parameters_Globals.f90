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
  CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE, SAVE :: input_data

  TYPE(keys), SAVE :: in_str
  TYPE(Param), DIMENSION(2), SAVE :: Secondary_Seq
  TYPE(PATCH), DIMENSION(:), ALLOCATABLE, SAVE :: patches  
END MODULE Parameters_globals
