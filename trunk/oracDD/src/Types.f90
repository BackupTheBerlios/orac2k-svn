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
MODULE TYPES

!!$***********************************************************************
!!$   Time-stamp: <2007-01-10 17:11:34 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Nov 16 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program ORAC ----*

  USE CONSTANTS
  TYPE List
     CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE :: g
  END TYPE List
  TYPE KHASH
     CHARACTER(len=max_char) :: Type
     CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE :: keys
     INTEGER :: n,Method=0
  END TYPE KHASH
  TYPE Param
     CHARACTER(len=max_char) :: Type
     CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE :: line
  END TYPE Param
  TYPE Docs
     CHARACTER(len=max_char) :: l
     CHARACTER(len=120), DIMENSION(:), POINTER :: C
  END TYPE Docs
  TYPE Commands_D
     CHARACTER(len=max_char) :: l
     CHARACTER(len=120), DIMENSION(:), POINTER :: C
     TYPE(Docs), DIMENSION(:), POINTER :: keys
  END TYPE Commands_D
  TYPE Environment_D
     CHARACTER(len=max_char) :: l
     CHARACTER(len=120), DIMENSION(:), POINTER :: C
     TYPE(commands_D), DIMENSION(:), POINTER :: Comm
  END TYPE Environment_D

  TYPE link_tble
     CHARACTER(len=max_char) :: key
     TYPE(link_tble), DIMENSION(:), POINTER :: n_key
  END TYPE link_tble

  TYPE keys
     CHARACTER(len=max_char) :: label
     CHARACTER(len=max_char), DIMENSION(:), POINTER :: line
     TYPE(keys), DIMENSION(:), POINTER :: next
  END TYPE keys
  TYPE PATCH
     CHARACTER(8) :: Type
     CHARACTER(8) :: New_Res,pres,res
     CHARACTER(8) :: Res_l(2)
     INTEGER :: one,two
  END TYPE PATCH
  TYPE Chain
     INTEGER, DIMENSION(:), ALLOCATABLE :: g
  END TYPE Chain
  TYPE ChainR8
     REAL(8), DIMENSION(:), ALLOCATABLE :: g
  END TYPE ChainR8
  TYPE ChainR8C
     CHARACTER(len=max_char) :: lab
     REAL(8), DIMENSION(:), ALLOCATABLE :: g
  END TYPE ChainR8C
END MODULE TYPES
