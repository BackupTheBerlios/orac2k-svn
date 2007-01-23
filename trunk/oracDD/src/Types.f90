!!$/---------------------------------------------------------------------\
!!$   Copyright  © 2006-2007 Massimo Marchi <Massimo.Marchi at cea.fr>   |
!!$                                                                      |
!!$    This software is a computer program named oracDD whose            |
!!$    purpose is to simulate and model complex molecular systems.       |
!!$    The code is written in fortran 95 compliant with Technical        |
!!$    Report TR 15581, and uses MPI-1 routines for parallel             |
!!$    coding.                                                           |
!!$                                                                      |
!!$    This software is governed by the CeCILL license under             |
!!$    French law and abiding by the rules of distribution of            |
!!$    free software.  You can  use, modify and/ or redistribute         |
!!$    the software under the terms of the CeCILL icense as              |
!!$    circulated by CEA, CNRS and INRIA at the following URL            |
!!$    "http://www.cecill.info".                                         |
!!$                                                                      |
!!$    As a counterpart to the access to the source code and rights      |
!!$    to copy, modify and redistribute granted by the license,          |
!!$    users are provided only with a limited warranty and the           |
!!$    software's author, the holder of the economic rights, and         |
!!$    the successive licensors have only limited liability.             |
!!$                                                                      |
!!$    The fact that you are presently reading this means that you       |
!!$    have had knowledge of the CeCILL license and that you accept      |
!!$    its terms.                                                        |
!!$                                                                      |
!!$    You should have received a copy of the CeCILL license along       |
!!$    with this program; if not, you can collect copies on the URL's    |
!!$    "http://www.cecill.info/licences/Licence_CeCILL_V2-en.html"       |
!!$    "http://www.cecill.info/licences/Licence_CeCILL_V2-fr.html"       |
!!$                                                                      |
!!$----------------------------------------------------------------------/
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
