MODULE TYPES

!!$***********************************************************************
!!$   Time-stamp: <2006-12-15 17:03:18 marchi>                           *
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
  TYPE Unit_Char
     CHARACTER(len=max_char) :: FField=' ',Type=' ',Residue=' '
     CHARACTER(len=max_char) , DIMENSION(:,:), ALLOCATABLE :: bonds,imph&
          &,acc,acc_,don,don_,dele,ends
     CHARACTER(len=max_char), DIMENSION(:,:), ALLOCATABLE :: mass
     TYPE(list), DIMENSION(:), ALLOCATABLE :: group
  END TYPE Unit_Char
  TYPE KHASH
     CHARACTER(len=max_char) :: Type
     CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE :: keys
     INTEGER :: n,Method=0
  END TYPE KHASH
  TYPE Resid
     CHARACTER(len=max_char) :: Type,Residue,FField
     CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE :: line
  END TYPE Resid
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
  TYPE Chain_M
     INTEGER :: pt 
     REAL(8), DIMENSION(:), ALLOCATABLE :: g
  END TYPE Chain_M
  TYPE ChainR8C
     CHARACTER(len=max_char) :: lab
     REAL(8), DIMENSION(:), ALLOCATABLE :: g
  END TYPE ChainR8C
END MODULE TYPES
