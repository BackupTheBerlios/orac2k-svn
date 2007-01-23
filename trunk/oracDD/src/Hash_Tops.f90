!!$/---------------------------------------------------------------------\
!!$   Copyright  � 2006-2007 Massimo Marchi <Massimo.Marchi at cea.fr>   |
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
MODULE Hash_Tops

!!$***********************************************************************
!!$   Time-stamp: <2007-01-10 15:29:32 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Jan  5 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*

  USE Constants
  USE Errors, ONLY: errmsg_f, Add_Errors=>Add, Print_Errors
  IMPLICIT none
  PRIVATE
  PUBLIC Hash_Tops_, Hash_Tops__Extract, Hash_Tops__Type
  TYPE :: Hash_Tops__Type
     CHARACTER(len=max_char) :: Type
     CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE :: keys
     INTEGER :: n,Method=0
  END TYPE Hash_Tops__Type
  TYPE(Hash_Tops__Type), DIMENSION(8), SAVE :: Hash

CONTAINS
  SUBROUTINE Hash_Tops_
    Hash(1)%Type='bond'
    Hash(2)%Type='imph'
    Hash(3)%Type='acc'
    Hash(4)%Type='acc_'
    Hash(5)%Type='don'
    Hash(6)%Type='don_'
    Hash(7)%Type='dele'
    Hash(8)%Type='term'
    ALLOCATE(Hash(1)%keys(3))
    Hash(1)%keys=(/'bond','doub','trip'/)
    Hash(1)%n=2
    ALLOCATE(Hash(2)%keys(2))
    Hash(2)%keys=(/'impr','imph'/)
    Hash(2)%n=4
    ALLOCATE(Hash(3)%keys(1))
    Hash(3)%n=2
    Hash(3)%keys=(/'acce'/)
    ALLOCATE(Hash(4)%keys(1))
    Hash(4)%n=1
    Hash(4)%keys=(/'acce'/)
    ALLOCATE(Hash(5)%keys(1))
    Hash(5)%n=2
    Hash(5)%keys=(/'dono'/)
    ALLOCATE(Hash(6)%keys(1))
    Hash(6)%n=1
    Hash(6)%keys=(/'dono'/)
    ALLOCATE(Hash(7)%keys(1))
    Hash(7)%Method=-1
    Hash(7)%n=1
    Hash(7)%keys=(/'dele'/)
    ALLOCATE(Hash(8)%keys(1))
    Hash(8)%n=2
    Hash(8)%keys=(/'term'/)
  END SUBROUTINE Hash_Tops_
  FUNCTION Hash_Tops__Extract(types) RESULT(out)
    TYPE(Hash_Tops__Type) :: out
    CHARACTER(len=*) :: types
    INTEGER :: n,IHash

    IHash=-1
    DO n=1,SIZE(Hash)
       IF(TRIM(Hash(n)%type) == TRIM(types)) THEN
          IHash=n
       END IF
    END DO
    IF(IHash == -1) THEN
       errmsg_f='ORAC internal: Type '&
            &//TRIM(types)//' not found in hash Tables'
       CALL Add_Errors(-1,errmsg_f)
       CALL Print_Errors()
    ELSE
       out=Hash(IHash)
    END IF
  END FUNCTION Hash_Tops__Extract

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE Hash_Tops
