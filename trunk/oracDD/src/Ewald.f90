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
MODULE Ewald
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Wed Feb 14 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*
  
  USE PI_Communicate
  USE FactorizeNo
  USE Potential
  IMPLICIT none
  PRIVATE
  PUBLIC Ewald__Validate
CONTAINS
  SUBROUTINE Ewald__Validate
    INTEGER, PARAMETER :: n_Max=7
    INTEGER :: n,nx,ny,nz,nprocs,npx,npy,npz

    IF(Ewald__Param % do_not_change) RETURN

    CALL PI__GetParameters(nprocs,npx,npy,npz)
    nx=Ewald__Param % nx
    ny=Ewald__Param % ny
    nz=Ewald__Param % nz
    CALL FindGrid(nx)
    CALL FindGrid(ny)
    nz=NINT(DBLE(nz)/DBLE(nprocs))*nprocs
    Ewald__Param % nx = nx
    Ewald__Param % ny = ny
    Ewald__Param % nz = nz
    WRITE(*,*) 'Ewald Grid has been changed. New grid:  X = ',nx,' Y = ',ny,' Z = ',nz
  CONTAINS
    SUBROUTINE FindGrid(nxy)
      INTEGER :: nxy
      LOGICAL :: ok,ok_i
      INTEGER, POINTER :: vect(:)=>NULL()

      ok=.FALSE.
      DO WHILE(.NOT. ok) 
         vect=>FactorizeNo_(nxy, n_Max)
         ok_i=.TRUE.
         DO n=1,SIZE(vect)
            IF(vect(n) > n_Max) ok_i=.FALSE.
         END DO
         IF(.NOT. ok_i) THEN
            nxy=nxy+1
         ELSE
            ok=.TRUE.
         END IF
      END DO
      nxy=PRODUCT(vect)
    END SUBROUTINE FindGrid
  END SUBROUTINE Ewald__Validate
END MODULE Ewald
