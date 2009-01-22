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
MODULE CellAtoms
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Jan 16 2009 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*

#include "CellAtoms.h"

  IMPLICIT none
  PRIVATE
  PUBLIC CellAtoms_at_,CellAtoms_gr_,CellAtoms_kn_,CellAtoms_at__&
       &,CellAtoms_gr__,CellAtoms_kn__ 

  REAL(8), ALLOCATABLE, SAVE :: at(:,:),gr(:,:)
  INTEGER, ALLOCATABLE, SAVE :: Idx(:)
  INTEGER, SAVE :: Idx_Max=0
  INTEGER, SAVE :: FCalls=0
CONTAINS
  FUNCTION CellAtoms_(n,m
  FUNCTION CellAtoms_at_(x,y,z,Ind) RESULT(out)
    REAL(8) :: x(:),y(:),z(:)
    INTEGER, OPTIONAL :: Ind(:)
    REAL(8), POINTER :: out(:,:)
    INTEGER :: n
    IF(PRESENT(Ind)) THEN
       n=SIZE(Ind)
       WRITE(*,*) ' at n = ',n
       ALLOCATE(at(3,n))
       at(1,:)=x(Ind(:))
       at(2,:)=y(Ind(:))
       at(3,:)=z(Ind(:))
    END IF
    out=>at
  END FUNCTION CellAtoms_at_
  FUNCTION CellAtoms_gr_(x,y,z,Ind) RESULT(out)
    REAL(8), OPTIONAL :: x(:),y(:),z(:)
    INTEGER, OPTIONAL :: Ind(:)
    REAL(8), POINTER :: out(:,:)
    INTEGER :: n

    IF(PRESENT(Ind)) THEN
       n=SIZE(Ind)
       WRITE(*,*) ' gr n = ',n
       ALLOCATE(gr(3,n))
       gr(1,:)=x(Ind(:))
       gr(2,:)=y(Ind(:))
       gr(3,:)=z(Ind(:))
    END IF
    out=>gr
  END FUNCTION CellAtoms_gr_
  FUNCTION CellAtoms_kn_(k,Ind) RESULT(out)
    INTEGER, OPTIONAL :: Ind(:),k(:)
    INTEGER, POINTER :: out(:)
    INTEGER :: n

    IF(PRESENT(Ind)) THEN
       n=SIZE(Ind)
       ALLOCATE(knwn(n))
       knwn(:)=k(Ind(:))
    END IF
    out=>knwn
  END FUNCTION CellAtoms_kn_
  FUNCTION CellAtoms_at__(x,y,z,Ind) RESULT(out)
    INTEGER :: Ind(:)
    REAL(8), POINTER :: out(:,:)
    REAL(8) :: x(:),y(:),z(:)
    
    x(Ind(:))=at(1,:)
    y(Ind(:))=at(2,:)
    z(Ind(:))=at(3,:)
    DEALLOCATE(at)
    out=>NULL()
  END FUNCTION CellAtoms_at__
  FUNCTION CellAtoms_gr__(x,y,z,Ind) RESULT(out)
    INTEGER :: Ind(:)
    REAL(8), POINTER :: out(:,:)
    REAL(8) :: x(:),y(:),z(:)

    x(Ind(:))=gr(1,:)
    y(Ind(:))=gr(2,:)
    z(Ind(:))=gr(3,:)
    DEALLOCATE(gr)
    out=>NULL()
  END FUNCTION CellAtoms_gr__
  FUNCTION CellAtoms_kn__(k,Ind) RESULT(out)
    INTEGER :: Ind(:),k(:)
    INTEGER, POINTER :: out(:)

    k(Ind(:))=knwn(:)
    out=>NULL()
    DEALLOCATE(knwn)
  END FUNCTION CellAtoms_kn__
END MODULE CellAtoms
