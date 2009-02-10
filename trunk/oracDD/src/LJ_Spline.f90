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
MODULE LJ_Spline
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Mon Jan 26 2009 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*

#include "LJ_Spline.h"                       

  USE Units
  IMPLICIT none
  PRIVATE
  PUBLIC  Spline_,s,taus,LJ_Switch,dx,xbeg
  REAL(8), SAVE :: dx=0.002_8,xbeg=0.01_8,xend,tol=1.0_8
  INTEGER, SAVE :: np
  REAL(8), ALLOCATABLE, SAVE :: s(:,:),taus(:)
  LOGICAL, SAVE :: LJ_Switch=.TRUE.
CONTAINS
  SUBROUTINE Spline_(rcut)
    REAL(8) :: rcut
    np=INT((rcut+Tol-xbeg)/dx)+1
    xend=DBLE(np-1)*dx+xbeg
    ALLOCATE(s(4,np),taus(np))
    CALL Fill_Tables
  CONTAINS
    SUBROUTINE Fill_Tables
      INTEGER :: n,ibeg,iend
      REAL(8) :: x 

      ibeg=1; iend=np

      __dlj(xbeg,s(2,1))
      __dlj(xend,s(2,np))

      DO n=1,np
         x=dx*DBLE(n-1)+xbeg
         taus(n)=x
         __lj(x,s(1,n))
      END DO
      CALL Cubspl(taus,s,np,1,1)
    END SUBROUTINE Fill_Tables
  END SUBROUTINE Spline_
  SUBROUTINE Test_it
    REAL(8) :: bin,offset,rend
    INTEGER :: length,n,i_c
    REAL(8) :: x,y,h,New_Lj,New_Dlj,yy
    bin=0.011
    offset=0.1_8
    rend=10.0_8
    length=INT((rend-offset)/bin)+1
    DO n=1,length
       x=DBLE(n-1)*bin+offset

       __lj(x,y)
       __dlj(x,yy)
       i_c=INT((x-xbeg)/dx)+1
       h=x-taus(i_c)

       New_Lj=_Lj(i_c,h)

       New_Dlj=_DLj(i_c,h)

       WRITE(*,'(2i6,5e15.7)') n,i_c,x,y,New_Lj,yy,New_Dlj
    END DO
  END SUBROUTINE Test_it
END MODULE LJ_Spline
