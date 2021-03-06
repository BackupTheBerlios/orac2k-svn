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
MODULE Erfc_Spline
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

#include "Erfc_Spline.h"                       

  USE Units
  IMPLICIT none
  PRIVATE
  PUBLIC  Erfc_,c,tau,Erfc_Switch,dx,Xbeg
  REAL(8), SAVE :: dx=0.01_8,xbeg=1.0_8,xend,tol=2.0_8
  INTEGER, SAVE :: np
  REAL(8), ALLOCATABLE, SAVE :: c(:,:),tau(:)
  LOGICAL, SAVE :: Erfc_Switch=.TRUE.
CONTAINS
  SUBROUTINE Erfc_(rcut,alphal)
    REAL(8) :: rcut,alphal
    np=INT((rcut+Tol-xbeg)/dx)+1
    xend=DBLE(np-1)*dx+xbeg
    ALLOCATE(c(4,np),tau(np))
    CALL Fill_Tables
  CONTAINS
    SUBROUTINE Fill_Tables
      INTEGER :: n,ibeg,iend
      REAL(8) :: x 

      ibeg=1; iend=np
      c(2,1) =(-2*alphal)/(EXP(alphal**2*xbeg**2)*Sqrt(Pi)*xbeg) - MyErfc(alphal*xbeg)/xbeg**2
      c(2,np)=(-2*alphal)/(EXP(alphal**2*xend**2)*Sqrt(Pi)*xend) - MyErfc(alphal*xend)/xend**2
      DO n=1,np
         x=dx*DBLE(n-1)+xbeg
         tau(n)=x
         c(1,n)=MyErfc(alphal*x)/x
      END DO
      CALL Cubspl(tau,c,np,1,1)
    END SUBROUTINE Fill_Tables
  END SUBROUTINE Erfc_
  SUBROUTINE Test_it(alphal)
    REAL(8) :: alphal
    REAL(8) :: bin,offset,rend
    INTEGER :: length,n,i_c
    REAL(8) :: x,y,h,New_ERfc,New_Derfc,yy
    bin=0.011
    offset=2.0_8
    rend=8.0_8
    length=INT((rend-offset)/bin)+1
    WRITE(*,*) alphal,MyErfc(alphal*2.0_8)/2.0_8
    DO n=1,length
       x=DBLE(n-1)*bin+offset
       y=ERFC(alphal*x)/x
       
       yy=(-2*alphal)/(EXP(alphal**2*x**2)*Sqrt(Pi)*x) - MyErfc(alphal*x)/x**2

       i_c=INT((x-xbeg)/dx)+1
       h=x-tau(i_c)


       New_Erfc=_Erfc(i_c,h)

       New_Derfc=_DErfc(i_c,h)

       WRITE(*,'(2i6,5e15.7)') n,i_c,x,y,yy,New_Erfc,New_Derfc
    END DO
  END SUBROUTINE Test_it
END MODULE Erfc_Spline
