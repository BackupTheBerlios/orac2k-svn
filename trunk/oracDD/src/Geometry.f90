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
MODULE Geometry
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Oct 10 2008 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*

  USE PI_
  IMPLICIT none
  PRIVATE
  PUBLIC Equation_Plane
CONTAINS
  FUNCTION Equation_Plane(v1,v2,v3) RESULT(out)
    REAL(8), DIMENSION(:) :: v1,v2,v3
    REAL(8) :: out(4)
    REAL(8) :: q(4)
    
    REAL(8) :: x1,y1,z1,x2,y2,z2,x3,y3,z3,Norm

    x1=v1(1) ; y1=v1(2) ; z1=v1(3)
    x2=v2(1) ; y2=v2(2) ; z2=v2(3)
    x3=v3(1) ; y3=v3(2) ; z3=v3(3)

    q(1)=y1*(z2 - z3) + y2*(z3 - z1) + y3*(z1 - z2) 
    q(2)=z1*(x2 - x3) + z2*(x3 - x1) + z3*(x1 - x2)
    q(3)=x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2)
    q(4)=(x1*(y2*z3 - y3*z2) + x2*(y3*z1 - y1*z3) + x3*(y1*z2 - y2*z1))
    Norm=SQRT(q(1)**2+q(2)**2+q(3)**2)
    q=q/Norm
    out=q
  END FUNCTION Equation_Plane
END MODULE Geometry
