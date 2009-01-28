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
#include "Geometry.h"
  USE PI_
  IMPLICIT none
  PRIVATE
  PUBLIC  Point, Plane, Equation_Plane,PointPlane_Dist,Line,operator(&
       &+),operator(-),operator(*),operator(==) 

  TYPE :: Point
     REAL(8) :: x,y,z
  END type Point
  type :: Plane
     type(point) :: p
     REAL(8) :: d
  end type Plane
  type :: line
     type(point) :: a,b
  END type line
  interface operator(+)
     module procedure add_points
  end interface
  interface operator(-)
     module procedure subtract_points
  end interface
  interface operator(*)
     module procedure multiply_points
     module procedure multiplyscalar_points
  end interface
  interface operator(==)
     module procedure equal_points
  end interface
CONTAINS
  FUNCTION Equation_Plane(v1,v2,v3) RESULT(out)
    type(plane) :: out
    type(point) :: v1,v2,v3
    TYPE(Point) :: p1,p2,p3,p4
    type(plane) :: q
    
    REAL(8) :: Norm,dist

    p1=v1 ; p2=v2; p3=v3

    __subtract_v(p2, p1, p2)
    __subtract_v(p3, p1, p3)
    __crossprod(p2,p3,p4)
    __normalize_v(p4,Norm)
    __dotprod(p4,p1,Dist)
    
    q % p = p4
    q % d = Dist
    out=q
  END FUNCTION Equation_Plane
  FUNCTION Equation_Line(v1,v2) RESULT(out)
    type(point) :: v1,v2
    type(point) :: p1,p2
    type(line) :: out
    
    
    p1=v1 ; p2=v2; __subtract_v(p2, p1, p2)
    out % a=p1; out % b=p2
    
  END FUNCTION Equation_Line
  FUNCTION  Foot(v1,l1) RESULT(out)
    type(line) :: l1
    type(point) :: v1,out
    type(point) :: p1,p2
    real(8) :: t0,dist1,dist2
    
    p1=v1

    p2=p1-l1%a

    __dotprod(p2,l1%b,Dist1)
    __dotprod(l1%b,l1%b,Dist2)

    t0=Dist1/Dist2

    out=l1%a+t0*l1%b

  END FUNCTION Foot
  FUNCTION PointPlane_Dist(pl,p) RESULT(out)
    REAL(8) :: out
    REAL(8) :: Dist
    Type(Point) :: p
    Type(Plane) :: pL

    __dotprod(pL % p, p, Dist)
    out=Dist-pL % d
  END FUNCTION PointPlane_Dist
  function add_points(a,b) result(out)
    type(point), intent(in) :: a,b
    type(point) :: out
    __add_v(a,b,out)
  end function add_points
  function subtract_points(a,b) result(out)
    type(point), intent(in) :: a,b
    type(point) :: out
    __subtract_v(a,b,out)
  end function subtract_points
  function multiply_points(a,b) result(out)
    type(point), intent(in) :: a,b
    type(point) :: out
    out%x=a%x*b%x
    out%y=a%y*b%y
    out%z=a%z*b%z
  end function multiply_points
  function multiplyscalar_points(a,b) result(out)
    type(point), intent(in) :: b
    real(8), intent(in) :: a
    type(point) :: out
    out%x=a*b%x
    out%y=a*b%y
    out%z=a*b%z
  end function multiplyscalar_points
  function equal_points(a,b) result(out)
    type(point), intent(in) :: a,b
    LOGICAL :: out
    real(8), save :: tol_p=1.0D-8

    out=(Abs(a%x-b%x) < tol_p) .And. (Abs(a%y-b%y) < tol_p) .And. (Abs(a%z-b%z) < tol_p) 
  end function equal_points
END MODULE Geometry
