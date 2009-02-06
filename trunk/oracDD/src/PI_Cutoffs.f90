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
MODULE PI_Cutoffs
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Nov  7 2008 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*


  USE Cell
  USE Geometry, ONLY: Plane, Point, Equation_Plane, PointPlane_Dist
  USE Forces, ONLY: Force, Radii
  USE PI_
  IMPLICIT none
  PRIVATE
  PUBLIC Thickness,Thick_Intra,rcut,ddx,ddy,ddz
  REAL(8), SAVE :: dx,dy,dz,ddx,ddy,ddz
  TYPE :: Cuts
     REAL(8) :: r(3)=0.0_8
  END type Cuts
  TYPE(Cuts), SAVE :: rcuts0(10)
  REAL(8), SAVE :: Thick(3),rcut(3)
CONTAINS
  SUBROUTINE Thickness(i_pa)
    INTEGER :: i_pa
    INTEGER :: i_p
    REAL(8) :: x1,x2,y1,y2,z1,z2
    type(point) :: p1,p2,p3,v1,v2,v3,r1,r2,r3,r0
    type(plane) :: qq
    
    i_p=i_pa-2
    IF(rcuts0(i_pa) % r(1) /= 0.0_8) THEN
       rcut(:)=rcuts0(i_pa) % r(:)
       RETURN
    END IF

    ddx=2.d0/PI_npx
    ddy=2.d0/PI_npy
    ddz=2.d0/PI_npz

    r0%x=0.0D0; r0%y=0.0D0; r0%z=0.0D0

    v1%x=1.0D0
    v1%y=0.0D0
    v1%z=0.0D0
    r1=Convert(v1)
    v2%x=0.0D0
    v2%y=1.0D0
    v2%z=0.0D0
    r2=Convert(v2)
    v3%x=0.0D0
    v3%y=0.0D0
    v3%z=1.0D0
    r3=Convert(v3)

!!$    
!!$---- Thickness along x
!!$

    qq=Equation_Plane(r0,r2,r3)

    Thick(1)=1.0D0/ABS(PointPlane_Dist(qq,r1))

!!$    
!!$---- Thickness along y
!!$
    qq=Equation_Plane(r0,r1,r3)
    Thick(2)=1.0D0/ABS(PointPlane_Dist(qq,r2))

!!$    
!!$---- Thickness along z
!!$
    qq=Equation_Plane(r0,r1,r2)
    Thick(3)=1.0D0/ABS(PointPlane_Dist(qq,r3))
!!$    
!!$---- Thickness along x
!!$
    rcuts0(i_pa) % r (:)=(Radii(i_p) % out+Radii(i_p)% update)*Thick(:)
    rcut(:)=rcuts0(i_pa) % r (:)
  END SUBROUTINE Thickness
  SUBROUTINE Thick_Intra(i_p,MyCutoff)
    REAL(8) :: MyCutoff
    INTEGER :: i_p
    REAL(8) :: x1,x2,y1,y2,z1,z2
    type(point) :: p1,p2,p3,v1,v2,v3,r1,r2,r3,r0
    type(plane) :: qq
    
    IF(rcuts0(i_p) % r(1) /= 0.0_8) THEN
       rcut(:)=0.0_8
    END IF

    ddx=2.d0/PI_npx
    ddy=2.d0/PI_npy
    ddz=2.d0/PI_npz

    r0%x=0.0D0; r0%y=0.0D0; r0%z=0.0D0

    v1%x=1.0D0
    v1%y=0.0D0
    v1%z=0.0D0
    r1=Convert(v1)
    v2%x=0.0D0
    v2%y=1.0D0
    v2%z=0.0D0
    r2=Convert(v2)
    v3%x=0.0D0
    v3%y=0.0D0
    v3%z=1.0D0
    r3=Convert(v3)

!!$    
!!$---- Thickness along x
!!$
    qq=Equation_Plane(r0,r2,r3)
    Thick(1)=1.0D0/ABS(PointPlane_Dist(qq,r1))
!!$    
!!$---- Thickness along y
!!$
    qq=Equation_Plane(r0,r1,r3)
    Thick(2)=1.0D0/ABS(PointPlane_Dist(qq,r2))
!!$    
!!$---- Thickness along z
!!$
    qq=Equation_Plane(r0,r1,r2)
    Thick(3)=1.0D0/ABS(PointPlane_Dist(qq,r3))
    
!!$    
!!$---- Thickness along x
!!$
    rcuts0(i_p) % r (:)=MyCutoff*Thick(:)
    rcut(:)=rcuts0(i_p) % r (:)
  END SUBROUTINE Thick_Intra
  FUNCTION Convert(v1) RESULT(out)
    type(point) :: v1,out
    type(point) :: r1
    r1%x=co(1,1)*v1%x+co(1,2)*v1%y+co(1,3)*v1%z
    r1%y=co(2,1)*v1%x+co(2,2)*v1%y+co(2,3)*v1%z
    r1%z=co(3,1)*v1%x+co(3,2)*v1%y+co(3,3)*v1%z
    out=r1
  END FUNCTION Convert

END MODULE PI_Cutoffs
