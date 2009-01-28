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
MODULE BoxGeometry
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Wed Jan 28 2009 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*

#include "BoxGeometry.h"
  USE Geometry
  IMPLICIT none
  PRIVATE
  PUBLIC BoxGeometry_
  Type :: face
     INTEGER :: Dir,Axis
     type(plane) :: pL
     type(point) :: a,b,c,d
     type(line) :: ab,bc,cd,da,ac
  End Type face
  Type(face), Save :: Surf(__MaxSide)
  Real(8), Save :: co(3,3),oc(3,3)
  Type(point) ::  MyOrigin, MySides
  Integer, Save :: counter=0
  Type(point), SAVE :: p001=point(0.0_8,0.0_8,1.0_8),p010=point(0.0_8&
       &,1.0_8,0.0_8),p100=point(1.0_8,0.0_8,0.0_8)
CONTAINS
  Subroutine BoxGeometry_(MyOrigina, MyLengtha, coa, oca, Dir,Axis)
    Integer :: Dir,Axis !- Axis is in the scaled frame
    Real(8), Optional :: MyOrigina(3),MyLengtha(3),coa(3,3), oca(3,3)
    Type(point) :: p1,p2,p3,p4
    Real(8) :: aux1

!!$--- Initialize the My Box, the box the CPU is working on
    If(Present(MyOrigina)) THEN
       MyOrigin % x =MyOrigina(1); MyOrigin % y =MyOrigina(2); MyOrigin % z =MyOrigina(3)
       MySides % x =MyLengtha(1); MySides % y =MyLengtha(2); MySides % z =MyLengtha(3)
       co=coa
       oc=oca
       Return
    End If
    counter=counter+1
    Surf(counter)%Dir=Dir
    Surf(counter)%Axis=Axis
    aux1=0.5D0*(1.0D0+DBLE(Dir))

    If(aux1 == 0.0_8) THEN
       select Case(Axis)
       Case(_x_) 
          p1=MyOrigin
          p2=p1+p001*MySides
          p3=p2+p010*MySides
          p4=p3-p001*MySides
       Case(_y_)
          p1=MyOrigin
          p2=p1+p100*MySides
          p3=p2+p001*MySides
          p4=p3-p100*MySides
       Case(_z_)
          p1=MyOrigin
          p2=p1+p010*MySides
          p3=p2+p100*MySides
          p4=p3-p010*MySides
       End select
    Else
       select Case(Axis)
       Case(_x_) 
          p1=MyOrigin+aux1*MySides*p100
          p2=p1+p010*MySides
          p3=p2+p001*MySides
          p4=p3-p010*MySides
       Case(_y_)
          p1=MyOrigin+aux1*MySides*p010
          p2=p1+p001*MySides
          p3=p2+p100*MySides
          p4=p3-p001*MySides
       Case(_z_)
          p1=MyOrigin+aux1*MySides*p001
          p2=p1+p100*MySides
          p3=p2+p010*MySides
          p4=p3-p100*MySides
       End select
    End If
    
    
    
  End Subroutine BoxGeometry_
END MODULE BoxGeometry
