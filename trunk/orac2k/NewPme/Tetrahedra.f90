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
MODULE Tetrahedra
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Feb 27 2009 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program orac2k ----*

  Implicit None
  Private
  Public Tetrahedra_
Contains
  Subroutine Tetrahedra_(x0,y0,z0,x1,y1,z1,x2,y2,z2,v)
    Real(8) :: x0,y0,z0,x1,y1,z1,x2,y2,z2,v(3,4)

    Real(8) :: v1(3),v2(3),v3(3),v4(3),ax(3),pax(3),xn(2),yn(2),zn(2)&
         &,Norm,Length,xh(3),yh(3),zh(3)


    Integer :: m,blist(2),iret,ia,nbl,nhb
    Character(80) :: errmsg
    Character(7) :: type

    ia=1; blist(1)=2; blist(2)=3;
    xh(1)=x0; yh(1)=y0; zh(1)=z0; 
    xh(2)=x1; yh(2)=y1; zh(2)=z1; 
    xh(3)=x2; yh(3)=y2; zh(3)=z2; 
    type='c'
    nbl=2; nhb=2
    Call addh(ia,nbl,nhb,blist,xh,yh,zh,type,xn,yn,zn,iret,errmsg)

    v(:,1)=(/xh(2),yh(2),zh(2)/)-(/xh(1),yh(1),zh(1)/)
    v(:,2)=(/xh(3),yh(3),zh(3)/)-(/xh(1),yh(1),zh(1)/)
    v(:,3)=(/xn(1),yn(1),zn(1)/)-(/xh(1),yh(1),zh(1)/)
    v(:,4)=(/xn(2),yn(2),zn(2)/)-(/xh(1),yh(1),zh(1)/)
    Do m=1,4
       Norm=Sqrt(Sum(v(:,m)**2))
       v(:,m)=v(:,m)/Norm
    End Do
  End Subroutine Tetrahedra_

End Module Tetrahedra
