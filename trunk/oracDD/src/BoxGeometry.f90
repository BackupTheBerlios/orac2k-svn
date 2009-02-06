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
#include "Geometry.h"
  USE Units
  USE PI_
  USE Geometry
  USE PI_Cutoffs, ONLY: ddx,ddy,ddz
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f, errmsg_w
  USE PI_Collectives, ONLY: PI_ErrSignal_
!!$  Use NeighCells, Only: NeiMaps,ncx,ncy,ncz
  IMPLICIT none
  PRIVATE
  PUBLIC BoxGeometry_, Distance_
  Type :: face
     Integer :: Dir,Axis
     type(plane) :: pL
     type(point) :: vab,vad
     type(point) :: a,b,c,d
     type(point) :: ab,bc,cd,da,db
     type(line)  :: lab,lbc,lcd,lda
     Real(8)     :: m_ab,m_bc,m_cd,m_da
  End Type face
  Type :: RecFace
     Integer :: Dir,Axis
     Real(8) :: d
  End type RecFace
  Real(8), Save :: co(3,3),oc(3,3)
  Type(point), Save ::  MyOrigin, MySides
  Integer, Save :: counter=0
  Type(point), SAVE :: p001=point(0.0_8,0.0_8,1.0_8),p010=point(0.0_8&
       &,1.0_8,0.0_8),p100=point(1.0_8,0.0_8,0.0_8)
  Type(RecFace), SAVE, Target :: Box(__MaxSide)
  Type(Face), SAVE, Target :: Faces(__MaxSide)
  Type(point), SAVE :: MyCenter,ItsCenter
  Type(point), PARAMETER :: p_zero=Point(0.0_8,0.0_8,0.0_8)
CONTAINS
  Subroutine BoxGeometry_(coa, oca)
    Integer :: Dir,Axis !- Axis is in the scaled frame
    Real(8) :: coa(3,3), oca(3,3)
    Type(point) :: p1,p2,p3,p4
    Type(point) :: c0,g1,g2,g3,g4,g5,n1,n2,n3,n4,n5,s1,s2,s3,s4
    Type(Line)  :: l1,l2,l3,l4,l5
    Type(Plane) :: pL
    Real(8) :: aux1,Norm,vout
    Integer :: n

!!$--- Initialize the My Box, the box the CPU is working on

    co=coa
    oc=oca

    MyOrigin=Point(Dble(PI__Ranks(PI_Node_Cart+1) % nx),&
         & Dble(PI__Ranks(PI_Node_Cart+1) % ny)&
         &,Dble(PI__Ranks(PI_Node_Cart+1) % nz)) 
    MySides=Point(ddx,ddy,ddz)
    MyOrigin=MyOrigin*MySides

    Box(1)%Dir=-1; Box(1)%Axis=1; Box(1)%d=MyOrigin%x+__Segment(Box(1)%Dir)
    Box(2)%Dir=1 ; Box(2)%Axis=1; Box(2)%d=MyOrigin%x+__Segment(Box(2)%Dir)
    
    Box(3)%Dir=-1; Box(3)%Axis=2; Box(3)%d=MyOrigin%y+__Segment(Box(3)%Dir)
    Box(4)%Dir=1 ; Box(4)%Axis=2; Box(4)%d=MyOrigin%y+__Segment(Box(4)%Dir)
    
    Box(5)%Dir=-1; Box(5)%Axis=3; Box(5)%d=MyOrigin%z+__Segment(Box(5)%Dir)
    Box(6)%Dir=1 ; Box(6)%Axis=3; Box(6)%d=MyOrigin%z+__Segment(Box(6)%Dir)
    MyCenter=MyOrigin+(__Half)*MySides
    
    Do n=1,6
       Faces(n)%Dir=Box(n)%Dir
       Faces(n)%Axis=Box(n)%Axis
       aux1=__Segment(Box(n)%Dir)
       
       If(aux1 == 0.0_8) THEN
          select Case(Box(n)%Axis)
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
          select Case(Box(n)%Axis)
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
       s1=co*p1;s2=co*p2;s4=co*p4
       g1=s2-s1; g2=s4-s1
       __normalize_v(g1,Norm) ; __normalize_v(g2,Norm) 
       Norm=__dotprod_v(g1,g2)
       Norm=Acos(Norm)*180.0D0/Pi
       IF(Norm < 90.0D0) THEN
          s1=p1; p1=p4; p4=p3;p3=p2;p2=s1
       END IF

       c0=p1+p2+p3+p4
       c0=0.25_8*c0
 
       l1%a=p1; l1%b=p2-p1
       l2%a=p2; l2%b=p3-p2
       l3%a=p3; l3%b=p4-p3
       l4%a=p4; l4%b=p1-p4
       l5%a=p1; l5%b=p3-p1

       g1=l1%a+l1%b*__Half
       n1=c0+(g1-c0)*(1.0D0+__Fraction)

       g2=l2%a+l2%b*__Half
       n2=c0+(g2-c0)*(1.0D0+__Fraction)

       g3=l3%a+l3%b*__Half
       n3=c0+(g3-c0)*(1.0D0+__Fraction)

       g4=l4%a+l4%b*__Half
       n4=c0+(g4-c0)*(1.0D0+__Fraction)

       n5=l5%a+l5%b*(__Half+__Fraction)
       l5%a=p4; l5%b=p2-p4

       
       p1=co*p1; p2=co*p2; p3=co*p3; p4=co*p4; 
       n1=co*n1; n2=co*n2; n3=co*n3; n4=co*n4; n5=co*n5;
       l1%a=co*l1%a; l2%a=co*l2%a; l3%a=co*l3%a; l4%a=co*l4%a; l5%a=co*l5%a; 
       l1%b=co*l1%b; l2%b=co*l2%b; l3%b=co*l3%b; l4%b=co*l4%b; l5%b=co*l5%b; 

       g1=Foot(n1,l1)
       g2=Foot(n2,l2)
       g3=Foot(n3,l3)
       g4=Foot(n4,l4)
       g5=Foot(n5,l5)

       n1=Versor(g1,n1)
       n2=Versor(g2,n2)
       n3=Versor(g3,n3)
       n4=Versor(g4,n4)
       n5=Versor(g5,n5)

       pL=Equation_Plane(p1,p2,p4)

       Faces(n)%Dir=Box(n)%Dir
       Faces(n)%Axis=Box(n)%Axis

       Faces(n)%a=p1
       Faces(n)%b=p2
       Faces(n)%c=p3
       Faces(n)%d=p4

       Faces(n)%ab=n1
       Faces(n)%bc=n2
       Faces(n)%cd=n3
       Faces(n)%da=n4
       Faces(n)%db=n5

       Faces(n)%lab=l1
       Faces(n)%lbc=l2
       Faces(n)%lcd=l3
       Faces(n)%lda=l4
       Faces(n)%m_ab=__dotprod_v(l1%b,l1%b)
       Faces(n)%m_bc=__dotprod_v(l2%b,l2%b)
       Faces(n)%m_cd=__dotprod_v(l3%b,l3%b)
       Faces(n)%m_da=__dotprod_v(l4%b,l4%b)
       
       Faces(n)%pL=pL
    End Do
  End Subroutine BoxGeometry_
  Function Distance_(xa,ya,za,ngroupa,Dir,Axis,dest) Result(out)
    Integer, Optional :: Dir, Axis,dest
    Real(8) :: out
    Real(8) :: xa,ya,za
    Integer :: ngroupa

    Integer, SAVE :: MyFaceNo
    Type(point), SAVE, Target :: pj,p1a,n1,n2,n3,n4,n5,p1,p2,p3,p4
    Type(Line), SAVE, Target :: l1,l2,l3,l4
    Real(8), SAVE, Target :: m_ab,m_bc,m_cd,m_da
    Type(RecFace), Pointer,Save :: MyBox
    Type(Face), Pointer, Save :: MyFace

    Type(point) :: p,p_a,tmp,vout,tmp1,tmp2,tmp3,tmp4
    Real(8) :: Dist_Pl, mag_vab, mag_vad, dot_vab, dot_vad
    Logical :: InsideFace
    Integer, Save :: pipo=0
    Type(point), Pointer :: pa,pb,pc,nac,nbc
    Type(point) :: pcp,pcb
    Type(line),  Pointer :: lac,lbc
    Real(8),     Pointer :: p_bc,p_ac
    Real(8) :: dpnbc,dbnbc,dcpcb,dpnac,dcnac,dpn5,dp2n5
    Type(point) :: dummy
    Real(8) :: Dists(27)
    Integer :: Region
    Integer, Save :: Axis_a,dest_a,dir_a
    Integer :: i,j,k,ia,ja,ka
    Real(8) :: dist_min
    Real(8) :: ang1,ang2,ang3
    Integer :: nx,ny,nz,numcell
    Real(8) :: x1,y1,z1,dx,dy,dz

    out=-1.0_8
    
    If(Present(Dir)) Then
       Axis_a=Axis
       dir_a=dir
       dest_a=dest
       IF(FindFace_(Dir,Axis,dest) < 0) Then
          out=-1.0D0
       End IF
       out=0.0D0
       Return
    End If
    
    If(Abs(co(1,2)) > 1.0D-6 .And. Abs(co(1,3)) > 1.0D-6 .And.&
         & Abs(co(2,3)) > 1.0D-6) Then
       out=1.0; Return
    End If

    p_a%x=xa; p_a%y=ya; p_a%z=za;       

    tmp=p_a-ItsCenter
    __Pbc(tmp,vout)
    p_a=p_a+vout


    p=co*p_a
    pj=FootDist_Plane(MyFace%pL,p, Dist_pl)


    dpn5=__dotprod_v(pj,n5)
    dp2n5=__dotprod_v(p2,n5)
    If(dpn5 < dp2n5) Then
       pa=>p2; pb=>p4; pc=>p1
       nac=>n1; nbc=>n4
       lac=>l1; lbc=>l4
       p_ac=>m_ab; p_bc=>m_da
    Else
       pa=>p4; pb=>p2; pc=>p3
       nac=>n3; nbc=>n2
       lac=>l3; lbc=>l2
       p_ac=>m_cd; p_bc=>m_bc
    End If
    
    dpnbc=__dotprod_v(pj,nbc)
    dbnbc=__dotprod_v(pb,nbc)
    Region=0
    If(dpnbc > dbnbc) Then
       pcp=pj-pc; pcb=pb-pc
       dcpcb=__dotprod_v(pcp,pcb)
       If(dcpcb < 0.0D0) Then
          Region=3
       Else
          Region=2
       End If
    Else
       dpnac=__dotprod_v(pj,nac)
       dcnac=__dotprod_v(pc,nac)
       If(dpnac > dcnac) Then
          Region=4
       Else
          Region=5
       End If
    End If

    Select Case(Region)
    Case(2)
       out=Region_2()
    Case(3)
       out=Region_3()
    Case(4)
       out=region_4()
    Case(5)
       out=Dist_pl
    Case(0)
       Errmsg_f='BoxGeometry Fails: Can''t attribute group to region'
       Call Add_Errors(-1,errmsg_f)
       Call Print_Errors()
    End Select
       
  Contains
    Function FindFace_(Dir,Axis,dest) Result(out)
      Integer :: Dir, Axis, Dest, out
      Integer :: n
      Logical :: ok
      Type(point) :: ItsOrigin
      Real(8) :: vv


      out=-1
      errmsg_f=' '
      ok=.FALSE.
      Do n=1,__MaxSide
         If(Dir == Box(n) % Dir .And. Axis == Box(n) % Axis) Then
            out=n
            ok=.True.
         End If
      End Do
      If(out == -1) Then
         Errmsg_f='Internal Error From BoxGeometry::FindFace: can find f&
              &ace in the list'
         Call Add_Errors(-1,Errmsg_f)
      End If
      Call PI_ErrSignal_(errmsg_f,ok)
      IF(.Not. ok) Then
         out=-1
         Return
      End IF

      ItsOrigin=Point(Dble(PI__Ranks(dest+1) % nx),&
           & Dble(PI__Ranks(dest+1) % ny)&
           &,Dble(PI__Ranks(dest+1) % nz)) 
      ItsOrigin=ItsOrigin*MySides

      Select Case(Axis) 
      Case(_x_)
         vv=ItsOrigin%x-MyOrigin%x
         If(vv*Dble(Dir) < 0.0D0) ItsOrigin%x=ItsOrigin%x+Dble(Dir)*2.0D0
      Case(_y_)
         vv=ItsOrigin%y-MyOrigin%y
         If(vv*Dble(Dir) < 0.0D0) ItsOrigin%y=ItsOrigin%y+Dble(Dir)*2.0D0
      Case(_z_)
         vv=ItsOrigin%z-MyOrigin%z
         If(vv*Dble(Dir) < 0.0D0) ItsOrigin%z=ItsOrigin%z+Dble(Dir)*2.0D0
      End Select

      ItsCenter=ItsOrigin+(__Half)*MySides

      MyFaceNo=out
      MyBox=>Box(MyFaceNo)
      MyFace=>Faces(MyFaceNo)
      n1=Faces(MyFaceNo)%ab
      n2=Faces(MyFaceNo)%bc
      n3=Faces(MyFaceNo)%cd
      n4=Faces(MyFaceNo)%da
      n5=Faces(MyFaceNo)%db

      p1=Faces(MyFaceNo)%a
      p2=Faces(MyFaceNo)%b
      p3=Faces(MyFaceNo)%c
      p4=Faces(MyFaceNo)%d
      
      l1=Faces(MyFaceNo)%lab
      l2=Faces(MyFaceNo)%lbc
      l3=Faces(MyFaceNo)%lcd
      l4=Faces(MyFaceNo)%lda

      m_ab=Faces(MyFaceNo)%m_ab
      m_bc=Faces(MyFaceNo)%m_bc
      m_cd=Faces(MyFaceNo)%m_cd
      m_da=Faces(MyFaceNo)%m_da

    End Function FindFace_
    Function Region_2() Result(out)
      Real(8) :: out
      Real(8) :: rsq,iop
      Type(point) :: MyP
      Type(line) ::  MyL
      out=-1.0D0
      If(dcpcb > p_bc) Then
         out=__dist_v(p,pb)
      Else
         MyP=Foot(pj,lbc)
         rsq=__dist2_v(pj,MyP)
         out=Sqrt(rsq+dist_pl**2)
      End If
    End Function Region_2
    Function Region_3() Result(out)
      Real(8) :: out
      Real(8) :: dcpca
      Real(8) :: rsq
      Type(point) :: MyP

      out=-1.0D0
      dcpca=__dotprod_v(pcp,lac%b)
      If(dcpca < 0.0D0) Then
         out=__dist_v(p,pc)
      Else
         MyP=Foot(pj,lac)         
         rsq=__dist2_v(pj,MyP)
         out=Sqrt(rsq+dist_pl**2)
      End If
    End Function Region_3
    Function Region_4() Result(out)
      Real(8) :: out
      Real(8) :: rsq
      Real(8) :: dcpca
      Type(point) :: MyP

      out=-1.0D0
      dcpca=__dotprod_v(pcp,lac%b)
      If(dcpca > p_ac) Then
         out=__dist_v(p,pa)
      Else
         MyP=Foot(pj,lac)
         rsq=__dist2_v(pj,MyP)
         out=Sqrt(rsq+dist_pl**2)
      End If
    End Function Region_4
  End Function Distance_
END MODULE BoxGeometry
