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
MODULE Decompose
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Jan 26 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*

#ifdef PARALLEL
  USE mpi
#endif
  USE PI_
  USE Node
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f
  USE Print_Defs
  Use FactorizeNo
  IMPLICIT none
  PRIVATE
  PUBLIC Decompose_,AssignAtomsToCells,Indbox_g_p, Indbox_a_p
  Integer, Allocatable, Save :: IndBox_g_p(:),IndBox_a_p(:)
CONTAINS
  SUBROUTINE Decompose_(a,b,c,npxa,npya,npza)
    Real(8) :: a,b,c
    INTEGER, OPTIONAL :: npxa,npya,npza
    INTEGER, POINTER :: vect(:)=>NULL()
    LOGICAL, ALLOCATABLE :: Mask(:)
    INTEGER :: o0, Ic,  Jc, Kc, n, count0, m, l
    REAL(8) :: cube
    INTEGER, ALLOCATABLE :: combs(:)
    INTEGER, POINTER :: vect1(:)=>NULL()
    LOGICAL :: ok
    INTEGER :: n_combs, m_combs, l_combs
    REAL(8) :: d_x, d_Min
    REAL(8) :: xyz(3,8) !-- vertices
    REAL(8) :: dx, dy, dz, fx, fy, fz, gx, gy, gz
    INTEGER :: n1,n2,n3

    IF(PRESENT(npxa)) THEN
       PI_npx=npxa
       PI_npy=npya
       PI_npz=npza
       IF(PI_nprocs /= PI_npx*PI_npy*PI_npz) THEN
          errmsg_f='PI Geometry does not conform to the number of processors '
          CALL Add_Errors(-1,errmsg_f)
          CALL Print_Errors()
       END IF
    ELSE
       vect=>FactorizeNo_(PI_nprocs, 7)
       Ic=NINT(((a**2/(b*c))*DBLE(PI_nprocs))**(1.0D0/3.0D0))
       Jc=NINT(Ic*b/a)
       Kc=NINT(Ic*c/a)
       o0=PRODUCT(vect)
       ALLOCATE(Mask(SIZE(vect)))
       Mask=.TRUE.
       IF(.NOT. Node_()) STOP
       CALL Combinations(o0, Mask, vect)
       
       count0=Node__Size()
       ALLOCATE(combs(count0))
       count0=0
       DO WHILE(Node__Pop(vect1))
          count0=count0+1
          ok=.TRUE.
          DO n=1,count0-1
             IF(combs(n) == vect1(1)) THEN
                count0=count0-1
                ok=.FALSE.
             END IF
          END DO
          IF(.NOT. ok) CYCLE
          combs(count0)=vect1(1)
       END DO
       d_Min=1.0D10
       DO l=1,count0
          DO m=1,count0
             DO n=1,count0
                IF(combs(n)*combs(m)*combs(l) == PI_nprocs) THEN
                   d_x=DSQRT(DBLE(Ic-combs(n))**2+DBLE(Jc-combs(m))**2+DBLE(Kc-combs(l))**2)
                   IF(d_x < d_Min) THEN
                      D_min=d_x
                      n_combs=combs(n)
                      m_combs=combs(m)
                      l_combs=combs(l)
                   END IF
                END IF
             END DO
          END DO
       END DO
       PI_npx=n_combs; PI_npy=m_combs; PI_npz=l_combs
    END IF
    IF(PI_Npx /= 1) THEN
       IF(PI_npy  == 1) THEN
          PI_npy=PI_npx
          PI_npx=1
       ELSE IF(PI_npz  == 1) THEN
          PI_npz=PI_npx
          PI_npx=1
       END IF
    END IF

    WRITE(kprint,100) PI_nprocs,PI_npx,PI_npy,PI_npz
    CALL PI__Setup_Cart
100 FORMAT(' ====> Running with ',i3,' Processors ',' nx = ',i2,' ny =&
         & ',i2,' nz = ',i2,' <====')
  CONTAINS
    RECURSIVE SUBROUTINE Combinations(Prod, Mask, vect)
      INTEGER :: vect(:)
      INTEGER :: Prod
      LOGICAL :: mask(:)
      INTEGER :: i
      LOGICAL, ALLOCATABLE :: mask_new(:)
      INTEGER :: Prod_New
      
      ALLOCATE(Mask_New(SIZE(Mask)) )
      mask_new=Mask

      DO i=1,SIZE(vect)
         IF(Mask(i)) THEN
            Mask_New(i)=.FALSE.
            Prod_New=Prod/vect(i)
            CALL Node__Push(Prod_New)
            Call Combinations(Prod_new, Mask_New, vect)
         END IF
      END DO
    END SUBROUTINE Combinations
  END SUBROUTINE Decompose_
  SUBROUTINE AssignAtomsToCells(natom,ngrp,grppt,xa,ya,za)
    Integer :: ngrp,grppt(2,*),natom
    Real(8) :: xa(*),ya(*),za(*)
    INTEGER :: n,m,nx,ny,nz,numcell,l,nmin,ntap&
         &,cp_0,cp_3, AtSt, AtEn,mm,ox,oy,oz,mpe,mp&
         &,count0,count1,i,numcell1,ncx,ncy,ncz
    REAL(8) :: x1,y1,z1,dx,dy,dz
    Integer, Allocatable :: Ind0(:),ind1(:)
    INTEGER :: i_p

    ncx=PI_npx
    ncy=PI_npy
    ncz=PI_npz

    dx=2.d0/DBLE(ncx)
    dy=2.d0/DBLE(ncy)
    dz=2.d0/DBLE(ncz)

    m=PI_Node_Cart+1
    numcell=PI__Ranks(m) % n

    Allocate(Ind0(ngrp),Ind1(natom))
    count0=0
    count1=0
    DO n=1,ngrp
       x1=xa(n)/dx
       y1=ya(n)/dy
       z1=za(n)/dz
       nx=INT(x1)+(SIGN(1.D0,x1-INT(x1))-1.)/2
       ny=INT(y1)+(SIGN(1.D0,y1-INT(y1))-1.)/2
       nz=INT(z1)+(sign(1.d0,z1-int(z1))-1.)/2
       nx=MOD(MOD(nx,ncx)+ncx,ncx)
       ny=MOD(MOD(ny,ncy)+ncy,ncy)
       nz=MOD(MOD(nz,ncz)+ncz,ncz)
       numcell1=nz+ncz*(ny+ncy*nx)+1
       IF(numcell == numcell1) THEN
          count0=count0+1
          Ind0(count0)=n
          Do i=grppt(1,n),grppt(2,n)
             count1=count1+1
             Ind1(count1)=i
          End Do
       END IF
    END DO
    Allocate(IndBox_g_p(count0),IndBox_a_p(count1))
    IndBox_g_p=Ind0(1:count0)
    IndBox_a_p=Ind1(1:count1)
    n=count0
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,n,1,MPI_INTEGER4,MPI_SUM,PI_Comm_Cart,ierr)
    WRITE(kprint,*) 'ngroup ',n
    n=count1
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,n,1,MPI_INTEGER4,MPI_SUM,PI_Comm_Cart,ierr)
    WRITE(kprint,*) 'natom ',n

  END SUBROUTINE AssignAtomsToCells
END MODULE Decompose
