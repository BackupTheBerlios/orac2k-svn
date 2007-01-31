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
MODULE PI_Communicate
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
  
  USE UNITS
  USE Node
  USE FactorizeNo
  USE Cell
  IMPLICIT none
  PRIVATE
  PUBLIC PI__PickDecomposition, Plane, PI__Topology, pe, PI__GetParameters
  INTEGER, SAVE :: nproc,npx,npy,npz
  TYPE :: Plane
     REAL(8) :: nm(3)  !-- Normal to the planes      
     REAL(8) :: x,y,z
  END TYPE Plane
  TYPE :: PI__Topology
     INTEGER :: nx,ny,nz !-- 3-D label
     TYPE(Plane) :: pln(6)
  END TYPE PI__Topology
  TYPE(PI__Topology), ALLOCATABLE, SAVE :: pe(:)
  REAL(8), PARAMETER :: nv(8,3)=RESHAPE((/&
       & 0.0D0, 1.0D0, 0.0D0, 0.0D0, 1.0D0, 1.0D0, 0.0D0, 1.0D0 &
       &,0.0D0, 0.0D0, 1.0D0, 0.0D0, 1.0D0, 0.0D0, 1.0D0, 1.0D0 &
       &,0.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0, 1.0D0, 1.0D0, 1.0D0 &
       &/),(/8, 3/))
CONTAINS
  SUBROUTINE PI__PickDecomposition(nproca, npxa,npya,npza)
    INTEGER, OPTIONAL :: npxa,npya,npza
    INTEGER :: nproca
    INTEGER, POINTER :: vect(:)=>NULL()
    LOGICAL, POINTER :: Mask(:)
    INTEGER :: o0, Ic,  Jc, Kc, n, count0, m, l
    REAL(8) :: cube
    INTEGER, ALLOCATABLE :: combs(:)
    INTEGER, POINTER :: vect1(:)=>NULL()
    REAL(8), POINTER :: vec(:,:)=> NULL()
    LOGICAL :: ok
    INTEGER :: n_combs, m_combs, l_combs
    REAL(8) :: d_x, d_Min
    REAL(8) :: xyz(3,8) !-- vertices
    REAL(8) :: dx, dy, dz, fx, fy, fz, gx, gy, gz
    INTEGER :: n1,n2,n3

    nproc=nproca
    IF(PRESENT(npxa)) THEN
       npx=npxa
       npy=npya
       npz=npza
       nproc=npx*npy*npz
    ELSE
       vect=>FactorizeNo_(nproc)
       Ic=NINT(((a**2/(b*c))*DBLE(nproc))**(1.0D0/3.0D0))
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
       DO n=1,count0
          DO m=1,count0
             DO l=1,count0
                IF(combs(n)*combs(m)*combs(l) == nproc) THEN
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
       npx=n_combs; npy=m_combs; npz=l_combs
    END IF
    WRITE(*,*) ' nproc, nx, ny, nz ',nproc,npx,npy,npz
    ALLOCATE(Pe(nproc))
    dx=boxl/DBLE(npx)
    dy=boxl/DBLE(npy)
    dz=boxl/DBLE(npz)
    count0=0
    ALLOCATE(vec(3,6))
    DO n1=1,npx
       fx=DBLE(n1-1)*dx
       DO n2=1,npy
          fy=DBLE(n2-1)*dy
          DO n3=1,npz
             fz=DBLE(n3-1)*dz
             count0=count0+1
             Pe(count0) % nx=n1
             Pe(count0) % ny=n2
             Pe(count0) % nz=n3
             DO m=1,8
                gx = fx + nv(m,1)*dx
                gy = fy + nv(m,2)*dy
                gz = fz + nv(m,3)*dz
                xyz(1,m) = co(1,1)*gx+co(1,2)*gy+co(1,3)*gz
                xyz(2,m) = co(2,1)*gx+co(2,2)*gy+co(2,3)*gz
                xyz(3,m) = co(3,1)*gx+co(3,2)*gy+co(3,3)*gz
             END DO
             vec(:,1)=xyz(:,2)-xyz(:,1)
             vec(:,2)=xyz(:,3)-xyz(:,1)
             vec(:,3)=xyz(:,4)-xyz(:,1)
             vec(:,4)=xyz(:,5)-xyz(:,8)
             vec(:,5)=xyz(:,6)-xyz(:,8)
             vec(:,6)=xyz(:,7)-xyz(:,8)
             vec=NormVect(vec)
             pe(count0) % pln(1) % nm(:)=CrossProd(vec(:,2),vec(:,1))             
             pe(count0) % pln(1) % x = xyz(1,1)
             pe(count0) % pln(1) % y = xyz(2,1)
             pe(count0) % pln(1) % z = xyz(3,1)
             pe(count0) % pln(2) % nm(:)=CrossProd(vec(:,1),vec(:,3))
             pe(count0) % pln(2) % x = xyz(1,1)
             pe(count0) % pln(2) % y = xyz(2,1)
             pe(count0) % pln(2) % z = xyz(3,1)
             pe(count0) % pln(3) % nm(:)=CrossProd(vec(:,3),vec(:,2))
             pe(count0) % pln(3) % x = xyz(1,1)
             pe(count0) % pln(3) % y = xyz(2,1)
             pe(count0) % pln(3) % z = xyz(3,1)
             pe(count0) % pln(4) % nm(:)=CrossProd(vec(:,5),vec(:,4))
             pe(count0) % pln(4) % x = xyz(1,8)
             pe(count0) % pln(4) % y = xyz(2,8)
             pe(count0) % pln(4) % z = xyz(3,8)
             pe(count0) % pln(5) % nm(:)=CrossProd(vec(:,6),vec(:,5))
             pe(count0) % pln(5) % x = xyz(1,8)
             pe(count0) % pln(5) % y = xyz(2,8)
             pe(count0) % pln(5) % z = xyz(3,8)
             pe(count0) % pln(6) % nm(:)=CrossProd(vec(:,4),vec(:,6))
             pe(count0) % pln(6) % x = xyz(1,8)
             pe(count0) % pln(6) % y = xyz(2,8)
             pe(count0) % pln(6) % z = xyz(3,8)
          END DO
       END DO
    END DO
  CONTAINS
    FUNCTION NormVect(vect) RESULT(out)
      REAL(8) :: vect(:,:)
      REAL(8), POINTER :: out(:,:)
      REAL(8) :: norm
      INTEGER :: n

      ALLOCATE(out(SIZE(vect,1),SIZE(vect,2)))
      DO n=1,SIZE(vect,2)
         norm=SUM(vect(:,n)**2)
         out(:,n)=vect(:,n)/SQRT(norm)
      END DO
    END FUNCTION NormVect
    FUNCTION CrossProd(v1,v2) RESULT(out)
      REAL(8) :: v1(3),v2(3)
      REAL(8) :: out(3)
      out(1)=v1(2)*v2(3)-v1(3)*v2(2)
      out(2)=v1(3)*v2(1)-v1(1)*v2(3)
      out(3)=v1(1)*v2(2)-v1(2)*v2(1)
      out=out/SQRT(SUM(out**2))
    END FUNCTION CrossProd
    RECURSIVE SUBROUTINE Combinations(Prod, Mask, vect)
      INTEGER :: vect(:)
      INTEGER :: Prod
      LOGICAL :: mask(:)
      INTEGER :: i
      LOGICAL, POINTER :: mask_new(:)
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
  END SUBROUTINE PI__PickDecomposition
  SUBROUTINE PI__GetParameters(nproca,npxa,npya,npza)
    INTEGER :: nproca,npxa,npya,npza
    nproca=nproc;npxa=npx;npya=npy;npza=npz
  END SUBROUTINE PI__GetParameters
END MODULE PI_Communicate
