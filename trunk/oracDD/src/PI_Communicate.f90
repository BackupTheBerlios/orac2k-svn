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
  USE NeighCells
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f
  IMPLICIT none
  PRIVATE
  PUBLIC PI__PickDecomposition, PI__Topology, pe, PI__GetParameters
  INTEGER, SAVE :: nprocs,npx,npy,npz
  TYPE :: PI__Topology
     INTEGER :: nx,ny,nz !-- 3-D label
  END TYPE PI__Topology
  TYPE(PI__Topology), ALLOCATABLE, SAVE :: pe(:)
  REAL(8), PARAMETER :: nv(8,3)=RESHAPE((/&
       & 0.0D0, 1.0D0, 0.0D0, 0.0D0, 1.0D0, 1.0D0, 0.0D0, 1.0D0 &
       &,0.0D0, 0.0D0, 1.0D0, 0.0D0, 1.0D0, 0.0D0, 1.0D0, 1.0D0 &
       &,0.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0, 1.0D0, 1.0D0, 1.0D0 &
       &/),(/8, 3/))
CONTAINS
  SUBROUTINE PI__PickDecomposition(rcut,nprocsa, npxa,npya,npza)
    REAL(8) :: rcut
    INTEGER, OPTIONAL :: npxa,npya,npza
    INTEGER :: nprocsa
    REAL(8) :: rcuut
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

    nprocs=nprocsa
    IF(PRESENT(npxa)) THEN
       npx=npxa
       npy=npya
       npz=npza
       nprocs=npx*npy*npz
    ELSE
       vect=>FactorizeNo_(nprocs)
       Ic=NINT(((a**2/(b*c))*DBLE(nprocs))**(1.0D0/3.0D0))
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
                IF(combs(n)*combs(m)*combs(l) == nprocs) THEN
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
    WRITE(*,*) ' nprocs, nx, ny, nz ',nprocs,npx,npy,npz
    ALLOCATE(Pe(nprocs))
    dx=boxl/DBLE(npx)
    dy=boxl/DBLE(npy)
    dz=boxl/DBLE(npz)
    count0=0
    DO n1=1,npx
       DO n2=1,npy
          DO n3=1,npz
             count0=count0+1
             Pe(count0) % nx=n1
             Pe(count0) % ny=n2
             Pe(count0) % nz=n3
          END DO
       END DO
    END DO
    IF(.NOT. NeighCells_(rcut,nprocs,npx,npy,npz)) CALL Print_Errors()
  CONTAINS
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
  SUBROUTINE PI__GetParameters(nprocsa,npxa,npya,npza)
    INTEGER :: nprocsa,npxa,npya,npza
    nprocsa=nprocs;npxa=npx;npya=npy;npza=npz
  END SUBROUTINE PI__GetParameters
END MODULE PI_Communicate
