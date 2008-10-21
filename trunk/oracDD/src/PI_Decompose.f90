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
MODULE PI_Decompose
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

#ifdef HAVE_MPI
  USE mpi
#endif
  USE PI_
  USE PI_Communicate
  USE UNITS
  USE Node
  USE FactorizeNo
  USE Cell
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f
  USE Print_Defs
  USE PI_Neighbors
  USE Atom
  USE Groups
  IMPLICIT none
  PRIVATE
  PUBLIC PI__Decomposition_NB, PI__GetParameters, PI__Type, PI__AssignAtomsToCells
  TYPE :: PI__Type
     LOGICAL :: ok=.FALSE.
     INTEGER, ALLOCATABLE :: exc(:,:)
  END TYPE PI__Type

  REAL(8), SAVE :: Ref_Radius=3.11D0
  REAL(8), PARAMETER :: nv(8,3)=RESHAPE((/&
       & 0.0D0, 1.0D0, 0.0D0, 0.0D0, 1.0D0, 1.0D0, 0.0D0, 1.0D0 &
       &,0.0D0, 0.0D0, 1.0D0, 0.0D0, 1.0D0, 0.0D0, 1.0D0, 1.0D0 &
       &,0.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0, 1.0D0, 1.0D0, 1.0D0 &
       &/),(/8, 3/))
CONTAINS
  SUBROUTINE PI__Decomposition_nb(rcut,npxa,npya,npza)
    REAL(8) :: rcut(:)
    INTEGER, OPTIONAL :: npxa,npya,npza
    INTEGER, POINTER :: vect(:)=>NULL()
    LOGICAL, ALLOCATABLE :: Mask(:)
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
       vect=>FactorizeNo_(PI_nprocs)
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
    CALL PI__Setup_SndRcv

100 FORMAT(' ====> Running with ',i3,' Processors ',' nx = ',i2,' ny = ',i2,' nz = ',i2,' <====')
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
  END SUBROUTINE PI__Decomposition_nb
  SUBROUTINE PI__GetParameters(nprocsa,npxa,npya,npza)
    INTEGER :: nprocsa,npxa,npya,npza
    nprocsa=PI_nprocs; npxa=PI_npx; npya=PI_npy; npza=PI_npz
  END SUBROUTINE PI__GetParameters
  SUBROUTINE PI__AssignAtomsToCells
    INTEGER :: npx,npy,npz,n,m,nx,ny,nz,mx,my,mz,numcell,l,nmin,ntap&
         &,ngrp,cp_0,cp_3, AtSt, AtEn,mm,ox,oy,oz,mpe,mp&
         &,count0
    LOGICAL, ALLOCATABLE :: Mask(:)
    INTEGER :: i_p

    npx=PI_npx
    npy=PI_npy
    npz=PI_npz

    ntap=SIZE(Atoms)
    ngrp=SIZE(Groupa)

    groupa(:) % knwn=3

    IF(.NOT. PI_Neighbors_(Groupa(:) % xa, Groupa(:) % ya, Groupa(:) &
         &% za, Groupa(:) % knwn)) CALL Print_Errors()

    m=PI_Node_Cart+1
    mx=PI__Ranks(m) % nx + 1
    my=PI__Ranks(m) % ny + 1
    mz=PI__Ranks(m) % nz + 1
    numcell=PI__Ranks(m) % n

    l=Head_xyz(numcell)
    nmin=0
    count0=0
    DO WHILE(l > nmin)
       AtSt=Groupa(l) % AtSt
       AtEn=Groupa(l) % AtEn
       count0=count0+AtEn-AtSt+1
       Groupa(l) % knwn = 1
       l=Chain_xyz(l) % p
    END DO

    DO n=1,SIZE(Groupa)
       IF(Groupa(n) % knwn /= 1 ) THEN
          Groupa(n) % knwn = 0
          Groupa(n) % xa=0.0D0
          Groupa(n) % ya=0.0D0
          Groupa(n) % za=0.0D0
          Groupa(n) % x=0.0D0
          Groupa(n) % y=0.0D0
          Groupa(n) % z=0.0D0
          AtSt=Groupa(n) % AtSt
          AtEn=Groupa(n) % AtEn
          DO mm=AtSt,AtEn
             Atoms(mm) % x = 0.0D0
             Atoms(mm) % y = 0.0D0
             Atoms(mm) % z = 0.0D0
             Atoms(mm) % xa = 0.0D0
             Atoms(mm) % ya = 0.0D0
             Atoms(mm) % za = 0.0D0
          END DO
       END IF          
    END DO

  END SUBROUTINE PI__AssignAtomsToCells
END MODULE PI_Decompose
