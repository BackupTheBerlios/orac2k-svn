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
  
  USE Node
  USE FactorizeNo
  USE Cell
  IMPLICIT none
  PRIVATE
  PUBLIC PI__PickDecomposition
  INTEGER, SAVE :: nproc,npx,npy,npz
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
    LOGICAL :: ok
    INTEGER :: n_combs, m_combs, l_combs
    REAL(8) :: d_x, d_Min

    nproc=nproca
    IF(PRESENT(npxa)) THEN
       npx=npxa
       npy=npya
       npz=npza
       nproc=npx*npy*npz
       RETURN
    END IF
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
    WRITE(*,*) n_combs,m_combs,l_combs
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
END MODULE PI_Communicate
