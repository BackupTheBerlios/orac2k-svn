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
  USE UNITS
  USE Node
  USE FactorizeNo
  USE Cell
  USE NeighCells
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f
  USE Print_Defs
  IMPLICIT none
  PRIVATE
  PUBLIC PI__Decomposition_NB, PI__Topology, PI__GetParameters, PI__Type
  TYPE :: PI__Type
     LOGICAL :: ok=.FALSE.
     INTEGER, ALLOCATABLE :: exc(:,:)
  END TYPE PI__Type
  TYPE :: PI__Topology
     TYPE(PI__Type), ALLOCATABLE :: tbl(:,:)
  END TYPE PI__Topology
  TYPE(PI__Topology), ALLOCATABLE, SAVE :: pe_nb(:)
  TYPE(PI__Topology), SAVE :: pe_bb

  REAL(8), PARAMETER :: nv(8,3)=RESHAPE((/&
       & 0.0D0, 1.0D0, 0.0D0, 0.0D0, 1.0D0, 1.0D0, 0.0D0, 1.0D0 &
       &,0.0D0, 0.0D0, 1.0D0, 0.0D0, 1.0D0, 0.0D0, 1.0D0, 1.0D0 &
       &,0.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0, 1.0D0, 1.0D0, 1.0D0 &
       &/),(/8, 3/))
CONTAINS
  SUBROUTINE PI__Decomposition_nb(rcut,npxa,npya,npza)
    REAL(8) :: rcut(3)
    INTEGER, OPTIONAL :: npxa,npya,npza
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
    WRITE(kprint,*) ' PI_nprocs, nx, ny, nz ',PI_nprocs,PI_npx,PI_npy,PI_npz

    CALL PI__Get_Ranks

    ALLOCATE(Pe_nb(3))
    DO n=1,3
       IF(rcut(n) > 0.0D0) THEN
          IF(.NOT. NeighCells_(rcut(1),rcut(n),PI_nprocs,PI_npx,PI_npy,PI_npz))&
               & CALL Print_Errors()
          CALL Make_Comm(n)
       END IF
    END DO
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
    SUBROUTINE Make_Comm(i_n)
      INTEGER :: i_n
      INTEGER :: n,nx,ny,nz,mp,m,np,i,j,k,ox,oy,oz&
           &,mpe,ncx,ncy,ncz,count0
      INTEGER, POINTER :: ind_x(:)=>NULL(),ind_y(:)=>NULL()

      np=PI_nprocs
      IF(ALLOCATED(Pe_nb(i_n) % tbl)) DEALLOCATE(Pe_nb(i_n) % tbl)
      ALLOCATE(Pe_nb(i_n) % tbl(np,np))
      CALL NeighCells__Param(mp,ncx,ncy,ncz)
      WRITE(kprint,*) ' mp, ncx, ncy, ncz ',mp,ncx,ncy,ncz
      ALLOCATE(ind_x(np), ind_y(mp))
      DO n=1,np
         ind_x=0
         ind_y=0
         mp=SIZE(Nei(n) % c,2)
         DO m=1,mp
            ox=Nei(n) % c(1,m)
            oy=Nei(n) % c(2,m)
            oz=Nei(n) % c(3,m)
            i=Ind_Small(ox,oy,oz) % nx - 1
            j=Ind_Small(ox,oy,oz) % ny - 1
            k=Ind_Small(ox,oy,oz) % nz - 1

            mpe=k+PI_npz*(i*PI_npy+j)+1
            ind_y(m)=mpe
            ind_x(mpe)=ind_x(mpe)+1
         END DO
         DO m=1,np
            IF(ind_x(m) /= 0) THEN
               ALLOCATE(Pe_nb(i_n) % tbl (n,m) % exc (3,ind_x(m)))
            END IF
         END DO
         ind_x=0
         DO m=1,mp
            ox=Nei(n) % c(1,m)
            oy=Nei(n) % c(2,m)
            oz=Nei(n) % c(3,m)
            mpe=ind_y(m)
            ind_x(mpe)=ind_x(mpe)+1
            Pe_nb(i_n) % tbl (n,mpe) % exc (:,ind_x(mpe)) = (/ ox, oy, oz /)
         END DO
      END DO
      count0=0
      DO n=1,np
         DO m=1,np
            IF(ALLOCATED(Pe_nb(i_n) % tbl (n,m) % exc)) THEN
               count0=count0+1
!!$               IF(i_n == 1) THEN
!!$                  IF((n == 1 .AND. m == 4) .OR. (n == 4 .AND. m == 1)) THEN
!!$                     mpe=SIZE(Pe_nb(i_n) % tbl (n, m) % exc,2)
!!$                     WRITE(kprint,*) 'n= ',n,'m= ',m,mpe
!!$                     DO k=1,mpe
!!$                        WRITE(kprint,'(3i8)') Pe_nb(i_n) % tbl (n, m) % exc(:,k)
!!$                     END DO
!!$                  END IF
!!$               END IF
            END IF
         END DO
      END DO
      WRITE(kprint,*) count0,' parallel communications are estimated'

      count0=0
      DO n=1,np
         DO m=1,np
            IF(ALLOCATED(Pe_nb(i_n) % tbl (n,m) % exc)) THEN
               mpe=SIZE(Pe_nb(i_n) % tbl (n, m) % exc,2)
               count0=count0+mpe
!!$               IF(i_n == 1) THEN
!!$                  IF((n == 1 .AND. m == 4) .OR. (n == 4 .AND. m == 1)) THEN
!!$                     mpe=SIZE(Pe_nb(i_n) % tbl (n, m) % exc,2)
!!$                     WRITE(kprint,*) 'n= ',n,'m= ',m,mpe
!!$                     DO k=1,mpe
!!$                        WRITE(kprint,'(3i8)') Pe_nb(i_n) % tbl (n, m) % exc(:,k)
!!$                     END DO
!!$                  END IF
!!$               END IF
            END IF
         END DO
      END DO
      WRITE(kprint,*) count0,' cells are transfered for cutoff = ',rcut(i_n)

!!$      DO m=1,np
!!$         IF(ALLOCATED(Pe_nb(i_n) % tbl (m,48) % exc)) THEN
!!$            WRITE(kprint,*) ' m = ',m
!!$            WRITE(kprint,*) (Pe_nb(i_n) % tbl (m,48) % exc(:,mpe),mpe=1,SIZE(Pe_nb(i_n) % tbl (m,48) %&
!!$                 & exc,2))
!!$         END IF
!!$      END DO

    END SUBROUTINE Make_Comm
  END SUBROUTINE PI__Decomposition_nb
  SUBROUTINE PI__GetParameters(nprocsa,npxa,npya,npza)
    INTEGER :: nprocsa,npxa,npya,npza
    nprocsa=PI_nprocs; npxa=PI_npx; npya=PI_npy; npza=PI_npz
  END SUBROUTINE PI__GetParameters
END MODULE PI_Decompose
