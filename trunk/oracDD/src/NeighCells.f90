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
MODULE NeighCells
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Jan 25 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*
  
  USE Constants
  USE Neighbors
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f
  USE Node
  USE Cell
  IMPLICIT none
  PRIVATE
  PUBLIC NeighCells_, NeighCells__Map, NeighCells__MapLarge, NeighCells__Neigh&
       &, Ind_Large, Ind_Small, Nei

  TYPE :: NeighCells__Map
     INTEGER :: nx,ny,nz
  END TYPE NeighCells__Map

  TYPE :: NeighCells__MapLarge
     TYPE(NeighCells__Map), ALLOCATABLE :: pt(:)
  END TYPE NeighCells__MapLarge
  TYPE :: NeighCells__Neigh
     INTEGER, ALLOCATABLE :: c(:,:)
  END TYPE NeighCells__Neigh

!!$
!!$--- List of small lattice units contained in each larger unit
!!$

  TYPE(NeighCells__MapLarge), ALLOCATABLE, SAVE :: Ind_Large(:,:,:)

!!$
!!$--- For each small lattice unit gives the corresponding large
!!$--- lattice unit it belongs to
!!$

  TYPE(NeighCells__Map), ALLOCATABLE, SAVE :: Ind_Small(:,:,:)

!!$
!!$--- For each large lattice unit give the list of neighboring small units
!!$

  TYPE(NeighCells__Neigh), ALLOCATABLE, SAVE :: Nei(:)
  INTEGER, SAVE :: ncx,ncy,ncz
CONTAINS
!!$
!!$--- Constructor
!!$
  FUNCTION NeighCells_(rcut,nprocs,npx,npy,npz) RESULT(out)
    INTEGER :: nprocs,npx,npy,npz
    LOGICAL :: out
    INTEGER :: nx,ny,nz,mx,my,mz,i,j,k,n,np,count0&
         &,count1,o,iv,jv,kv
    REAL(8) :: rcut
    INTEGER, POINTER :: Ind_x(:,:,:)
    REAL(8) :: dx,dy,dz,ddx,ddy,ddz,x1,y1,z1
    CHARACTER(len=max_char) :: labs0
    LOGICAL, POINTER :: mask(:,:,:)=>NULL()
    INTEGER :: vec0(3)
    INTEGER, POINTER :: vec(:)=>NULL()
    LOGICAL, SAVE :: first_time=.TRUE.
    out=.TRUE.
    IF(First_Time) THEN
       nx=NINT(a/rcut)
       ny=NINT(b/rcut)
       nz=NINT(c/rcut)
       IF(nx < npx .OR. ny < npy .OR. nz < npz) THEN
          IF(nx < npx) THEN
             WRITE(labs0,'(i3,i3)') nx,npx
             errmsg_f='Number of neighbor list cells along X is smaller than&
                  & the number of processors: There are '//labs0(1:3)//' neighbor &
                  &list cells and '//labs0(4:6)//' No. processors on X'
          ELSE IF(ny < npy) THEN
             WRITE(labs0,'(i3,i3)') ny,npy
             errmsg_f='Number of neighbor list cells along Y is smaller than&
                  & the number of processors: There are '//labs0(1:3)//' neighbor &
                  &list cells and '//labs0(4:6)//' No. processors on Y'
          ELSE IF(nz < npz) THEN
             WRITE(labs0,'(i3,i3)') nz,npz
             errmsg_f='Number of neighbor list cells along Z is smaller than&
                  & the number of processors: There are '//labs0(1:3)//' neighbor &
                  &list cells and '//labs0(4:6)//' No. processors on Z'
          END IF
          CALL Add_Errors(-1,errmsg_f)
          out=.FALSE.
          RETURN
       END IF
       
       IF(DBLE(MOD(nx,npx))/DBLE(npx) <= 0.5D0) THEN
          ncx=nx-MOD(nx,npx)
       ELSE
          ncx=nx+MOD(npx-MOD(nx,npx),npx)
       END IF
       
       IF(DBLE(MOD(ny,npy))/DBLE(npy) <= 0.5D0) THEN
          ncy=ny-MOD(ny,npy)
       ELSE
          ncy=ny+MOD(npy-MOD(ny,npy),npy)
       END IF
       
       IF(DBLE(MOD(nz,npz))/DBLE(npz) <= 0.5D0) THEN
          ncz=nz-MOD(nz,npz)
       ELSE
          ncz=nz+MOD(npz-MOD(nz,npz),npz)
       END IF
       ALLOCATE(Ind_Large(npx,npy,npz))
       ALLOCATE(Ind_X(npx,npy,npz))
       ALLOCATE(Ind_Small(ncx,ncy,ncz))
       dx=2.d0/ncx
       dy=2.d0/ncy
       dz=2.d0/ncz
       ddx=2.d0/npx
       ddy=2.d0/npy
       ddz=2.d0/npz
       Ind_X=0
       DO i=1,ncx
          x1=(i-1)*dx/ddx
          DO j=1,ncy
             y1=(j-1)*dy/ddy
             DO k=1,ncz
                z1=(k-1)*dz/ddz
                mx=INT(x1)+(SIGN(1.D0,x1-INT(x1))-1.)/2
                my=INT(y1)+(SIGN(1.D0,y1-INT(y1))-1.)/2
                mz=INT(z1)+(sign(1.d0,z1-int(z1))-1.)/2
                mx=MOD(MOD(mx,npx)+npx,npx)+1
                my=MOD(MOD(my,npy)+npy,npy)+1
                mz=MOD(MOD(mz,npz)+npz,npz)+1
                Ind_Small(i,j,k) % nx = mx
                Ind_Small(i,j,k) % ny = my
                Ind_Small(i,j,k) % nz = mz
                Ind_X(mx,my,mz) = Ind_X(mx,my,mz) + 1
             END DO
          END DO
       END DO
       DO i=1,npx
          DO j=1,npy
             DO k=1,npz
                ALLOCATE(Ind_Large(i,j,k) % pt(Ind_X(i,j,k)))
             END DO
          END DO
       END DO
       Ind_X=0
       DO i=1,ncx
          x1=(i-1)*dx/ddx
          DO j=1,ncy
             y1=(j-1)*dy/ddy
             DO k=1,ncz
                z1=(k-1)*dz/ddz
                mx=INT(x1)+(SIGN(1.0D0,x1-INT(x1))-1.0D0)/2
                my=INT(y1)+(SIGN(1.0D0,y1-INT(y1))-1.0D0)/2
                mz=INT(z1)+(sign(1.0D0,z1-int(z1))-1.0D0)/2
                mx=MOD(MOD(mx,npx)+npx,npx)+1
                my=MOD(MOD(my,npy)+npy,npy)+1
                mz=MOD(MOD(mz,npz)+npz,npz)+1
                Ind_X(mx,my,mz) = Ind_X(mx,my,mz) +1
                Ind_Large(mx,my,mz) % pt (Ind_X(mx,my,mz)) % nx = i
                Ind_Large(mx,my,mz) % pt (Ind_X(mx,my,mz)) % ny = j
                Ind_Large(mx,my,mz) % pt (Ind_X(mx,my,mz)) % nz = k
             END DO
          END DO
       END DO
       First_Time=.FALSE.
    END IF
    IF(.NOT. Neighbors_(rcut-1.5D0,ncx,ncy,ncz)) CALL Print_Errors()
       

    IF(ALLOCATED(Nei)) DEALLOCATE(Nei)
    ALLOCATE(Nei(nprocs)); ALLOCATE(mask(ncx,ncy,ncz))
    DO mx=1,npx
       DO my=1,npy
          DO mz=1,npz
             np=SIZE(Ind_Large(mx,my,mz) % pt)
             count0=(mx-1)*npy*npz+(my-1)*npz+mz
             mask=.TRUE.
             DO n=1,np
                i=Ind_Large(mx,my,mz) % pt (n) % nx 
                j=Ind_Large(mx,my,mz) % pt (n) % ny 
                k=Ind_Large(mx,my,mz) % pt (n) % nz
                mask(i,j,k)=.FALSE.
             END DO
             IF(.NOT. Node_()) STOP
             DO n=1,np
                i=Ind_Large(mx,my,mz) % pt (n) % nx-1
                j=Ind_Large(mx,my,mz) % pt (n) % ny-1
                k=Ind_Large(mx,my,mz) % pt (n) % nz-1
                DO o=1,SIZE(Ind_xyz)
                   iv=Ind_xyz(o) % i
                   jv=Ind_xyz(o) % j
                   kv=Ind_xyz(o) % k
                   nx=mod(mod(i+iv,ncx)+ncx,ncx)+1
                   ny=mod(mod(j+jv,ncy)+ncy,ncy)+1
                   nz=mod(mod(k+kv,ncz)+ncz,ncz)+1
                   IF(.NOT. mask(nx,ny,nz)) CYCLE
                   mask(nx,ny,nz)=.FALSE.
                   vec0(1)=nx; vec0(2)=ny; vec0(3)=nz
                   CALL Node__Push(vec0)
                END DO
             END DO
             count1=Node__Size()
             ALLOCATE(Nei(count0) % c (3,count1))
             count1=0
             DO WHILE(Node__Pop(vec))
                count1=count1+1
                Nei(count0) % c(:,count1) = vec
             END DO
          END DO
       END DO
    END DO

  END FUNCTION NeighCells_

END MODULE NeighCells
