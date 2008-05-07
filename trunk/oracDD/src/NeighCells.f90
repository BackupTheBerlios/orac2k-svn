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

  USE PI_
  USE Constants
  USE Neighbors_S
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f
  USE Node
  USE Cell
  IMPLICIT none
  PRIVATE
  PUBLIC NeighCells_, NeighCells__Map, NeighCells__MapLarge, NeighCells__Neigh&
       &, Ind_Large, Ind_Small, Nei, NeighCells__Param

  TYPE :: NeighCells__Map
     INTEGER :: nx,ny,nz
  END TYPE NeighCells__Map

  TYPE :: NeighCells__MapLarge
     TYPE(NeighCells__Map), ALLOCATABLE :: pt(:)
  END TYPE NeighCells__MapLarge
  TYPE :: NeighCells__Neigh
     INTEGER, ALLOCATABLE :: c(:,:)
  END TYPE NeighCells__Neigh
  TYPE :: NodeN
     REAL(8) :: rcut=0.0D0
     INTEGER :: npx=0,npy=0,npz=0,ncx=0,ncy=0,ncz=0
     TYPE(NeighCells__Map), ALLOCATABLE :: Ind_S(:,:,:)
     TYPE(NeighCells__MapLarge), ALLOCATABLE :: Ind_L(:,:,:)
     TYPE(NodeN), POINTER :: next
  END TYPE NodeN
  TYPE(NodeN), POINTER :: Root, Current, new_node

!!$
!!$--- List of small lattice units contained in each larger unit
!!$

  TYPE(NeighCells__MapLarge), POINTER :: Ind_Large(:,:,:)=>NULL()

!!$
!!$--- For each small lattice unit gives the corresponding large
!!$--- lattice unit it belongs to
!!$

  TYPE(NeighCells__Map), POINTER :: Ind_Small(:,:,:)=>NULL()

!!$
!!$--- For each large lattice unit give the list of neighboring small units
!!$

  TYPE(NeighCells__Neigh), ALLOCATABLE, SAVE :: Nei(:)
  INTEGER, SAVE :: ncx,ncy,ncz
CONTAINS
!!$
!!$--- Constructor
!!$
  FUNCTION NeighCells_(rcut0,rcut,nprocs,npax,npay,npaz) RESULT(out)
    INTEGER :: nprocs,npax,npay,npaz
    LOGICAL :: out
    INTEGER :: nx,ny,nz,mx,my,mz,i,j,k,n,np,count1,o,iv,jv,kv
    REAL(8) :: rcut,rcut0
    INTEGER, POINTER :: Ind_x(:,:,:)
    REAL(8) :: dx,dy,dz,ddx,ddy,ddz,x1,y1,z1
    CHARACTER(len=max_char) :: labs0
    LOGICAL, POINTER :: mask(:,:,:)=>NULL()
    INTEGER :: vec0(3)
    INTEGER, POINTER :: vec(:)=>NULL()
    INTEGER, SAVE :: count0=0, iter_max=20
    INTEGER :: npx,npy,npz,i_n,nvalues(3),iter_inst
    LOGICAL :: okkk

    out=.TRUE.

    npx=npax; npy=npay; npz=npaz

!!$
!!$--- If it the first time allocate space for root
!!$

    IF(count0 == 0) THEN
       ALLOCATE(root)
       NULLIFY(root%next)
    ELSE
!!$
!!$--- Search if a cutoff match
!!$
       current=>root
       Ind_Large=>NULL()
       Ind_Small=>NULL()
       count1=0
       current=>current % next
       DO WHILE(ASSOCIATED(current))
          IF(current % rcut == rcut) THEN
             Ind_Large=>current % Ind_L
             Ind_Small=>current % Ind_S
             ncx=current % ncx
             ncy=current % ncy
             ncz=current % ncz
             npx=current % npx
             npy=current % npy
             npz=current % npz
             EXIT
          END IF
          count1=count1+1
          current=>current % next
       END DO
    END IF

!!$
!!$--- If a cutoff does not match, Ind_Large and Ind_Small are not associated
!!$

    IF(.NOT. ASSOCIATED(Ind_Small)) THEN
       ALLOCATE(new_node)
       NULLIFY(new_node % next)
       new_node % rcut = rcut
       okkk=.FALSE.
       nx=NINT(a/rcut0)
       ny=NINT(b/rcut0)
       nz=NINT(c/rcut0)
       iter_inst=0
       DO WHILE(.NOT. okkk)
          iter_inst=iter_inst+1
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
          nvalues=Neighbors_S_Check(i_n,rcut,ncx,ncy,ncz)
          IF(nvalues(1) == 0 .AND. nvalues(2) == 0 .AND. nvalues(3) == 0) THEN
             okkk=.TRUE.
          ELSE
             nx=nx+nvalues(1)
             ny=ny+nvalues(2)
             nz=nz+nvalues(3)
          END IF
!!$
!!$--- Abort after iter_Max unsuccessfull iterations
!!$
          IF(iter_inst > iter_Max) okkk=.TRUE.
       END DO


       ALLOCATE(new_node % Ind_L(npx,npy,npz))
       ALLOCATE(new_node % Ind_S(ncx,ncy,ncz))
       ALLOCATE(Ind_X(npx,npy,npz))
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
                new_node % Ind_S(i,j,k) % nx = mx
                new_node % Ind_S(i,j,k) % ny = my
                new_node % Ind_S(i,j,k) % nz = mz
                Ind_X(mx,my,mz) = Ind_X(mx,my,mz) + 1
             END DO
          END DO
       END DO
       DO i=1,npx
          DO j=1,npy
             DO k=1,npz
                ALLOCATE(new_node % Ind_L(i,j,k) % pt(Ind_X(i,j,k)))
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
                new_node % Ind_L(mx,my,mz) % pt (Ind_X(mx,my,mz)) % nx = i
                new_node % Ind_L(mx,my,mz) % pt (Ind_X(mx,my,mz)) % ny = j
                new_node % Ind_L(mx,my,mz) % pt (Ind_X(mx,my,mz)) % nz = k
             END DO
          END DO
       END DO

       current=>root
       DO WHILE(ASSOCIATED(current % next))
          current=>current % next
       END DO
       count0=count0+1
       current % next => new_node
       current => current % next

       Ind_Large=>current % Ind_L
       Ind_Small=>current % Ind_S
       current % ncx=ncx
       current % ncy=ncy
       current % ncz=ncz
       current % npx=npx
       current % npy=npy
       current % npz=npz
    END IF
       
    IF(.NOT. Neighbors_S_(i_n,rcut,ncx,ncy,ncz)) CALL Print_Errors()
     
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
                DO o=1,SIZE(clst(i_n) % Ind_xyz)
                   iv=clst(i_n) % Ind_xyz(o) % i
                   jv=clst(i_n) % Ind_xyz(o) % j
                   kv=clst(i_n) % Ind_xyz(o) % k
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
  SUBROUTINE NeighCells__Param(np,ncxa,ncya,ncza)
    INTEGER :: np,ncxa,ncya,ncza
    np=ncx*ncy*ncz; ncxa=ncx; ncya=ncy; ncza=ncz
  END SUBROUTINE NeighCells__Param
END MODULE NeighCells
