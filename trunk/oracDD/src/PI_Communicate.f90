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
  
#define _X_ 1
#define _Y_ 2
#define _Z_ 3
#define  _PLUS_  1
#define _MINUS_ -1

#ifdef HAVE_MPI
  USE mpi
#endif
  USE Forces, ONLY: Force
  USE UNITS
  USE PI_
  USE NeighCells
  USE Groups
  USE Atom
  USE Cell
  USE LittleBoxes
  
  USE Errors,ONLY: Add_errors=>Add, Print_Errors, errmsg_f
  USE Neighbors_S, ONLY: Neighbors_S__Particles,Neighbors_S__nc, nc,&
       & Neighbors_S__Chain, Chain_xyz, Head_xyz
  IMPLICIT none
  PRIVATE
  PUBLIC PI__Shift,PI__Fold_F
  REAL(8), ALLOCATABLE, SAVE :: Buffer(:)
  REAL(8), SAVE :: dx,dy,dz,ddx,ddy,ddz,vp(3),vd(3)
  INTEGER, SAVE :: ncx,ncy,ncz,npx,npy,npz,nsend,Nrec
  TYPE(NeighCells__Neigh), POINTER, SAVE :: Nei(:)
CONTAINS
  SUBROUTINE PI__Shift(i_p,pme)
    INTEGER, OPTIONAL :: pme
    INTEGER :: i_p
    INTEGER, SAVE :: ShiftTime,source, dest
    INTEGER :: ox,oy,oz,numcell,mpe,mp,m,n
    INTEGER :: nmin,i,j,k,np,AtSt,AtEn,l,q
    LOGICAL, POINTER :: Mask(:)
    INTEGER, pointer :: Nrecs(:),Nsends(:)
    INTEGER :: iv(3),Axis
      
    ncx=nc(i_p) % x
    ncy=nc(i_p) % y
    ncz=nc(i_p) % z
    npx=PI_npx
    npy=PI_npy
    npz=PI_npz
    dx=2.0D0/DBLE(npx)
    dy=2.0D0/DBLE(npy)
    dz=2.0D0/DBLE(npz)
    vd(1)=dx
    vd(2)=dy
    vd(3)=dz
    ddx=2.0D0/DBLE(ncx)
    ddy=2.0D0/DBLE(ncy)
    ddz=2.0D0/DBLE(ncz)
    vp(1)=DBLE(PI__Ranks(PI_Node_Cart+1) % nx)*dx
    vp(2)=DBLE(PI__Ranks(PI_Node_Cart+1) % ny)*dy
    vp(3)=DBLE(PI__Ranks(PI_Node_Cart+1) % nz)*dz

    
    IF(PRESENT(pme)) THEN
       Nei=>Neighc_pme_(i_p) % Nei
    ELSE
       Nei=>Neighc_(i_p) % Nei
    END IF

!!$
!!$ --- Atoms to send
!!$
    
    Nrec=0
    Nsend=0

    iv(1)=npx
    iv(2)=npy
    iv(3)=npz

    DO Axis=1,3
       SELECT CASE(iv(Axis))
       CASE(1)
          CYCLE
       CASE(2)
          CALL Shift(i_p,Axis,_PLUS_)
       CASE DEFAULT
          CALL Shift(i_p,Axis,_PLUS_)
          CALL Shift(i_p,Axis,_MINUS_)
       END SELECT
    END DO
    


!!$    CALL Shift(i_p,_X_,_PLUS_)
!!$    CALL Shift(i_p,_Y_,_PLUS_)
!!$    CALL Shift(i_p,_Z_,_PLUS_)
!!$    CALL Shift(i_p,_X_,_MINUS_)
!!$    CALL Shift(i_p,_Y_,_MINUS_)
!!$    CALL Shift(i_p,_Z_,_MINUS_)
!!$
!!$    ELSE
!!$       SELECT CASE(Direction)
!!$       CASE(1)
!!$          CALL Shift(i_p,_X_,_PLUS_)
!!$          CALL Shift(i_p,_Y_,_PLUS_)
!!$          CALL Shift(i_p,_Z_,_PLUS_)
!!$       CASE(-1)
!!$          CALL Shift(i_p,_X_,_MINUS_)
!!$          CALL Shift(i_p,_Y_,_MINUS_)
!!$          CALL Shift(i_p,_Z_,_MINUS_)
!!$       END SELECT
!!$    END IF

!!$    DO n=1,SIZE(Atoms)
!!$       IF(groupa(Atoms(n) % Grp_No) % knwn == 1) THEN
!!$          Atoms(n) % x =0.0D0
!!$          Atoms(n) % y =0.0D0
!!$          Atoms(n) % z =0.0D0
!!$          Atoms(n) % xa =0.0D0
!!$          Atoms(n) % ya =0.0D0
!!$          Atoms(n) % za =0.0D0
!!$       END IF
!!$    END DO
  END SUBROUTINE PI__Shift
  SUBROUTINE Shift(i_p,Axis,Dir)
    INTEGER :: Axis,Dir,i_p
    INTEGER :: nn,n,m,l,count0,mx,my,mz,numcell,ox,oy,oz,mpe,mp&
         &,nmin,i,j,k
    INTEGER :: NoAtm_s,NoAtm_r,AtSt,AtEn,NoAtm_s3,NoAtm_r3,q,grp_no&
         &,np,startime,endtime,timea,nind_o
    INTEGER :: source,dest
    REAL(8) :: vc(3)
    REAL(8), POINTER :: Buff_s(:,:),Buff_r(:,:)
    INTEGER, POINTER :: iBuff_s(:),iBuff_r(:)
    INTEGER, POINTER :: ind_o(:)
    INTEGER, POINTER :: NoAtmss(:)

    CALL MPI_CART_SHIFT(PI_Comm_Cart,Axis-1,Dir,source,dest,ierr)

    IF(.NOT. Neighbors_S__Particles(i_p)) CALL Print_Errors()

    
    numcell=PI__Ranks(dest+1) % n
    mp=SIZE(Nei(numcell) % c,2)
    
    ALLOCATE(ind_o(mp))

    count0=0
    DO m=1,mp
       ox=Nei(numcell) % c(1,m)-1
       oy=Nei(numcell) % c(2,m)-1
       oz=Nei(numcell) % c(3,m)-1
       count0=count0+1
       mpe=oz+ncz*(ox*ncy+oy)+1
       ind_o(count0)=mpe
    END DO
    nind_o=count0
    
    count0=0
    DO m=1,nind_o
       mpe=ind_o(m)
       l=Head_xyz (mpe)
       nmin=0
       DO WHILE(l > nmin)
          AtSt=Groupa(l) % AtSt
          AtEn=Groupa(l) % AtEn
          count0=count0+(AtEn-AtSt+1)
          l=Chain_xyz(l) % p
       END DO
    END DO
    NoAtm_s=count0

    CALL MPI_SENDRECV(NoAtm_s,1,MPI_INTEGER4,dest,0,NoAtm_r&
         &,1,MPI_INTEGER4,source,0,PI_Comm_Cart,STATUS,ierr)

    ALLOCATE(Buff_s(3,NoAtm_s))
    ALLOCATE(iBuff_s(NoAtm_s))
    ALLOCATE(Buff_r(3,NoAtm_r))
    ALLOCATE(iBuff_r(NoAtm_r))
        
    count0=0
    DO m=1,nind_o
       mpe=ind_o(m)
       l=Head_xyz (mpe)
       nmin=0
       DO WHILE(l > nmin)
          AtSt=Groupa(l) % AtSt
          AtEn=Groupa(l) % AtEn
          DO q=AtSt,AtEn
             count0=count0+1
             Buff_s(1,count0)=Atoms(q) % x
             Buff_s(2,count0)=Atoms(q) % y
             Buff_s(3,count0)=Atoms(q) % z
             IBuff_s(count0) = q
          END DO
          l=Chain_xyz(l) % p
       END DO
       
    END DO
    NoAtm_s3=NoAtm_s*3
    NoAtm_r3=NoAtm_r*3
    CALL MPI_SENDRECV(iBuff_s,NoAtm_s,MPI_INTEGER4,dest,1,iBuff_r&
         &,NoAtm_r,MPI_INTEGER4,source,1,PI_Comm_Cart,STATUS,ierr)
    CALL MPI_SENDRECV(Buff_s,NoAtm_s3,MPI_REAL8,dest,2,Buff_r&
         &,NoAtm_r3,MPI_REAL8,source,2,PI_Comm_Cart,STATUS,ierr)
    
    DO nn=1,NoAtm_r
       n=iBuff_r(nn)
       atoms(n) % x=Buff_r(1,nn)
       atoms(n) % y=Buff_r(2,nn)
       atoms(n) % z=Buff_r(3,nn)
       Atoms(n) % xa = oc(1,1)*Atoms(n) % x+oc(1,2)*Atoms(n) % y+oc(1,3)*Atoms(n) % z    
       Atoms(n) % ya = oc(2,1)*Atoms(n) % x+oc(2,2)*Atoms(n) % y+oc(2,3)*Atoms(n) % z    
       Atoms(n) % za = oc(3,1)*Atoms(n) % x+oc(3,2)*Atoms(n) % y+oc(3,3)*Atoms(n) % z    
       Grp_No=Atoms(n) % Grp_No
       groupa(Grp_No) % Knwn = 2
    END DO
    IF(.NOT. Groups__Update_Knwn()) CALL Print_Errors()
  END SUBROUTINE Shift
!!$
!!$---- Fold forces
!!$
  SUBROUTINE PI__Fold_F(fp,i_p)
    TYPE(Force) :: fp(:)
    INTEGER :: i_p
    INTEGER, SAVE :: ShiftTime,source, dest
    INTEGER :: ox,oy,oz,numcell,mpe,mp,m,n
    INTEGER :: nmin,i,j,k,np,AtSt,AtEn,l,q,nn
    LOGICAL, POINTER :: Mask(:)
    TYPE(Force), POINTER :: fp0(:)
    INTEGER :: iv(3),Axis
    
    ALLOCATE(fp0(SIZE(Atoms)))
    fp0(:) % x=0.0D0
    fp0(:) % y=0.0D0
    fp0(:) % z=0.0D0
    fp0(IndBox_t(:))=fp(:)
    
    ncx=nc(i_p) % x
    ncy=nc(i_p) % y
    ncz=nc(i_p) % z
    npx=PI_npx
    npy=PI_npy
    npz=PI_npz
    
    
    Nei=>Neighc_pme_(i_p) % Nei
    
!!$
!!$ --- Fold Forces
!!$
    IF(.NOT. Neighbors_S__Particles(i_p)) CALL Print_Errors()

    iv(1)=npx
    iv(2)=npy
    iv(3)=npz

    DO Axis=1,3
       SELECT CASE(iv(Axis))
       CASE(1)
          CYCLE
       CASE(2)
          CALL Fold_F(i_p,Axis,_PLUS_)
       CASE DEFAULT
          CALL Fold_F(i_p,Axis,_PLUS_)
          CALL Fold_F(i_p,Axis,_MINUS_)
       END SELECT
    END DO

    DO nn=1,SIZE(IndBox_t)
       n=Indbox_t(nn)
       fp(nn)=fp0(n)
    END DO
  CONTAINS
    
    SUBROUTINE Fold_F(i_p,Axis,Dir)
      INTEGER :: Axis,Dir,i_p
      INTEGER :: nn,n,m,l,count0,mx,my,mz,numcell,ox,oy,oz,mpe,mp&
           &,nmin,i,j,k
      INTEGER :: NoAtm_s,NoAtm_r,AtSt,AtEn,NoAtm_s3,NoAtm_r3,q,grp_no&
           &,np,startime,endtime,timea,nind_o
      INTEGER :: source,dest
      REAL(8) :: vc(3),aux
      REAL(8), POINTER :: Buff_s(:,:),Buff_r(:,:)
      INTEGER, POINTER :: iBuff_s(:),iBuff_r(:)
      INTEGER, POINTER :: ind_o(:)
      LOGICAL, POINTER :: Mask(:,:,:)
      
      
      
      CALL MPI_CART_SHIFT(PI_Comm_Cart,Axis-1,Dir,source,dest,ierr)

      ALLOCATE(Mask(ncx,ncy,ncz)) ; mask=.FALSE.

      numcell=PI__Ranks(dest+1) % n
      mx=PI__Ranks(dest+1) % nx+1
      my=PI__Ranks(dest+1) % ny+1
      mz=PI__Ranks(dest+1) % nz+1
      mp=SIZE(Ind_Large(mx,my,mz) % pt)

      DO n=1,mp
         ox=Ind_Large(mx,my,mz) % pt (n) % nx
         oy=Ind_Large(mx,my,mz) % pt (n) % ny
         oz=Ind_Large(mx,my,mz) % pt (n) % nz
         Mask(ox,oy,oz)=.TRUE.
      END DO

      mp=SIZE(Nei(numcell) % c,2)
      DO n=1,mp
         ox=Nei(numcell) % c(1,n)
         oy=Nei(numcell) % c(2,n)
         oz=Nei(numcell) % c(3,n)
         Mask(ox,oy,oz)=.TRUE.
      END DO

      numcell=PI__Ranks(PI_Node_Cart+1) % n
      mp=SIZE(Nei(numcell) % c,2)
      
      ALLOCATE(ind_o(mp))
      
      count0=0
      DO m=1,mp
         ox=Nei(numcell) % c(1,m)
         oy=Nei(numcell) % c(2,m)
         oz=Nei(numcell) % c(3,m)
         IF(Mask(ox,oy,oz)) THEN
            count0=count0+1
            mpe=(ox-1)*ncy*ncz+(oy-1)*ncz+oz
            ind_o(count0)=mpe
         END IF
      END DO
      nind_o=count0

      count0=0
      DO m=1,nind_o
         mpe=ind_o(m)
         l=Head_xyz (mpe)
         nmin=0
         DO WHILE(l > nmin)
            AtSt=Groupa(l) % AtSt
            AtEn=Groupa(l) % AtEn
            count0=count0+(AtEn-AtSt+1)
            l=Chain_xyz(l) % p
         END DO
      END DO
      NoAtm_s=count0
      CALL MPI_SENDRECV(NoAtm_s,1,MPI_INTEGER4,dest,0,NoAtm_r&
           &,1,MPI_INTEGER4,source,0,PI_Comm_Cart,STATUS,ierr)
      
      ALLOCATE(Buff_s(3,NoAtm_s))
      ALLOCATE(iBuff_s(NoAtm_s))
      ALLOCATE(Buff_r(3,NoAtm_r))
      ALLOCATE(iBuff_r(NoAtm_r))
      
      count0=0
      DO m=1,nind_o
         mpe=ind_o(m)
         l=Head_xyz (mpe)
         nmin=0
         DO WHILE(l > nmin)
            AtSt=Groupa(l) % AtSt
            AtEn=Groupa(l) % AtEn
            DO q=AtSt,AtEn
               count0=count0+1
               Buff_s(1,count0)=fp0(q) % x
               Buff_s(2,count0)=fp0(q) % y
               Buff_s(3,count0)=fp0(q) % z
               IBuff_s(count0) = q
            END DO
            l=Chain_xyz(l) % p
         END DO
         
      END DO
      NoAtm_s3=NoAtm_s*3
      NoAtm_r3=NoAtm_r*3
      CALL MPI_SENDRECV(iBuff_s,NoAtm_s,MPI_INTEGER4,dest,1,iBuff_r&
           &,NoAtm_r,MPI_INTEGER4,source,1,PI_Comm_Cart,STATUS,ierr)
      
      CALL MPI_SENDRECV(Buff_s,NoAtm_s3,MPI_REAL8,dest,2,Buff_r&
           &,NoAtm_r3,MPI_REAL8,source,2,PI_Comm_Cart,STATUS,ierr)

      DO nn=1,NoAtm_r
         n=iBuff_r(nn)
         fp0(n) % x=fp0(n) % x+Buff_r(1,nn)
         fp0(n) % y=fp0(n) % y+Buff_r(2,nn)
         fp0(n) % z=fp0(n) % z+Buff_r(3,nn)
      END DO

    END SUBROUTINE Fold_F
  END SUBROUTINE PI__Fold_F
END MODULE PI_Communicate
