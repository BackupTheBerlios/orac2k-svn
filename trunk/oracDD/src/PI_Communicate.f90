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
#define _2NDHALF_ 1

#ifdef HAVE_MPI
  USE mpi
#endif
  USE Geometry
  USE Forces, ONLY: Force, Radii
  USE UNITS
  USE PI_
  USE Groups
  USE Atom
  USE Cell
  USE IndBox
  
  USE Errors,ONLY: Add_errors=>Add, Print_Errors, errmsg_f
  IMPLICIT none
  PRIVATE
  PUBLIC PI__Shift,PI__Fold_F,PI__Write_Stats, PI__ZeroSecondary,&
       & PI__ZeroPrimary 
  TYPE :: Communicate__Chain
     INTEGER :: i,j,k
     INTEGER :: p
  END TYPE Communicate__Chain

  REAL(8), ALLOCATABLE, SAVE :: Buffer(:)
  REAL(8), SAVE :: dx,dy,dz,ddx,ddy,ddz,vp(3),vd(3)
  INTEGER, SAVE :: ncx,ncy,ncz,npx,npy,npz,nsend,Nrec

  TYPE(Communicate__Chain), ALLOCATABLE,SAVE :: Chain_xyz(:)
  INTEGER, ALLOCATABLE,SAVE :: Head_xyz(:)
  REAL(8), SAVE :: Thick(3),rcut(3)
  LOGICAL, SAVE :: ok_pme
  REAL(8), PARAMETER :: one=1.0D0,two=2.0D0,half=0.5D0
  INTEGER, SAVE :: Calls=0
  TYPE :: Statistics
     REAL(8) :: KByte_S=0.0_8,KByte_R=0.0_8
     REAL(8) :: Atoms_S=0.0_8,Atoms_R=0.0_8
  END type Statistics
  TYPE(Statistics), SAVE :: Comms
  TYPE :: Cuts
     REAL(8) :: r(3)=0.0_8
  END type Cuts
  TYPE(Cuts), SAVE :: rcuts0(10)
CONTAINS
  SUBROUTINE Thickness(i_p)
    INTEGER :: i_p
    REAL(8) :: x1,x2,y1,y2,z1,z2,v1(3),v2(3),v3(3),r1(3),r2(3)&
           &,r3(3),r0(3),qq(4)
    
    IF(rcuts0(i_p) % r(1) /= 0.0_8) THEN
       rcut(:)=rcuts0(i_p) % r(:)
       RETURN
    END IF

    ddx=2.d0/npx
    ddy=2.d0/npy
    ddz=2.d0/npz

    r0=0.0D0
    v1(1)=1.0D0
    v1(2)=0.0D0
    v1(3)=0.0D0
    r1=Convert(v1)
    v2(1)=0.0D0
    v2(2)=1.0D0
    v2(3)=0.0D0
    r2=Convert(v2)
    v3(1)=0.0D0
    v3(2)=0.0D0
    v3(3)=1.0D0
    r3=Convert(v3)

!!$    
!!$---- Thickness along x
!!$
    qq=Equation_Plane(r0,r2,r3)
    Thick(1)=1.0D0/ABS(qq(1)*r1(1)+qq(2)*r1(2)+qq(3)*r1(3)-qq(4))
!!$    
!!$---- Thickness along y
!!$
    qq=Equation_Plane(r0,r1,r3)
    Thick(2)=1.0D0/ABS(qq(1)*r2(1)+qq(2)*r2(2)+qq(3)*r2(3)-qq(4))
!!$    
!!$---- Thickness along z
!!$
    qq=Equation_Plane(r0,r1,r2)
    Thick(3)=1.0D0/ABS(qq(1)*r3(1)+qq(2)*r3(2)+qq(3)*r3(3)-qq(4))
    
!!$    
!!$---- Thickness along x
!!$
    rcuts0(i_p) % r (:)=(Radii(i_p) % out+Radii(i_p)% update)*Thick(:)
    rcut(:)=rcuts0(i_p) % r (:)
  CONTAINS
    FUNCTION Convert(v1) RESULT(out)
      REAL(8) :: v1(3),out(3)
      REAL(8) :: r1(3)
      r1(1)=co(1,1)*v1(1)+co(1,2)*v1(2)+co(1,3)*v1(3)
      r1(2)=co(2,1)*v1(1)+co(2,2)*v1(2)+co(2,3)*v1(3)
      r1(3)=co(3,1)*v1(1)+co(3,2)*v1(2)+co(3,3)*v1(3)
      out=r1
    END FUNCTION Convert
  END SUBROUTINE Thickness
  SUBROUTINE PI__Shift(i_p,pme)
    INTEGER, OPTIONAL :: pme
    INTEGER :: i_p
    INTEGER, SAVE :: ShiftTime,source, dest
    INTEGER :: ox,oy,oz,numcell,mpe,mp,m,n
    INTEGER :: nmin,i,j,k,np,AtSt,AtEn,l,q
    INTEGER :: iv(3),Axis
    npx=PI_npx
    npy=PI_npy
    npz=PI_npz
    CALL Thickness(i_p)
    
    ok_pme=PRESENT(pme)

!!$
!!$ --- Atoms to send
!!$
    
    Nrec=0
    Nsend=0

    iv(1)=npx
    iv(2)=npy
    iv(3)=npz
!!$
    IF(PI_Nprocs == 1) RETURN

    IF(ok_pme) THEN
       IF(iv(1) == 1 .AND. iv(2) /= 1) THEN
          CALL Shift(i_p,2,_PLUS_)
          CALL Shift(i_p,3,_PLUS_)
          CALL Shift(i_p,2,_MINUS_)
          CALL Shift(i_p,3,_MINUS_)
       ELSE IF(iv(1) == 1 .AND. iv(2) == 1) THEN
          CALL Shift(i_p,3,_PLUS_)
          CALL Shift(i_p,3,_MINUS_)
       ELSE
          CALL Shift(i_p,1,_MINUS_)
          CALL Shift(i_p,2,_PLUS_)
          CALL Shift(i_p,3,_PLUS_)
          CALL Shift(i_p,1,_PLUS_)
          CALL Shift(i_p,2,_MINUS_)
          CALL Shift(i_p,3,_MINUS_)
       END IF
    ELSE
       IF(iv(1) == 1 .AND. iv(2) /= 1) THEN
          CALL Shift(i_p,2,_PLUS_)
          CALL Shift(i_p,3,_PLUS_)
          CALL Shift(i_p,3,_MINUS_,_2NDHALF_)
       ELSE IF(iv(1) == 1 .AND. iv(2) == 1) THEN
          CALL Shift(i_p,3,_PLUS_)
       ELSE
          CALL Shift(i_p,1,_MINUS_)
          CALL Shift(i_p,2,_PLUS_)
          CALL Shift(i_p,3,_PLUS_)
          CALL Shift(i_p,3,_MINUS_,_2NDHALF_)
          CALL Shift(i_p,2,_MINUS_,_2NDHALF_)
       END IF
    END IF
    Calls=Calls+1

  END SUBROUTINE PI__Shift
  SUBROUTINE Shift(i_p,Axis,Dir,scnd_half)
    INTEGER, OPTIONAL :: scnd_half
    INTEGER :: Axis,Dir,i_p
    INTEGER :: nn,n,m,l,count0,mx,my,mz,numcell,ox,oy,oz,mpe,mp&
         &,nmin,i,j,k,MyCell,count1,nind_f,nx,ny,nz
    INTEGER :: NoAtm_s,NoAtm_r,AtSt,AtEn,NoAtm_s3,NoAtm_r3,q,grp_no&
         &,np,nind_o
    INTEGER :: source,dest
    REAL(8) :: x,y,z,qq(4),out,xc,yc,zc,xa,ya,za,xd,yd,zd
    REAL(8) :: v1(3),v0,v2(3),rsq,aux1,aux2
    REAL(8) :: point(3)
    REAL(8) :: vc(3),tx,ty,tz
    REAL(8), ALLOCATABLE :: Buff_s(:,:),Buff_r(:,:)
    INTEGER, ALLOCATABLE :: iBuff_s(:),iBuff_r(:)
    INTEGER, ALLOCATABLE, SAVE :: ind_o(:)
    LOGICAL :: oks
    REAL(8) :: Margin(3),Margin2_1,Margin2_3,Margin2_2
    REAL(8) :: Margin1(3),Margin2(3),Xmin,Xmax,Ymin,Ymax,Zmin,zmax
    LOGICAL :: ok_X,ok_Y,ok_Z
    INTEGER, SAVE :: MyCalls=0
    REAL(8) :: timeb,startime,endtime,timea,Axis_L,Axis_R,X_L,X_R
    
    IF(MyCalls == 0) THEN
       ALLOCATE(ind_o(SIZE(Groupa)))
       timea=0.0D0
    END IF

    oks=PRESENT(scnd_half) .AND. (.NOT. ok_PME)
    MyCalls=MyCalls+1
    CALL MPI_CART_SHIFT(PI_Comm_Cart,Axis-1,Dir,source,dest,ierr)
    
    aux1=0.5D0*(1.0D0+DBLE(Dir))

    ox=PI__Ranks(PI_Node_Cart+1) % nx+aux1
    oy=PI__Ranks(PI_Node_Cart+1) % ny+aux1
    oz=PI__Ranks(PI_Node_Cart+1) % nz+aux1

    Margin(1)=DBLE(ox)*ddx
    Margin(2)=DBLE(oy)*ddy
    Margin(3)=DBLE(oz)*ddz

    Margin2_1=DBLE(PI__Ranks(PI_Node_Cart+1) % nx+1.0D0)*ddx
    Margin2_2=DBLE(PI__Ranks(PI_Node_Cart+1) % ny)*ddy
    Margin2_3=DBLE(PI__Ranks(PI_Node_Cart+1) % nz)*ddz


    count0=0
    count1=0

    Axis_L=Margin(Axis)-Dir*rcut(Axis)
    Axis_R=Margin(Axis)
    X_L=Margin2_1
    X_R=Margin2_1+rcut(1)
    DO n=1,SIZE(Groupa)
       IF(Groupa(n) % knwn == 0) CYCLE
       IF(Groupa(n) % knwn == 1 .AND. oks) CYCLE
       IF(oks .AND. Axis == 2) THEN
          v1(1)=Groupa(n) % xa
          v1(1)=v1(1)-Two*ANINT(Half*(v1(1)-1.0D0))

          aux1=v1(1)-X_L
          aux1=aux1-Two*ANINT(Half*aux1)
          aux2=v1(1)-X_R
          aux2=aux2-Two*ANINT(Half*aux2)
          IF(aux1 > 0.0D0 .AND. aux2 < 0.0D0) THEN
             AtSt=Groupa(n) % AtSt
             AtEn=Groupa(n) % AtEn
             count1=count1+1
             ind_o(count1)=n
             count0=count0+(AtEn-AtSt+1)
          END IF
          CYCLE
       END IF
       v1(1)=Groupa(n) % xa
       v1(2)=Groupa(n) % ya
       v1(3)=Groupa(n) % za
       v1(Axis)=v1(Axis)-Two*ANINT(Half*(v1(Axis)-1.0D0))
       aux1=v1(Axis)-Axis_L
       aux1=aux1-Two*ANINT(Half*aux1)
       aux2=v1(Axis)-Axis_R
       aux2=aux2-Two*ANINT(Half*aux2)
       IF(Dir*aux1 > 0.0D0 .AND. Dir*aux2 < 0.0D0 ) THEN
          AtSt=Groupa(n) % AtSt
          AtEn=Groupa(n) % AtEn
          count1=count1+1
          ind_o(count1)=n
          count0=count0+(AtEn-AtSt+1)
       END IF
    END DO
    
    NoAtm_s=count0
    nind_o=count1
    NoAtm_r=0


    CALL MPI_SENDRECV(NoAtm_s,1,MPI_INTEGER4,dest,0,NoAtm_r&
         &,1,MPI_INTEGER4,source,0,PI_Comm_Cart,STATUS,ierr)

    ALLOCATE(Buff_s(3,NoAtm_s))
    ALLOCATE(iBuff_s(NoAtm_s))
    ALLOCATE(Buff_r(3,NoAtm_r))
    ALLOCATE(iBuff_r(NoAtm_r))
        
    count0=0
    DO m=1,nind_o
       l=ind_o(m)
       AtSt=Groupa(l) % AtSt
       AtEn=Groupa(l) % AtEn
       DO q=AtSt,AtEn
          count0=count0+1
          Buff_s(1,count0)=Atoms(q) % x
          Buff_s(2,count0)=Atoms(q) % y
          Buff_s(3,count0)=Atoms(q) % z
          IBuff_s(count0) = q
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

    Comms % Atoms_S=Comms % Atoms_S+DBLE(NoAtm_s)
    Comms % Atoms_R=Comms % Atoms_R+DBLE(NoAtm_r)
    Comms % KByte_S=Comms % KByte_S+DBLE(NoAtm_s*3*8)/1024.0_8
    Comms % KByte_R=Comms % KByte_R+DBLE(NoAtm_r*3*8)/1024.0_8
    
  END SUBROUTINE Shift
!!$
!!$---- Fold forces
!!$
  SUBROUTINE PI__Fold_F(fp,i_p)
    TYPE(Force) :: fp(:)
    INTEGER :: i_p
    INTEGER, SAVE :: ShiftTime,source, dest
    INTEGER :: ox,oy,oz,numcell,mpe,mp,m,n
    INTEGER :: nmin,i,j,k,np,AtSt,AtEn,l,q,nn,grp_no
    TYPE(Force), ALLOCATABLE :: fp0(:)
    INTEGER :: iv(3),Axis
    
    ALLOCATE(fp0(SIZE(Atoms)))
    fp0(:) % x=0.0D0
    fp0(:) % y=0.0D0
    fp0(:) % z=0.0D0
    fp0(IndBox_a_t(:))=fp(:)
    
    CALL Thickness(i_p)
    npx=PI_npx
    npy=PI_npy
    npz=PI_npz
    
    
    
!!$
!!$ --- Fold Forces
!!$

    iv(1)=npx
    iv(2)=npy
    iv(3)=npz

    IF(PI_Nprocs == 1) RETURN

    IF(ok_pme) THEN
       IF(iv(1) == 1 .AND. iv(2) /= 1) THEN
          CALL Fold_F(fp0,i_p,3,_PLUS_)
          CALL Fold_F(fp0,i_p,3,_MINUS_)
          CALL Fold_F(fp0,i_p,2,_PLUS_)
          CALL Fold_F(fp0,i_p,2,_MINUS_)
       ELSE IF(iv(1) == 1 .AND. iv(2) == 1) THEN
          CALL Fold_F(fp0,i_p,3,_PLUS_)
          CALL Fold_F(fp0,i_p,3,_MINUS_)
       ELSE
          CALL Fold_F(fp0,i_p,1,_PLUS_)
          CALL Fold_F(fp0,i_p,1,_MINUS_)
          CALL Fold_F(fp0,i_p,2,_PLUS_)
          CALL Fold_F(fp0,i_p,2,_MINUS_)
          CALL Fold_F(fp0,i_p,3,_MINUS_)
          CALL Fold_F(fp0,i_p,3,_PLUS_)
       END IF
    ELSE
       IF(iv(1) == 1 .AND. iv(2) /= 1) THEN
          CALL Fold_F(fp0,i_p,3,_PLUS_)
          CALL Fold_F(fp0,i_p,3,_MINUS_)
          CALL Fold_F(fp0,i_p,2,_MINUS_)
       ELSE IF(iv(1) == 1 .AND. iv(2) == 1) THEN
          CALL Fold_F(fp0,i_p,3,_MINUS_)
       ELSE
          CALL Fold_F(fp0,i_p,3,_PLUS_)
          CALL Fold_F(fp0,i_p,2,_PLUS_)
          CALL Fold_F(fp0,i_p,3,_MINUS_)
          CALL Fold_F(fp0,i_p,2,_MINUS_)
          CALL Fold_F(fp0,i_p,1,_MINUS_)
       END IF
    END IF
    Calls=Calls+1

    DO nn=1,SIZE(IndBox_a_t)
       n=Indbox_a_t(nn)
       fp(nn)=fp0(n)
       Grp_No=Atoms(n) % Grp_No
       IF(Groupa(Grp_No) % knwn == 3) Groupa(Grp_No) % knwn = 2
    END DO
  END SUBROUTINE PI__Fold_F
  SUBROUTINE Fold_F(fp0,i_p,Axis,Dir)
    TYPE(Force) :: fp0(:)
    INTEGER :: Axis,Dir,i_p
    INTEGER :: nn,n,m,l,count0,mx,my,mz,numcell,ox,oy,oz,mpe,mp&
         &,nmin,i,j,k,MyCell,count1,nind_f,nx,ny,nz
    INTEGER :: NoAtm_s,NoAtm_r,AtSt,AtEn,NoAtm_s3,NoAtm_r3,q,grp_no&
         &,np,nind_o
    INTEGER :: source,dest
    REAL(8) :: x,y,z,qq(4),out,xc,yc,zc,xa,ya,za,xd,yd,zd
    REAL(8) :: v1(3),v0,v2(3),rsq,aux1,aux2
    REAL(8) :: point(3)
    REAL(8) :: vc(3),tx,ty,tz
    REAL(8), ALLOCATABLE :: Buff_s(:,:),Buff_r(:,:)
    INTEGER, ALLOCATABLE :: iBuff_s(:),iBuff_r(:)
    INTEGER, ALLOCATABLE, SAVE :: ind_o(:)
    REAL(8) :: Margin(3),Margin2_1,Margin2_3,Margin2_2
    REAL(8) :: Margin1(3),Margin2(3),Xmin,Xmax,Ymin,Ymax,Zmin,zmax
    LOGICAL :: ok_X,ok_Y,ok_Z
    INTEGER, SAVE :: MyCalls=0
    REAL(8) :: timeb,startime,endtime,timea,Axis_L,Axis_R,X_L,X_R
    
    IF(MyCalls == 0) THEN
       ALLOCATE(ind_o(SIZE(Groupa)))
       timea=0.0D0
    END IF

    MyCalls=MyCalls+1
  
    CALL MPI_CART_SHIFT(PI_Comm_Cart,Axis-1,Dir,source,dest,ierr)
    
    aux1=0.5D0*(1.0D0+DBLE(Dir))

    ox=PI__Ranks(PI_Node_Cart+1) % nx+aux1
    oy=PI__Ranks(PI_Node_Cart+1) % ny+aux1
    oz=PI__Ranks(PI_Node_Cart+1) % nz+aux1

    Margin(1)=DBLE(ox)*ddx
    Margin(2)=DBLE(oy)*ddy
    Margin(3)=DBLE(oz)*ddz


    count0=0
    count1=0


    Axis_L=Margin(Axis)
    Axis_R=Margin(Axis)+Dir*rcut(Axis)
    DO n=1,SIZE(Groupa)
       IF(Groupa(n) % knwn == 2) THEN
          v1(1)=Groupa(n) % xa
          v1(2)=Groupa(n) % ya
          v1(3)=Groupa(n) % za
          v1(Axis)=v1(Axis)-Two*ANINT(Half*(v1(Axis)-1.0D0))
          aux1=v1(Axis)-Axis_L
          aux1=aux1-Two*ANINT(Half*aux1)
!!$          aux2=v1(Axis)-Axis_R
!!$          aux2=aux2-Two*ANINT(Half*aux2)
          IF(Dir*aux1 > 0.0D0) THEN
             AtSt=Groupa(n) % AtSt
             AtEn=Groupa(n) % AtEn
             count1=count1+1
             ind_o(count1)=n
             count0=count0+(AtEn-AtSt+1)
          END IF
       END IF
    END DO
    
    NoAtm_s=count0
    nind_o=count1
    NoAtm_r=0

    CALL MPI_SENDRECV(NoAtm_s,1,MPI_INTEGER4,dest,0,NoAtm_r&
         &,1,MPI_INTEGER4,source,0,PI_Comm_Cart,STATUS,ierr)

    ALLOCATE(Buff_s(3,NoAtm_s))
    ALLOCATE(iBuff_s(NoAtm_s))
    ALLOCATE(Buff_r(3,NoAtm_r))
    ALLOCATE(iBuff_r(NoAtm_r))
        
    count0=0
    DO m=1,nind_o
       l=ind_o(m)
       AtSt=Groupa(l) % AtSt
       AtEn=Groupa(l) % AtEn
       DO q=AtSt,AtEn
          count0=count0+1
          Buff_s(1,count0)=fp0(q) % x
          Buff_s(2,count0)=fp0(q) % y
          Buff_s(3,count0)=fp0(q) % z
          IBuff_s(count0) = q
          Grp_No=Atoms(q) % Grp_No
          groupa(Grp_No) % Knwn = 3
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
    Comms % Atoms_S=Comms % Atoms_S+DBLE(NoAtm_s)
    Comms % Atoms_R=Comms % Atoms_R+DBLE(NoAtm_r)
    Comms % KByte_S=Comms % KByte_S+DBLE(NoAtm_s*3*8)/1024.0_8
    Comms % KByte_R=Comms % KByte_R+DBLE(NoAtm_r*3*8)/1024.0_8    
  END SUBROUTINE Fold_F

  SUBROUTINE PI__ZeroSecondary
    INTEGER :: n,AtSt,AtEn,mm
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
    
  END SUBROUTINE PI__ZeroSecondary
  SUBROUTINE PI__ZeroPrimary
    INTEGER :: n,AtSt,AtEn,mm
    DO n=1,SIZE(Groupa)
       IF(Groupa(n) % knwn /= 2 ) THEN
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
    
  END SUBROUTINE PI__ZeroPrimary
  FUNCTION PI__Write_Stats RESULT(out)
    LOGICAL :: out
    REAL(8) :: Kbyte_s, Kbyte_r, Atoms_s, Atoms_r

    out=.TRUE.
    IF(PI_Nprocs == 1) RETURN
#ifdef HAVE_MPI
    Kbyte_s=Comms % KByte_s/DBLE(Calls)
    Kbyte_r=Comms % KByte_r/DBLE(Calls)
    Atoms_s=Comms % Atoms_s/DBLE(Calls)
    Atoms_r=Comms % Atoms_r/DBLE(Calls)
    CALL MPI_ALLREDUCE(Kbyte_s,Kbyte_s,1,MPI_REAL8,MPI_SUM,PI_Comm_Cart,ierr)
    CALL MPI_ALLREDUCE(Kbyte_r,Kbyte_r,1,MPI_REAL8,MPI_SUM,PI_Comm_Cart,ierr)
    CALL MPI_ALLREDUCE(Atoms_s,Atoms_s,1,MPI_REAL8,MPI_SUM,PI_Comm_Cart,ierr)
    CALL MPI_ALLREDUCE(Atoms_r,Atoms_r,1,MPI_REAL8,MPI_SUM,PI_Comm_Cart,ierr)
    IF(PI_Node_Cart == 0) THEN
       Kbyte_s=Kbyte_s/DBLE(PI_Nprocs)
       Kbyte_r=Kbyte_r/DBLE(PI_Nprocs)
       Atoms_s=Atoms_s/DBLE(PI_Nprocs)
       Atoms_r=Atoms_r/DBLE(PI_Nprocs)
       WRITE(*,100) Atoms_s,Kbyte_s,Atoms_r,Kbyte_r
    END IF
100 FORMAT(/'=====>    Average data transfer by each CPU per full shift    <====='&
         &/'=====>    ',f12.2,' atoms (',f12.4,' KB ) sent          <====='&
         &/'=====>    ',f12.2,' atoms (',f12.4,' KB ) received      <=&
         &===='/)
#endif
  END FUNCTION PI__Write_Stats
END MODULE PI_Communicate
