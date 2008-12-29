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
MODULE PI_Fold
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Nov 13 2008 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*

#ifdef HAVE_MPI
  USE mpi
#endif
  USE PI_
  USE PI_Cutoffs
  USE PI_Statistics, ONLY: PI__Write_Stats=>Write_It, PI__Time_It&
       &=>Time_It, PI__Sample_Exchange=>Sample_Exchange, PI__Add_Calls&
       &=>Add_Calls
  USE Forces, ONLY: Force, Radii
  USE Groups
  USE Atom
  USE Cell
  USE Constants, ONLY: max_pars,max_data,max_char
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, error_args, errmsg_f
  IMPLICIT none
  PRIVATE
  PUBLIC Setup,iFold_init, Buff_Fold
  INTEGER, PARAMETER :: NSHell_Max=3
  TYPE :: Indx
     INTEGER :: NoAtm_s,NoAtm_r
     INTEGER :: NoGrp_s,NoGrp_r
     INTEGER, ALLOCATABLE :: ibuff_r(:)
     INTEGER, ALLOCATABLE :: ibuff_s(:)
  END type Indx
  TYPE :: iBuffer
     TYPE(Indx) :: sh(Nshell_max)
  END type iBuffer
  TYPE(iBuffer), SAVE, TARGET :: iFold(6)
  INTEGER, SAVE :: Calls=0
  REAL(8), PARAMETER :: one=1.0D0,two=2.0D0,half=0.5D0
  REAL(8), SAVE :: startime,endtime,startime0,endtime0
  LOGICAL, SAVE :: ok_pme
CONTAINS
  SUBROUTINE Setup(ok_pmea)
    LOGICAL :: ok_pmea
    Calls=0
    ok_pme=ok_pmea
  END SUBROUTINE Setup
  SUBROUTINE iFold_init(fp0,i_p,Axis,Dir)
    TYPE(Force) :: fp0(:)
    INTEGER :: Axis,Dir,i_p
    INTEGER :: nn,n,m,l,count0,mx,my,mz,numcell,ox,oy,oz,mpe,mp&
         &,nmin,i,j,k,MyCell,count1,nind_f,nx,ny,nz
    INTEGER :: NoAtm_s,NoAtm_r,AtSt,AtEn,NoAtm_s3,NoAtm_r3,q,grp_no&
         &,np,nind_o,NoGrp_s,NoGrp_r
    INTEGER :: source,dest,nmax
    REAL(8) :: x,y,z,qq(4),out,xc,yc,zc,xa,ya,za,xd,yd,zd
    REAL(8) :: v1(3),v0,v2(3),rsq,aux1,aux2
    REAL(8) :: point(3)
    REAL(8) :: vc(3),tx,ty,tz
    INTEGER, POINTER :: iBuff_s(:),iBuff_r(:)
    REAL(8), ALLOCATABLE :: Buff_s(:,:),Buff_r(:,:)
    INTEGER, ALLOCATABLE, SAVE :: ind_o(:)
    REAL(8) :: Margin(3),Margin2_1,Margin2_3,Margin2_2
    REAL(8) :: Margin1(3),Margin2(3),Xmin,Xmax,Ymin,Ymax,Zmin,zmax
    LOGICAL :: ok_X,ok_Y,ok_Z
    INTEGER, SAVE :: MyCalls=0
    REAL(8) :: Axis_L,Axis_R,X_L,X_R
    CHARACTER(len=max_char) :: lab0
    
    Calls=Calls+1
    IF(Calls > SIZE(iFold)) THEN
       WRITE(lab0,'(i1)') Calls
       errmsg_f='Folds are set to 6 per box, but they were '&
            &//TRIM(lab0)
       CALL Add_Errors(-1,errmsg_f)
       CALL Print_Errors()
    END IF

    IF(MyCalls == 0) THEN
       ALLOCATE(ind_o(SIZE(Groupa)))
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
    nmax=SIZE(Groupa)

    DO n=1,nmax
       IF(Groupa(n) % knwn == 2) THEN
          v1(1)=Groupa(n) % xa
          v1(2)=Groupa(n) % ya
          v1(3)=Groupa(n) % za
          v1(Axis)=v1(Axis)-Two*ANINT(Half*(v1(Axis)-1.0D0))
          aux1=v1(Axis)-Axis_L
          aux1=aux1-Two*ANINT(Half*aux1)
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
    NoGrp_s=count1
    NoGrp_r=0

    iFold(Calls) % sh(i_p) % NoGrp_s=NoGrp_s
    iFold(Calls) % sh(i_p) % NoAtm_s=NoAtm_s
    iFold(Calls) % sh(i_p) % NoGrp_r=0
    iFold(Calls) % sh(i_p) % NoAtm_r=0


    startime=MPI_WTIME()
    CALL MPI_SENDRECV(NoGrp_s,1,MPI_INTEGER4,dest,3,NoGrp_r&
         &,1,MPI_INTEGER4,source,3,PI_Comm_Cart,STATUS,ierr)
    CALL MPI_SENDRECV(NoAtm_s,1,MPI_INTEGER4,dest,0,NoAtm_r&
         &,1,MPI_INTEGER4,source,0,PI_Comm_Cart,STATUS,ierr)
    endtime=MPI_WTIME()
    CALL PI__Time_It(startime,endtime)

    iFold(Calls) % sh(i_p) % NoGrp_r=NoGrp_r
    iFold(Calls) % sh(i_p) % NoAtm_r=NoAtm_r

    IF(ALLOCATED(iFold(Calls) % sh(i_p) % iBuff_S))&
         & DEALLOCATE(iFold(Calls) % sh(i_p) % iBuff_S)
    IF(ALLOCATED(iFold(Calls) % sh(i_p) % iBuff_R))&
         & DEALLOCATE(iFold(Calls) % sh(i_p) % iBuff_R)

    ALLOCATE(iFold(Calls) % sh(i_p) % iBuff_S(NoGrp_S))
    ALLOCATE(iFold(Calls) % sh(i_p) % iBuff_R(NoGrp_R))

    iFold(Calls) % sh(i_p) % iBuff_S=ind_o(1:NoGrp_S)
    iBuff_s=>iFold(Calls) % sh(i_p) % iBuff_S
    iBuff_r=>iFold(Calls) % sh(i_p) % iBuff_R
    
    ALLOCATE(Buff_s(3,NoAtm_s))
    ALLOCATE(Buff_r(3,NoAtm_r))
        
    count0=0
    DO m=1,NoGrp_s
       l=iBuff_s(m)
       AtSt=Groupa(l) % AtSt
       AtEn=Groupa(l) % AtEn
       Groupa(l) % knwn=3
       DO q=AtSt,AtEn
          count0=count0+1
          Buff_s(1,count0)=fp0(q) % x
          Buff_s(2,count0)=fp0(q) % y
          Buff_s(3,count0)=fp0(q) % z
       END DO
    END DO
       
    NoAtm_s3=NoAtm_s*3
    NoAtm_r3=NoAtm_r*3

    startime=MPI_WTIME()
    CALL MPI_SENDRECV(iBuff_s,NoGrp_s,MPI_INTEGER4,dest,1,iBuff_r&
         &,NoGrp_r,MPI_INTEGER4,source,1,PI_Comm_Cart,STATUS,ierr)
    CALL MPI_SENDRECV(Buff_s,NoAtm_s3,MPI_REAL8,dest,2,Buff_r&
         &,NoAtm_r3,MPI_REAL8,source,2,PI_Comm_Cart,STATUS,ierr)
    endtime=MPI_WTIME()
    CALL PI__Time_It(startime,endtime)

    nn=0
    DO q=1,NoGrp_r
       l=iBuff_r(q)
       AtSt=Groupa(l) % AtSt
       AtEn=Groupa(l) % AtEn
       DO n=AtSt,AtEn
          nn=nn+1
          fp0(n) % x=fp0(n) % x+Buff_r(1,nn)
          fp0(n) % y=fp0(n) % y+Buff_r(2,nn)
          fp0(n) % z=fp0(n) % z+Buff_r(3,nn)
       END DO
    END DO
    CALL PI__Sample_Exchange(NoAtm_S,NoAtm_R)
  END SUBROUTINE IFold_init
  SUBROUTINE Buff_Fold(fp0,i_p,Axis,Dir)
    TYPE(Force) :: fp0(:)
    INTEGER :: Axis,Dir,i_p
    INTEGER :: nn,n,m,l,count0,mx,my,mz,numcell,ox,oy,oz,mpe,mp&
         &,nmin,i,j,k,MyCell,count1,nind_f,nx,ny,nz
    INTEGER :: NoAtm_s,NoAtm_r,AtSt,AtEn,NoAtm_s3,NoAtm_r3,q,grp_no&
         &,np,nind_o,NoGrp_s,NoGrp_r
    INTEGER :: source,dest
    REAL(8) :: x,y,z,qq(4),out,xc,yc,zc,xa,ya,za,xd,yd,zd
    REAL(8) :: v1(3),v0,v2(3),rsq,aux1,aux2
    REAL(8) :: point(3)
    REAL(8) :: vc(3),tx,ty,tz
    REAL(8), ALLOCATABLE :: Buff_s(:,:),Buff_r(:,:)
    INTEGER, POINTER :: iBuff_s(:),iBuff_r(:)
    REAL(8) :: Margin(3),Margin2_1,Margin2_3,Margin2_2
    REAL(8) :: Margin1(3),Margin2(3),Xmin,Xmax,Ymin,Ymax,Zmin,zmax
    LOGICAL :: ok_X,ok_Y,ok_Z
    INTEGER, SAVE :: MyCalls=0
    REAL(8) :: Axis_L,Axis_R,X_L,X_R
    CHARACTER(len=max_char) :: lab0

    Calls=Calls+1
    IF(Calls > SIZE(iFold)) THEN
       WRITE(lab0,'(i1)') Calls
       errmsg_f='Folds are set to 6 per box, but they were '&
            &//TRIM(lab0)
       CALL Add_Errors(-1,errmsg_f)
       CALL Print_Errors()
    END IF
    
  
    CALL MPI_CART_SHIFT(PI_Comm_Cart,Axis-1,Dir,source,dest,ierr)
    NoGrp_s=iFold(Calls) % sh(i_p) % NoGrp_s
    NoGrp_r=iFold(Calls) % sh(i_p) % NoGrp_r
    NoAtm_s=iFold(Calls) % sh(i_p) % NoAtm_s
    NoAtm_r=iFold(Calls) % sh(i_p) % NoAtm_r

    iBuff_s=>iFold(Calls) % sh(i_p) % iBuff_s
    iBuff_r=>iFold(Calls) % sh(i_p) % iBuff_r

    ALLOCATE(Buff_s(3,NoAtm_s))
    ALLOCATE(Buff_r(3,NoAtm_r))
        
    count0=0
    DO m=1,NoGrp_s
       l=iBuff_s(m)
       AtSt=Groupa(l) % AtSt
       AtEn=Groupa(l) % AtEn
       Groupa(l) % knwn = 3
       DO q=AtSt,AtEn
          count0=count0+1
          Buff_s(1,count0)=fp0(q) % x
          Buff_s(2,count0)=fp0(q) % y
          Buff_s(3,count0)=fp0(q) % z
       END DO
    END DO
       
    NoAtm_s3=NoAtm_s*3
    NoAtm_r3=NoAtm_r*3

    startime=MPI_WTIME()
    CALL MPI_SENDRECV(Buff_s,NoAtm_s3,MPI_REAL8,dest,2,Buff_r&
         &,NoAtm_r3,MPI_REAL8,source,2,PI_Comm_Cart,STATUS,ierr)
    endtime=MPI_WTIME()
    CALL PI__Time_It(startime,endtime)

    nn=0
    DO q=1,NoGrp_r
       l=iBuff_r(q)
       AtSt=Groupa(l) % AtSt
       AtEn=Groupa(l) % AtEn
       DO n=AtSt,AtEn
          nn=nn+1
          fp0(n) % x=fp0(n) % x+Buff_r(1,nn)
          fp0(n) % y=fp0(n) % y+Buff_r(2,nn)
          fp0(n) % z=fp0(n) % z+Buff_r(3,nn)
       END DO
    END DO

    CALL PI__Sample_Exchange(NoAtm_S,NoAtm_R)
  END SUBROUTINE Buff_Fold

END MODULE PI_Fold
