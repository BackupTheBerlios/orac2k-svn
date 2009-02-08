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

  SUBROUTINE PI__ShiftIntra(i_p,init,Ind)
    INTEGER :: init,i_p
    INTEGER, OPTIONAL :: Ind

    INTEGER, SAVE :: ShiftTime,source, dest
    INTEGER :: ox,oy,oz,numcell,mpe,mp,m,n
    INTEGER :: nmin,i,j,k,np,AtSt,AtEn,l,q
    INTEGER :: iv(3),Axis
    INTEGER :: ng
    startime=MPI_WTIME()
    npx=PI_npx
    npy=PI_npy
    npz=PI_npz

    CALL Thick_Intra(i_p,MyCutoff(i_p))

!!$
!!$ --- Atoms to send
!!$
    

    iv(1)=npx
    iv(2)=npy
    iv(3)=npz

    IF(PI_Nprocs == 1) RETURN

    CALL Setup

    SELECT CASE(init)
    CASE(_INIT_EXCHANGE_)
       CALL Shift_it(iShift_init)
    CASE(_EXCHANGE_ONLY_)
       CALL Shift_it(Buff_Shift)
    END SELECT


!!$
!!$--- Add to Calls counter: One data exchange has occurred
!!$
    CALL PI__Add_Calls
    endtime=MPI_WTIME()
    CALL PI__Time_It(startime,endtime) 
  CONTAINS
    SUBROUTINE Shift_it(Routine)
      INTERFACE
         SUBROUTINE Routine(i_p,Axis,Dir,scnd_half)
           INTEGER, OPTIONAL :: scnd_half
           INTEGER :: Axis,Dir,i_p
         END SUBROUTINE Routine
      END INTERFACE
      IF(PRESENT(Ind)) THEN
         IF(iv(1) == 1 .AND. iv(2) /= 1) THEN
            CALL Routine(i_p,2,_PLUS_)
            CALL Routine(i_p,3,_PLUS_)
            CALL Routine(i_p,3,_MINUS_,_2NDHALF_)
         ELSE IF(iv(1) == 1 .AND. iv(2) == 1) THEN
            CALL Routine(i_p,3,_PLUS_)
         ELSE
            CALL Routine(i_p,1,_MINUS_)
            CALL Routine(i_p,2,_PLUS_)
            CALL Routine(i_p,3,_PLUS_)
            CALL Routine(i_p,3,_MINUS_,_2NDHALF_)
            CALL Routine(i_p,2,_MINUS_,_2NDHALF_)
         END IF
      ELSE
         IF(iv(1) == 1 .AND. iv(2) /= 1) THEN
            CALL Routine(i_p,2,_PLUS_)
            CALL Routine(i_p,3,_PLUS_)
            CALL Routine(i_p,2,_MINUS_)
            CALL Routine(i_p,3,_MINUS_)
         ELSE IF(iv(1) == 1 .AND. iv(2) == 1) THEN
            CALL Routine(i_p,3,_PLUS_)
            CALL Routine(i_p,3,_MINUS_)
         ELSE
            CALL Routine(i_p,1,_MINUS_)
            CALL Routine(i_p,2,_PLUS_)
            CALL Routine(i_p,3,_PLUS_)
            CALL Routine(i_p,1,_PLUS_)
            CALL Routine(i_p,2,_MINUS_)
            CALL Routine(i_p,3,_MINUS_)
         END IF
      END IF
    END SUBROUTINE Shift_it
  END SUBROUTINE PI__ShiftIntra

  SUBROUTINE iShift_init(i_p,Axis,Dir,scnd_half)
    INTEGER, OPTIONAL :: scnd_half
    INTEGER :: Axis,Dir,i_p
    INTEGER :: nn,n,m,l,count0,mx,my,mz,numcell,ox,oy,oz,mpe,mp&
         &,nmin,i,j,k,MyCell,count1,nind_f,nx,ny,nz
    INTEGER :: NoAtm_s,NoAtm_r,AtSt,AtEn,NoGrp_s3,NoGrp_r3,q,grp_no&
         &,np,nind_o,NoGrp_s,NoGrp_r,nmax
    INTEGER :: source,dest
    REAL(8) :: x,y,z,qq(4),out,xc,yc,zc,xa,ya,za,xd,yd,zd
    REAL(8) :: v1(3),v0,v2(3),rsq,aux1,aux2
    REAL(8) :: point(3)
    REAL(8) :: vc(3),tx,ty,tz
    INTEGER, POINTER :: iBuff_s(:),iBuff_r(:)
    REAL(8), ALLOCATABLE :: Buff_s(:,:),Buff_r(:,:)
    INTEGER, ALLOCATABLE, SAVE :: ind_g(:),ind_o(:)
    INTEGER, ALLOCATABLE :: Struct_s(:),Struct_r(:)
    INTEGER :: nstruct_s,nstruct_r
    LOGICAL :: oks
    REAL(8) :: Margin(3),Margin2_1,Margin2_3,Margin2_2
    REAL(8) :: Margin1(3),Margin2(3),Xmin,Xmax,Ymin,Ymax,Zmin,zmax
    LOGICAL :: ok_X,ok_Y,ok_Z
    INTEGER, SAVE :: MyCalls=0
    REAL(8) :: Axis_L,Axis_R,X_L,X_R,tmass,xmass,xpga,ypga,zpga,xpg,ypg,zpg
    CHARACTER(len=max_char) :: lab0
    INTEGER :: iv_s(3),iv_r(3),g0

    Calls=Calls+1
    IF(Calls > SIZE(iShift)) THEN
       WRITE(lab0,'(i1)') Calls
       errmsg_f='Exchanges are set to 6 per box, but they were '&
            &//TRIM(lab0)
       CALL Add_Errors(-1,errmsg_f)
       CALL Print_Errors()
    END IF
    
    IF(MyCalls == 0) THEN
       ALLOCATE(ind_g(SIZE(Groupa)))
       ALLOCATE(ind_o(SIZE(Atoms)))
       ALLOCATE(MyMaps(SIZE(Groupa)))
    END IF

    oks=PRESENT(scnd_half)
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



    Axis_L=Margin(Axis)-Dir*rcut(Axis)*0.6D0
    Axis_R=Margin(Axis)
    X_L=Margin2_1
    X_R=Margin2_1+rcut(1)*0.6D0

    IF(i_p == 1 .AND. Calls == 1) THEN
       CALL Delete_Maps
       CALL Copy_Maps(Map_n0)
    ELSE IF(i_p == 2 .AND. Calls == 1) THEN
       CALL Delete_Maps
       CALL Copy_Maps(Map_n1)
    END IF

    count1=0
    nmax=nMyMaps
    DO nn=1,nmax
       n=MyMaps(nn) % Grp_No
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
             count1=count1+1
             ind_g(count1)=nn
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
          count1=count1+1
          ind_g(count1)=nn
       END IF
    END DO

    NoGrp_s=count1
    Nstruct_s=0

    If(NoGrp_s /= 0) Nstruct_s=SUM(MyMaps(Ind_g(1:count1)) % g0)+2*count1
    
    Allocate(Struct_s(Nstruct_s))

    count0=0
    count1=0
    DO nn=1,NoGrp_s
       n=Ind_g(nn)
       Struct_s(count1+1)=MyMaps(n) % Grp_no
       Struct_s(count1+2)=MyMaps(n) % g0
       count1=count1+2

       DO m=1,MyMaps(n) % g0
          Struct_s(count1+m)=MyMaps(n) % idx(m)
          ind_o(count0+m)=MyMaps(n) % idx(m)
       END DO
       count1=count1+MyMaps(n) % g0
       count0=count0+MyMaps(n) % g0
    END DO
    NoAtm_s=count0
    nstruct_r=0
    NoAtm_r=0
    NoGrp_r=0

    iShift(Calls) % sh(i_p) % NoAtm_s=NoAtm_s
    iShift(Calls) % sh(i_p) % NoAtm_r=0

    
    iv_s(1)=NoGrp_s;iv_s(2)=NoAtm_s;iv_s(3)=nstruct_s
    CALL MPI_SENDRECV(iv_s,3,MPI_INTEGER4,dest,0,iv_r&
         &,3,MPI_INTEGER4,source,0,PI_Comm_Cart,STATUS,ierr)

    NoGrp_r=iv_r(1);NoAtm_r=iv_r(2);Nstruct_r=iv_r(3)
!!$
!!$--- Set ishifts
!!$


    iShift(Calls) % sh(i_p) % NoAtm_r=NoAtm_r
    IF(ALLOCATED(iShift(Calls) % sh(i_p) % iBuff_S))&
         & DEALLOCATE(iShift(Calls) % sh(i_p) % iBuff_S)
    IF(ALLOCATED(iShift(Calls) % sh(i_p) % iBuff_R))&
         & DEALLOCATE(iShift(Calls) % sh(i_p) % iBuff_R)

    ALLOCATE(iShift(Calls) % sh(i_p) % iBuff_S(NoAtm_S))
    ALLOCATE(iShift(Calls) % sh(i_p) % iBuff_R(NoAtm_R))

    iShift(Calls) % sh(i_p) % iBuff_S(1:NoAtm_S)=ind_o(1:NoAtm_S)


    ALLOCATE(Struct_r(Nstruct_r))
    ALLOCATE(Buff_s(3,NoGrp_s))
    ALLOCATE(Buff_r(3,NoGrp_r))
    n=0
    m=0
    DO WHILE(n < Nstruct_s)
       m=m+1
       l=Struct_s(n+1)
       Buff_s(1,m)=Groupa(l) % xa
       Buff_s(2,m)=Groupa(l) % ya
       Buff_s(3,m)=Groupa(l) % za
       g0=Struct_s(n+2)
       n=n+2+g0
    END DO

    NoGrp_s3=NoGrp_s*3
    NoGrp_r3=NoGrp_r*3

    CALL MPI_SENDRECV(Struct_s,Nstruct_s,MPI_INTEGER4,dest,2,Struct_r&
         &,Nstruct_r,MPI_INTEGER4,source,2,PI_Comm_Cart,STATUS,ierr)

    CALL MPI_SENDRECV(Buff_s,NoGrp_s3,MPI_REAL8,dest,4,Buff_r&
         &,NoGrp_r3,MPI_REAL8,source,4,PI_Comm_Cart,STATUS,ierr)

    n=0
    m=0
    count1=0
    Do While(n < Nstruct_r)
       m=m+1
       l =Struct_r(n+1)
       g0=Struct_r(n+2)
       Groupa(l) % xa=Buff_r(1,m)
       Groupa(l) % ya=Buff_r(2,m)
       Groupa(l) % za=Buff_r(3,m)
       groupa(l) % Knwn = 2

       nMyMaps=nMyMaps+1
       MyMaps(nMyMaps) % Grp_no=l
       MyMaps(nMyMaps) % g0=g0
       If(g0 /= 0) Allocate(MyMaps(nMyMaps) % idx(g0))
       Do q=1,g0
          MyMaps(nMyMaps) % idx(q)=Struct_r(n+2+q)
          iShift(Calls) % sh(i_p) % iBuff_r(count1+q)=Struct_r(n+2+q)
          Atoms(Struct_r(n+2+q)) % knwn = 2
       End Do
       n=n+2+g0
       count1=count1+g0
    End Do
  Contains
    SUBROUTINE Copy_Maps(Map_n)
      TYPE(MyMAps_) :: Map_n(:)
      INTEGER :: n,g0,m
      
      nMyMaps=SIZE(Map_n)
      DO n=1,nMyMaps
         MyMaps(n) % Grp_no=Map_n(n) % Grp_no
         g0=Map_n(n) % g0
         MyMaps(n) % g0=g0
         If(g0 /= 0) Then
            Allocate(MyMaps(n) % idx(g0))
            MyMaps(n)%idx=Map_n(n) % idx
         End If
      END DO
    END SUBROUTINE Copy_Maps
    SUBROUTINE Delete_Maps
      INTEGER :: n,g0
      DO n=1,nMyMaps
         IF(ALLOCATED(MyMaps(n) % idx)) DEALLOCATE(MyMaps(n) % idx)
         MyMaps(n) % g0=0
         MyMaps(n) % Grp_No=0
      END DO
      nMyMaps=0
    END SUBROUTINE Delete_Maps
 END SUBROUTINE IShift_init
  SUBROUTINE Buff_Shift(i_p,Axis,Dir,scnd_half)
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
    INTEGER, POINTER :: iBuff_s(:),iBuff_r(:)
    REAL(8), ALLOCATABLE :: Buff_s(:,:),Buff_r(:,:)
    INTEGER, ALLOCATABLE, SAVE :: ind_o(:)
    LOGICAL :: oks
    REAL(8) :: Margin(3),Margin2_1,Margin2_3,Margin2_2
    REAL(8) :: Margin1(3),Margin2(3),Xmin,Xmax,Ymin,Ymax,Zmin,zmax
    LOGICAL :: ok_X,ok_Y,ok_Z
    INTEGER, SAVE :: MyCalls=0
    REAL(8) :: Axis_L,Axis_R,X_L,X_R,tmass,xmass,xpga,ypga,zpga,xpg,ypg,zpg
    CHARACTER(len=max_char) :: lab0

    
    Calls=Calls+1
    IF(Calls > SIZE(iShift)) THEN
       WRITE(lab0,'(i1)') Calls
       errmsg_f='Exchanges are set to 6 per box, but they were '&
            &//TRIM(lab0)
       CALL Add_Errors(-1,errmsg_f)
       CALL Print_Errors()
    END IF

    CALL MPI_CART_SHIFT(PI_Comm_Cart,Axis-1,Dir,source,dest,ierr)

    NoAtm_s=iShift(Calls) % sh(i_p) % NoAtm_s
    NoAtm_r=iShift(Calls) % sh(i_p) % NoAtm_r

    iBuff_s=>iShift(Calls) % sh(i_p) % iBuff_s
    iBuff_r=>iShift(Calls) % sh(i_p) % iBuff_r

    ALLOCATE(Buff_s(3,NoAtm_s))
    ALLOCATE(Buff_r(3,NoAtm_r))
    DO n=1,NoAtm_s
       q=iBuff_s(n)
       Buff_s(1,n)=Atoms(q) % x
       Buff_s(2,n)=Atoms(q) % y
       Buff_s(3,n)=Atoms(q) % z
    END DO

    NoAtm_s3=NoAtm_s*3
    NoAtm_r3=NoAtm_r*3

    CALL MPI_SENDRECV(Buff_s,NoAtm_s3,MPI_REAL8,dest,3,Buff_r&
         &,NoAtm_r3,MPI_REAL8,source,3,PI_Comm_Cart,STATUS,ierr)

    DO n=1,NoAtm_r
       q=iBuff_r(n)
       xc=Buff_r(1,n)
       yc=Buff_r(2,n)
       zc=Buff_r(3,n)
       atoms(q) % x=xc
       atoms(q) % y=yc
       atoms(q) % z=zc
    END DO

  END SUBROUTINE Buff_Shift
  SUBROUTINE Setup
    Calls=0
  END SUBROUTINE Setup
