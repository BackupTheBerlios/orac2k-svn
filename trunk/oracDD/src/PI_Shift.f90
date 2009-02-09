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
Module PI_Shift
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Nov  6 2008 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*

#ifdef HAVE_MPI
  Use mpi
#endif
  Use BoxGeometry
  Use Pi_
  Use Forces, Only: Radii
  Use Pi_cutoffs, Only: Thickness, ddx,ddy,ddz, rcut
  Use Pi_statistics, Only: Pi__write_stats=>Write_It, Pi__time_it&
       &=>Time_It, Pi__sample_exchange=>Sample_Exchange, Pi__add_calls&
       &=>Add_Calls
  Use Groups
  Use Atom
  Use Cell
  Use Constants, Only: max_pars,max_data,max_char
  Use Errors, Only: Add_Errors=>Add, Print_Errors, error_args, errmsg_f
  Implicit none
  Private
  Public Setup,iShift_init, Buff_Shift, iShift, Indx, iBuffer
  Integer, Parameter :: Nshell_max=3
  Type :: Indx
     Integer :: NoAtm_s,NoAtm_r
     Integer :: NoGrp_s,NoGrp_r
     Integer, Allocatable :: ibuff_r(:)
     Integer, Allocatable :: ibuff_s(:)
  End type Indx
  Type :: iBuffer
     Type(Indx) :: sh(Nshell_max)
  End type iBuffer
  Type(iBuffer), Save, Target :: iShift(6)
  Integer, Save :: Calls=0
  Real(8), Parameter :: one=1.0D0,two=2.0D0,half=0.5D0
  Real(8), Save :: startime,endtime,startime0,endtime0
  Logical, Save :: ok_pme
Contains
  Subroutine Setup(ok_pmea)
    Logical :: ok_pmea
    Calls=0
    ok_pme=ok_pmea
  End Subroutine Setup

  Subroutine iShift_init(i_p,Axis,Dir,scnd_half)
    Integer, Optional :: scnd_half
    Integer :: Axis,Dir,i_p
    Integer :: nn,n,m,l,count0,mx,my,mz,numcell,ox,oy,oz,mpe,mp&
         &,nmin,i,j,k,MyCell,count1,nind_f,nx,ny,nz
    Integer :: NoAtm_s,NoAtm_r,AtSt,AtEn,NoGrp_s3,NoGrp_r3,q,grp_no&
         &,np,nind_o,NoGrp_s,NoGrp_r,nmax
    Integer :: source,dest
    Real(8) :: x,y,z,qq(4),out,xc,yc,zc,xa,ya,za,xd,yd,zd
    Real(8) :: v1(3),v0,v2(3),rsq,aux1,aux2,aux3
    Real(8) :: point(3)
    Real(8) :: vc(3),tx,ty,tz
    Integer, Pointer :: iBuff_s(:),iBuff_r(:)
    Real(8), Allocatable :: Buff_s(:,:),Buff_r(:,:)
    Integer, Allocatable, Save :: ind_o(:)
    Logical :: oks
    Real(8) :: Margin(3),Margin2_1,Margin2_3,Margin2_2
    Real(8) :: Margin1(3),Margin2(3),Xmin,Xmax,Ymin,Ymax,Zmin,zmax
    Logical :: ok_X,ok_Y,ok_Z,Change_Nmax
    Integer, Save :: MyCalls=0
    Real(8) :: Axis_L,Axis_R,X_L,X_R,tmass,xmass,xpga,ypga,zpga,xpg,ypg,zpg
    Character(len=max_char) :: lab0
    Integer :: iv_s(2),iv_r(2),iv1(3),ActualFace
    Real(8) :: MyCut

    MyCut=Radii(i_p) % out+Radii(i_p)% update

    Calls=Calls+1
    If(Calls > Size(iShift)) Then
       Write(lab0,'(i1)') Calls
       errmsg_f='Exchanges are set to 6 per box, but they were '&
            &//Trim(lab0)
       Call Add_Errors(-1,errmsg_f)
       Call Print_Errors()
    End If
    
    If(MyCalls == 0) Then
       Allocate(ind_o(Size(Groupa)))
    End If

    oks=Present(scnd_half) .And. (.Not. ok_Pme)
    MyCalls=MyCalls+1

    Call Mpi_cart_shift(Pi_comm_cart,Axis-1,Dir,source,dest,ierr)

    If(Distance_(v1(1),v1(2),v1(3),MyCut,n,Dir,Axis,dest)) Call &
         & Print_Errors()
    
    aux1=0.5D0*(1.0D0+Dble(Dir))

    iv1=(/Pi__ranks(Pi_node_cart+1) % nx,Pi__ranks(Pi_node_cart+1) &
         &% ny,Pi__ranks(Pi_node_cart+1) % nz/)
    iv1(Axis)=iv1(Axis)+aux1

    ox=iv1(1) ; oy=iv1(2) ; oz=iv1(3) 

    Margin(1)=Dble(ox)*ddx
    Margin(2)=Dble(oy)*ddy
    Margin(3)=Dble(oz)*ddz

    Margin2_1=Dble(Pi__ranks(Pi_node_cart+1) % nx+1.0D0)*ddx
    Margin2_2=Dble(Pi__ranks(Pi_node_cart+1) % ny)*ddy
    Margin2_3=Dble(Pi__ranks(Pi_node_cart+1) % nz)*ddz


    count0=0
    count1=0

    Axis_L=Margin(Axis)-Dir*rcut(Axis)
    Axis_R=Margin(Axis)
    X_L=Margin2_1
    X_R=Margin2_1+rcut(1)

    nmax=Size(Groupa)
    Do n=1,nmax
       If(Groupa(n) % knwn == 0) Cycle
       If(Groupa(n) % knwn == 1 .And. oks) Cycle
       If(oks .And. Axis == 2) Then
          v1(1)=Groupa(n) % xa
          v1(1)=v1(1)-Two*Anint(Half*(v1(1)-1.0D0))
          
          aux1=v1(1)-X_L
          aux1=aux1-Two*Anint(Half*aux1)
          aux2=v1(1)-X_R
          aux2=aux2-Two*Anint(Half*aux2)
          If(aux1 > 0.0D0 .And. aux2 < 0.0D0) Then
             If(Distance_(Groupa(n) %xa, Groupa(n) %ya, Groupa(n) %za&
                  &,MyCut,n)) Then 
                AtSt=Groupa(n) % AtSt
                AtEn=Groupa(n) % AtEn
                count1=count1+1
                ind_o(count1)=n
                count0=count0+(AtEn-AtSt+1)
             End If
          End If
          Cycle
       End If
       v1(1)=Groupa(n) % xa
       v1(2)=Groupa(n) % ya
       v1(3)=Groupa(n) % za
       v1(Axis)=v1(Axis)-Two*Anint(Half*(v1(Axis)-1.0D0))
       aux1=v1(Axis)-Axis_L
       aux1=aux1-Two*Anint(Half*aux1)
       aux2=v1(Axis)-Axis_R
       aux2=aux2-Two*Anint(Half*aux2)
       If(Dir*aux1 > 0.0D0 .And. Dir*aux2 < 0.0D0 ) Then
          If(Distance_(Groupa(n) %xa, Groupa(n) %ya, Groupa(n) %za,MyCut,n)) Then
             AtSt=Groupa(n) % AtSt
             AtEn=Groupa(n) % AtEn
             count1=count1+1
             ind_o(count1)=n
             count0=count0+(AtEn-AtSt+1)
          End If
       End If
    End Do
    NoAtm_s=count0
    nind_o=count1
    NoGrp_s=count1
    NoAtm_r=0
    NoGrp_r=0
    iShift(Calls) % sh(i_p) % NoGrp_s=NoGrp_s
    iShift(Calls) % sh(i_p) % NoAtm_s=NoAtm_s
    iShift(Calls) % sh(i_p) % NoGrp_r=0
    iShift(Calls) % sh(i_p) % NoAtm_r=0

    iv_s(1)=NoGrp_s;iv_s(2)=NoAtm_s

    Call Mpi_sendrecv(iv_s,2,Mpi_integer4,dest,0,iv_r&
         &,2,Mpi_integer4,source,0,Pi_comm_cart,Status,ierr)

    NoGrp_r=iv_r(1);NoAtm_r=iv_r(2)

    iShift(Calls) % sh(i_p) % NoGrp_r=NoGrp_r
    iShift(Calls) % sh(i_p) % NoAtm_r=NoAtm_r

    If(Allocated(iShift(Calls) % sh(i_p) % iBuff_S))&
         & Deallocate(iShift(Calls) % sh(i_p) % iBuff_S)
    If(Allocated(iShift(Calls) % sh(i_p) % iBuff_R))&
         & Deallocate(iShift(Calls) % sh(i_p) % iBuff_R)

    Allocate(iShift(Calls) % sh(i_p) % iBuff_S(NoGrp_S))
    Allocate(iShift(Calls) % sh(i_p) % iBuff_R(NoGrp_R))

    iShift(Calls) % sh(i_p) % iBuff_S(1:NoGrp_S)=ind_o(1:NoGrp_S)



    iBuff_s=>iShift(Calls) % sh(i_p) % iBuff_S
    iBuff_r=>iShift(Calls) % sh(i_p) % iBuff_R


    Allocate(Buff_s(3,NoGrp_s))
    Allocate(Buff_r(3,NoGrp_r))

    Do m=1,NoGrp_s
       l=iBuff_s(m)
       Buff_s(1,m)=Groupa(l) % xa
       Buff_s(2,m)=Groupa(l) % ya
       Buff_s(3,m)=Groupa(l) % za
    End Do


    NoGrp_s3=NoGrp_s*3
    NoGrp_r3=NoGrp_r*3

    Call Mpi_sendrecv(iBuff_s,NoGrp_s,Mpi_integer4,dest,2,iBuff_r&
         &,NoGrp_r,Mpi_integer4,source,2,Pi_comm_cart,Status,ierr)
    Call Mpi_sendrecv(Buff_s,NoGrp_s3,Mpi_real8,dest,4,Buff_r&
         &,NoGrp_r3,Mpi_real8,source,4,Pi_comm_cart,Status,ierr)

    Do m=1,NoGrp_r
       l=iBuff_r(m)
       Groupa(l) % xa = Buff_r(1,m)
       Groupa(l) % ya = Buff_r(2,m)
       Groupa(l) % za = Buff_r(3,m)
       groupa(l) % Knwn = 2
       AtSt=Groupa(l) % AtSt
       AtEn=Groupa(l) % AtEn
       Do n=AtSt,AtEn
          Atoms(n) % knwn = 2
       End Do
    End Do
 End Subroutine Ishift_init
  Subroutine Buff_Shift(i_p,Axis,Dir,scnd_half)
    Integer, Optional :: scnd_half
    Integer :: Axis,Dir,i_p
    Integer :: nn,n,m,l,count0,mx,my,mz,numcell,ox,oy,oz,mpe,mp&
         &,nmin,i,j,k,MyCell,count1,nind_f,nx,ny,nz
    Integer :: NoAtm_s,NoAtm_r,AtSt,AtEn,NoAtm_s3,NoAtm_r3,q,grp_no&
         &,np,nind_o,NoGrp_s,NoGrp_r
    Integer :: source,dest
    Real(8) :: x,y,z,qq(4),out,xc,yc,zc,xa,ya,za,xd,yd,zd
    Real(8) :: v1(3),v0,v2(3),rsq,aux1,aux2
    Real(8) :: point(3)
    Real(8) :: vc(3),tx,ty,tz
    Integer, Pointer :: iBuff_s(:),iBuff_r(:)
    Real(8), Allocatable :: Buff_s(:,:),Buff_r(:,:)
    Integer, Allocatable, Save :: ind_o(:)
    Logical :: oks
    Real(8) :: Margin(3),Margin2_1,Margin2_3,Margin2_2
    Real(8) :: Margin1(3),Margin2(3),Xmin,Xmax,Ymin,Ymax,Zmin,zmax
    Logical :: ok_X,ok_Y,ok_Z
    Integer, Save :: MyCalls=0
    Real(8) :: Axis_L,Axis_R,X_L,X_R,tmass,xmass,xpga,ypga,zpga,xpg,ypg,zpg
    Character(len=max_char) :: lab0
    Logical, Allocatable :: Mask(:)
    Integer :: Myreq_s,MyReq_r

    Calls=Calls+1
    If(Calls > Size(iShift)) Then
       Write(lab0,'(i1)') Calls
       errmsg_f='Exchanges are set to 6 per box, but they were '&
            &//Trim(lab0)
       Call Add_Errors(-1,errmsg_f)
       Call Print_Errors()
    End If
    If(MyCalls == 0) Then
       Allocate(Mask(Size(Groupa)))
       Mask=.False.
    End If
    
    Call Mpi_cart_shift(Pi_comm_cart,Axis-1,Dir,source,dest,ierr)

    NoGrp_s=iShift(Calls) % sh(i_p) % NoGrp_s
    NoGrp_r=iShift(Calls) % sh(i_p) % NoGrp_r
    NoAtm_s=iShift(Calls) % sh(i_p) % NoAtm_s
    NoAtm_r=iShift(Calls) % sh(i_p) % NoAtm_r
    Call jBuff(iShift(Calls) % sh,i_p,iBuff_r,iBuff_s,NoGrp_s,NoGrp_r&
         &,NoAtm_S,NoAtm_r) 

!!$    iBuff_s=>iShift(Calls) % sh(i_p) % iBuff_s
!!$    iBuff_r=>iShift(Calls) % sh(i_p) % iBuff_r

    Allocate(Buff_s(3,NoAtm_s))
    Allocate(Buff_r(3,NoAtm_r))
    count0=0
    Do m=1,NoGrp_s
       l=iBuff_s(m)
       AtSt=Groupa(l) % AtSt
       AtEn=Groupa(l) % AtEn
       Do q=AtSt,AtEn
          count0=count0+1
          Buff_s(1,count0)=Atoms(q) % x
          Buff_s(2,count0)=Atoms(q) % y
          Buff_s(3,count0)=Atoms(q) % z
       End Do
    End Do


    NoAtm_s3=NoAtm_s*3
    NoAtm_r3=NoAtm_r*3

    ierr=PointToPoint(buff_s,NoAtm_s3,buff_r,NoAtm_r3,source,dest,Pi_comm_cart)

    nn=0
    Do m=1,NoGrp_r
       l=iBuff_r(m)
       AtSt=Groupa(l) % AtSt
       AtEn=Groupa(l) % AtEn
       xpga=0.0D0
       ypga=0.0D0
       zpga=0.0D0
       xpg=0.0D0
       ypg=0.0D0
       zpg=0.0D0
       Do n=AtSt,AtEn
          xmass=Atoms(n) % pmass
          nn=nn+1
          xc=Buff_r(1,nn)
          yc=Buff_r(2,nn)
          zc=Buff_r(3,nn)
          atoms(n) % x=xc
          atoms(n) % y=yc
          atoms(n) % z=zc
          Atoms(n) % xa = oc(1,1)*xc+oc(1,2)*yc+oc(1,3)*zc    
          Atoms(n) % ya = oc(2,1)*xc+oc(2,2)*yc+oc(2,3)*zc    
          Atoms(n) % za = oc(3,1)*xc+oc(3,2)*yc+oc(3,3)*zc
          xpga = xpga + xmass*Atoms(n) % xa
          ypga = ypga + xmass*Atoms(n) % ya
          zpga = zpga + xmass*Atoms(n) % za
          xpg = xpg + xmass*Atoms(n) % x
          ypg = ypg + xmass*Atoms(n) % y
          zpg = zpg + xmass*Atoms(n) % z
       End Do
       Groupa(l) % xa = xpga 
       Groupa(l) % ya = ypga
       Groupa(l) % za = zpga
       Groupa(l) % x = xpg
       Groupa(l) % y = ypg
       Groupa(l) % z = zpg
    End Do
    Do m=1,iShift(Calls) % sh(i_p) % NoGrp_r
       l=iShift(Calls) % sh(i_p) % iBuff_r(m)
       groupa(l) % Knwn = 2
       AtSt=Groupa(l) % AtSt
       AtEn=Groupa(l) % AtEn
       Do n=AtSt,AtEn
          Atoms(n) % knwn = 2
       End Do
    End Do

  Contains
    Subroutine jBuff(sh,i_p,iBuff_r,iBuff_s,NoGrp_s,NoGrp_r,NoAtm_S&
         &,NoAtm_r)
      Type(Indx), Target :: sh(:)
      Integer :: i_p,NoGrp_s,NoGrp_r,NoAtm_S,NoAtm_r
      Integer, Pointer :: iBuff_s(:),iBuff_r(:)
      Integer, Pointer :: jBuff_s(:),jBuff_r(:)
      Integer, Pointer :: kBuff_s(:),kBuff_r(:)
      Integer, Allocatable, Save, Target :: MyBuff_s(:),MyBuff_r(:)
      Integer :: i_pm

      
      If(i_p == 1) Then
         iBuff_s=>sh(i_p) % iBuff_s
         iBuff_r=>sh(i_p) % iBuff_r
         NoGrp_s=sh(i_p) % NoGrp_s
         NoGrp_r=sh(i_p) % NoGrp_r
         NoAtm_s=sh(i_p) % NoAtm_s
         NoAtm_r=sh(i_p) % NoAtm_r
         Return
      End If

      i_pm=i_p-1
      If(.Not. Allocated(sh(i_pm) % iBuff_s)) Then
         iBuff_s=>sh(i_p) % iBuff_s
         iBuff_r=>sh(i_p) % iBuff_r
         NoGrp_s=sh(i_p) % NoGrp_s
         NoGrp_r=sh(i_p) % NoGrp_r
         NoAtm_s=sh(i_p) % NoAtm_s
         NoAtm_r=sh(i_p) % NoAtm_r
         Return
      End If
      If(Allocated(MyBuff_s)) Then
         Deallocate(Mybuff_s,Mybuff_r)
      End If
      jBuff_s=>sh(i_p) % iBuff_s
      jBuff_r=>sh(i_p) % iBuff_r

      Mask(jBuff_s(:))=.True.
      i_pm=i_p-1
      Do While (i_pm /=0)
         kBuff_s=>sh(i_pm) % iBuff_s
         Mask(kBuff_s(:))=.False.
         i_pm=i_pm-1
      End Do
      NoGrp_s=Count(Mask)
      Allocate(MyBuff_s(NoGrp_s))
      MyBuff_s=Pack(jBuff_s(:),Mask(jBuff_s(:)))
      iBuff_s=>MyBuff_s


      Mask(jBuff_s(:))=.False. !$-- ReSet to .False. 

      Mask(jBuff_r(:))=.True.
      i_pm=i_p-1
      Do While (i_pm /=0)
         kBuff_r=>sh(i_pm) % iBuff_r
         Mask(kBuff_r(:))=.False.
         i_pm=i_pm-1
      End Do
      NoGrp_r=Count(Mask)
      Allocate(MyBuff_r(NoGrp_r))
      MyBuff_r=Pack(jBuff_r(:),Mask(jBuff_r(:)))
      iBuff_r=>MyBuff_r

      Mask(jBuff_r(:))=.False. !$-- ReSet to .False. 

      NoAtm_s=Sum(Groupa(MyBuff_s(:)) % AtEn-Groupa(MyBuff_s(:)) % AtSt+1)
      NoAtm_r=Sum(Groupa(MyBuff_r(:)) % AtEn-Groupa(MyBuff_r(:)) % AtSt+1)

    End Subroutine jBuff


  End Subroutine Buff_Shift

  Function PointToPoint(buff_s,n_s,buff_r,n_r,source,dest,comm) Result(out)
    Integer :: out,n_s,n_r,source,dest,comm
    Real(8) :: buff_s(:,:),buff_r(:,:)

    Integer :: Myreq_s,Myreq_r

    out=0
#ifdef __NonBlocked
    Call Mpi_isend(buff_s,n_s,Mpi_real8,dest,3,comm,MyReq_s,ierr)
    If(ierr /= 0) out=ierr
    Call Mpi_irecv(buff_r,n_r,Mpi_real8,source,3,comm,MyReq_r,ierr)
    If(ierr /= 0) out=ierr
    Call Mpi_wait(MyReq_r,Status,ierr)
    If(ierr /= 0) out=ierr
    Call Mpi_wait(MyReq_s,Status,ierr)
    If(ierr /= 0) out=ierr
#else
    Call Mpi_sendrecv(Buff_s,n_s,Mpi_real8,dest,3,Buff_r&
         &,n_r,Mpi_real8,source,3,comm,Status,ierr)
    If(ierr /= 0) out=ierr    
#endif
  End Function PointToPoint
End Module Pi_shift
