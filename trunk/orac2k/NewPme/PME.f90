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
Module Pme
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Tue Jul 29 2008 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*

#include "PmeParameters.h"

  
  Use Ewald
  Use Print_Defs
  Use PmeRfft3d
#ifdef PARALLEL
  Use Mpi
#endif
  Use Pi_
  Use Errors, Only: Add_Errors=>Add, Print_Errors, errmsg_f, errmsg_w
  Implicit none
  Private
  Public Pme_
!!$
!!$--- Pme Data
!!$
  Real(8), Allocatable, Save :: Cq(:,:,:)
  Real(8), Allocatable :: Cq_s(:,:,:),Cq_r(:,:,:)
  Real(8), Allocatable :: Cq_sb(:,:,:)
  Integer, Save :: nx,ny,nz
  Integer, Save :: nfft1,nfft2,nfft3
  Integer, Save :: ndim_fftw1,ndim_fftw2,ndim_fftw3
  Integer, Save :: nfftw1,nfftw2,nfftw2_start,nfftw2_end,nfftw2_local&
       &,nfftw3,nfftw3_start,nfftw3_end,nfftw3_local 
  Integer, Save :: myfft1,myfft2,myfft3
  Integer, Save :: my1_start,my1_end,my2_start&
       &,my2_end,my3_start,my3_end
  Integer, Save :: kstart,kend,jstart,jend,istart,iend,kwstart,kwend

  Integer, Save :: mx,my,mz,mzz
  Integer, Save :: rk_x,rk_y,rk_z
  Real(8), Allocatable, Save :: bsp_mod1(:),bsp_mod2(:),bsp_mod3(:)
  Real(8), Allocatable :: theta1(:,:),theta2(:,:),theta3(:,:)&
       &,dtheta1(:,:),dtheta2(:,:),dtheta3(:,:)
  Integer, Allocatable :: MyIndBox(:)
!!$
!!$--- System Data
!!$
  Real(8), Allocatable :: fx(:),fy(:),fz(:),chg(:),fr1(:),fr2(:),fr3(:)&
       &,phi(:)
  Real(8), Allocatable, Save :: MyCharges(:)
  Integer, Save :: natom,alltoall_dim,ntot_atom
  Integer, Pointer, Save :: order
  Real(8), Save :: recip(3,3),vir(3,3),eer,energy,eer_i,energy_i
  Real(8), Save :: Tol_q=1.0D-5
  Real(8), Save :: volume,alphal,rkcut
  Real(8), Parameter :: pi=3.1415926535897931D0, twopi=2.0D0*pi
Contains
  Subroutine Pme_First(order_a,alphal_a,rkcut_a,charges_a,natom_a,n1&
       &,n2,n3,nprocs)
    Integer :: order_a,natom_a,n1,n2,n3,nprocs
    Real(8) :: alphal_a,charges_a(*),rkcut_a

    order=order_a
    ntot_atom=natom_a
    alphal=alphal_a
    rkcut=rkcut_a
    Allocate(MyCharges(ntot_atom))
    MyCharges=Charges_a(1:ntot_atom)
    
!!$    
!!$--- Initialize fftw when needed and allocate pointers
!!$
    Call Ewald__Validate(n1,n2,n3,nprocs)

    nfft1=n1;nfft2=n2;nfft3=n3;

    If(.Not. Initialize_()) Call Print_Errors()
    Allocate(Cq_r(ndim_fftw1,ndim_fftw2,ndim_fftw3))
    Allocate(Cq_s(myfft1,myfft2,myfft3))
    Allocate(Cq_sb(myfft1,myfft2,myfft3))
    Call PI__
    Call PI__Setup_Cart
  End Subroutine Pme_First

  Subroutine Pme_(xa,ya,za,co,oc,volumea)
    Real(8) :: xa(*),ya(*),za(*),co(3,3),oc(3,3),volumea
    Integer :: i_pa,i,j,k,n,m,count0,AtSt,AtEn,ia,ja,ka,count1,i0,j0,k0
    Real(8) :: startime,endtime,timea,ts1,te1,ts2,te2

!!$    
!!$--- Copy coordinates and charges to local arrays
!!$

    Do i=1,3
       Do j=1,3
          recip(i,j) = 0.5*oc(j,i)
       End Do
    End Do

    natom=ntot_atom

    If(Allocated(theta1)) Then
       Deallocate(theta1,theta2,theta3)
       Deallocate(dtheta1,dtheta2,dtheta3)
       Deallocate(chg,fr1,fr2,fr3,MyIndBox)
       Deallocate(fx,fy,fz,phi)
    End If

    Call PurgeOutsideAtoms(natom)

    Allocate(theta1(order, natom),theta2(order, natom),theta3(order, natom))
    Allocate(dtheta1(order, natom),dtheta2(order, natom),dtheta3(order, natom))
    Allocate(fx(natom),fy(natom),fz(natom),phi(natom))

    fx=0.0D0; fy=0.0D0; fz=0.0D0


    Call Get_Bsplines

    Call Charges_onGrid(Cq_s)
!!$
!!$--- Change Cpu partition along Z axis
!!$

    Call Transpose_Cart2Fftw
    
    If(.Not. do_rfft3d(1,Cq_r)) Return

    Call ScalarSum_Transposed(Cq_r,nfft1,nfft2,nfft3)

    If(.Not. do_rfft3d(-1,Cq_r)) Return

!!$
!!$--- Back to Cartesian Cpu partition
!!$

    Call Transpose_Fftw2cart

    Call Grad_Sum(Cq_s,Energy)

!!$
!!$    fp(MyIndBox(:)) % x = fp(MyIndBox(:)) % x + fx(:)
!!$    fp(MyIndBox(:)) % y = fp(MyIndBox(:)) % y + fy(:)
!!$    fp(MyIndBox(:)) % z = fp(MyIndBox(:)) % z + fz(:)
!!$
    eer_i=eer
!!$    Call En_coul_rec_(eer)
!!$    
!!$    Energy_i=Energy
  Contains
    Subroutine PurgeOutsideAtoms(natom)
      Integer :: natom
      Logical :: ok_i,ok_j,ok_k
      Integer :: nn
      Integer, Allocatable :: ind(:)
      Real(8), Allocatable :: tr1(:),tr2(:),tr3(:),tchg(:)
      
      Real(8) :: gfr1,gfr2,gfr3
      Real(8) :: w1,w2,w3,x,y,z,v1,v2,v3,chg0
      Integer :: ith1,ith2,ith3

      Allocate(ind(natom),tr1(natom),tr2(natom),tr3(natom),tchg(natom))
      nn=0
      Do m=1,natom
         chg0=MyCharges(m)
         x=xa(m)
         y=ya(m)
         z=za(m)
         w1 = x*recip(1,1)+y*recip(2,1)+z*recip(3,1)
         w2 = x*recip(1,2)+y*recip(2,2)+z*recip(3,2)
         w3 = x*recip(1,3)+y*recip(2,3)+z*recip(3,3)
         gfr1 = nfft1*(w1-Anint(w1)+(0.5D0-Sign(0.5D0,w1-Anint(w1))))
         gfr2 = nfft2*(w2-Anint(w2)+(0.5D0-Sign(0.5D0,w2-Anint(w2))))
         gfr3 = nfft3*(w3-Anint(w3)+(0.5D0-Sign(0.5D0,w3-Anint(w3))))

         ok_k=.False.
         k0 = int(gfr3) - order
         Do ith3 = 1,order
            k0 = k0 + 1
            k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
            If(k < kstart .Or. k > kend) Cycle
            ok_k=.True.
         End Do
         If(.Not. ok_k) Cycle
         ok_j=.False.
         j0 = int(gfr2) - order
         Do ith2 = 1,order
            j0 = j0 + 1
            j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2
            If(j < jstart .Or. j > jend) Cycle
            ok_j=.True.
         End Do
         If(.Not. ok_j) Cycle
         ok_i=.False.
         i0 = int(gfr1) - order
         Do ith1 = 1,order
            i0 = i0 + 1
            i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
            If(i < istart .Or. i > iend) Cycle
            ok_i=.True.
         End Do
         If(ok_i) Then
            nn=nn+1
            ind(nn)=m
            tr1(nn)=gfr1
            tr2(nn)=gfr2
            tr3(nn)=gfr3
            tchg(nn)=chg0
         End If
      End Do
      natom=nn
      Allocate(MyIndBox(natom))
      Allocate(chg(natom),fr1(natom),fr2(natom),fr3(natom))
      MyIndBox=ind(1:natom)
      chg=tchg(1:natom)
      fr1=tr1(1:natom)
      fr2=tr2(1:natom)
      fr3=tr3(1:natom)

    End Subroutine PurgeOutsideAtoms
  End Subroutine Pme_
  Function Initialize_() Result(out)
    Real(8), Pointer :: dummy(:,:,:)
    Logical :: out
    out=.True.


    myfft1=nfft1/Pi_npx
    myfft2=nfft2/Pi_npy
    myfft3=nfft3/Pi_npz
    alltoall_dim=myfft1*myfft2*myfft3/(Pi_npx*Pi_npy)

    rk_x=Pi__ranks_fftw (Pi_node_fftw+1) % nx
    rk_y=Pi__ranks_fftw (Pi_node_fftw+1) % ny
    rk_z=Pi__ranks_fftw (Pi_node_fftw+1) % nz

    my1_start=rk_x*myfft1+1
    my1_end=(rk_x+1)*myfft1
    my2_start=rk_y*myfft2+1
    my2_end=(rk_y+1)*myfft2
    my3_start=rk_z*myfft3+1
    my3_end=(rk_z+1)*myfft3
    istart=my1_start; iend=my1_end
    jstart=my2_start; jend=my2_end
    kstart=my3_start; kend=my3_end

    If(Pi_node_fftw ==0) Write(*,100)

    If(.Not. do_rfft3d(0,dummy,nfft1,nfft2,nfft3,nfftw3_start,nfftw3_local&
         &,nfftw2_start,nfftw2_local,ndim_fftw1,ndim_fftw2&
         &,ndim_fftw3,Pi_comm_fftw)) Return

    nfftw1=nfft1
    nfftw2=nfft2
    nfftw3=ndim_fftw3

    kwstart=nfftw3_start
    kwend  =nfftw3_start+nfftw3_local-1
    Call Bsp_moduli
100 Format(/22x,'Finding optimal parameters for Fftws.'/&
     &     22x,'     This will take a while...'/ /)     
  Contains
    Subroutine Bsp_moduli
      Allocate(bsp_mod1(nfft1),bsp_mod2(nfft2),bsp_mod3(nfft3))
      Call load_bsp_moduli(bsp_mod1,bsp_mod2,bsp_mod3,nfft1,nfft2&
           &,nfft3,order)
    End Subroutine Bsp_moduli
  End Function Initialize_

  Subroutine Transpose_Cart2Fftw
    Integer :: o,p,q,r,s,t,rk2_x,rk2_y,n,mzz

    
    Call Mpi_alltoall(Cq_s,alltoall_dim,Mpi_real8, Cq_sb,alltoall_dim&
         &,Mpi_real8,Pi_comm_z,ierr)


    mzz=myfft3/Pi_nprocs_z
    Do n=1,Pi_nprocs_z
       rk2_x=Pi__ranks_z(n) % nx
       rk2_y=Pi__ranks_z(n) % ny
       Do t=1,mzz
          q=(n-1)*mzz+t
          Do p=1,myfft2
             s=rk2_y*myfft2+p
             Do o=1,myfft1
                r=rk2_x*myfft1+o
                Cq_r(r,s,t)=Cq_sb(o,p,q)
             End Do
          End Do
       End Do
    End Do
  End Subroutine Transpose_Cart2Fftw
  Subroutine Transpose_Fftw2cart
    Integer :: o,p,q,r,s,t,rk2_x,rk2_y,n,mzz

    mzz=myfft3/Pi_nprocs_z
    Do n=1,Pi_nprocs_z
       rk2_x=Pi__ranks_z(n) % nx
       rk2_y=Pi__ranks_z(n) % ny
       Do t=1,mzz
          q=(n-1)*mzz+t
          Do p=1,myfft2
             s=rk2_y*myfft2+p
             Do o=1,myfft1
                r=rk2_x*myfft1+o
                Cq_sb(o,p,q)=Cq_r(r,s,t)
             End Do
          End Do
       End Do
    End Do

    Call Mpi_alltoall(Cq_sb,alltoall_dim,Mpi_real8, Cq_s,alltoall_dim&
         &,Mpi_real8,Pi_comm_z,ierr)

  End Subroutine Transpose_Fftw2cart

#include "PME__Sources.f90"
#include "PME_Routines.f90"
End Module Pme
