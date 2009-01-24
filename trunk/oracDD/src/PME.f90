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
MODULE PME
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

#include "config.h"


  USE Energies, ONLY: EN_Coul_Rec_=>Coul_Rec_
  USE Print_Defs
  USE Units, ONLY: Pi, TwoPi
  USE Rfft3d
  USE Cell, ONLY: oc,co, Volume
  USE Ewald
  USE Forces, ONLY: Force,fp_m, fp_l, fp_h, fp_ew, FORCE_Pick=>Pick
#ifdef HAVE_MPI
  USE mpi
#endif
  USE PI_
  USE Groups
  USE Atom
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f, errmsg_w
  USE IndBox, ONLY: IndBoxP_, IndBoxP_a_t
  USE PI_Communicate
  IMPLICIT none
  PRIVATE
  PUBLIC PME_
!!$
!!$--- PME Data
!!$
  REAL(8), ALLOCATABLE, SAVE :: Cq(:,:,:)
  REAL(8), POINTER :: Cq_s(:,:,:),Cq_r(:,:,:)
  REAL(8), POINTER :: Cq_sb(:,:,:)
  INTEGER, SAVE :: nx,ny,nz
  INTEGER, SAVE :: nfft1,nfft2,nfft3
  INTEGER, SAVE :: ndim_fftw1,ndim_fftw2,ndim_fftw3
  INTEGER, SAVE :: nfftw1,nfftw2,nfftw2_start,nfftw2_end,nfftw2_local&
       &,nfftw3,nfftw3_start,nfftw3_end,nfftw3_local 
  INTEGER, SAVE :: myfft1,myfft2,myfft3
  INTEGER, SAVE :: my1_start,my1_end,my2_start&
       &,my2_end,my3_start,my3_end
  INTEGER, SAVE :: kstart,kend,jstart,jend,istart,iend,kwstart,kwend

  INTEGER, SAVE :: mx,my,mz,mzz
  INTEGER, SAVE :: rk_x,rk_y,rk_z
  REAL(8), ALLOCATABLE, SAVE :: bsp_mod1(:),bsp_mod2(:),bsp_mod3(:)
  REAL(8), POINTER :: theta1(:,:),theta2(:,:),theta3(:,:)&
       &,dtheta1(:,:),dtheta2(:,:),dtheta3(:,:)
!!$
!!$--- System Data
!!$
  REAL(8), POINTER :: fx(:),fy(:),fz(:),chg(:),fr1(:),fr2(:),fr3(:)&
       &,phi(:)
  INTEGER, SAVE :: natom,alltoall_dim,ntot_atom
  INTEGER, POINTER, SAVE :: order
  REAL(8), SAVE :: recip(3,3),vir(3,3),eer,energy,eer_i,energy_i
  
CONTAINS
  SUBROUTINE PME_(i_pa)
    INTEGER :: i_pa,i,j,k,n,m,count0,AtSt,AtEn,ia,ja,ka,count1,i0,j0,k0
    INTEGER, SAVE :: No_Calls=0
    REAL(8) :: startime,endtime,timea,ts1,te1,ts2,te2
    TYPE(Force), POINTER :: fp(:)
    INTEGER :: i_p

!!$    
!!$--- Initialize fftw when needed and allocate pointers
!!$

    i_p=i_pa-2
    order=>Ewald__Param % order


    IF(No_Calls == 0) THEN       
       IF(.NOT. Initialize_()) CALL Print_Errors()
       No_Calls=No_Calls+1
       RETURN
    END IF

    startime=MPI_WTIME()
    ALLOCATE(Cq_r(ndim_fftw1,ndim_fftw2,ndim_fftw3))
    ALLOCATE(Cq_s(myfft1,myfft2,myfft3))
    ALLOCATE(Cq_sb(myfft1,myfft2,myfft3))

!!$    
!!$--- Copy coordinates and charges to local arrays
!!$

    IF(.NOT. IndBoxP_(Groupa(:) % knwn,Groupa(:) % AtSt,Groupa(:) %&
         & AtEn)) CALL Print_Errors() 

    natom=SIZE(IndBoxP_a_t)
    ALLOCATE(theta1(order, natom),theta2(order, natom),theta3(order, natom))
    ALLOCATE(dtheta1(order, natom),dtheta2(order, natom),dtheta3(order, natom))
    ALLOCATE(chg(natom),fr1(natom),fr2(natom),fr3(natom))
    ALLOCATE(fx(natom),fy(natom),fz(natom),phi(natom))

    fx=0.0D0; fy=0.0D0; fz=0.0D0

    DO i=1,3
       DO j=1,3
          recip(i,j) = 0.5*oc(j,i)
       END DO
    END DO

    CALL Fractionals
    CALL Get_Bsplines

    CALL Charges_onGrid(Cq_s)
!!$
!!$--- Change CPU partition along Z axis
!!$

    CALL Transpose_Cart2FFTW
    
    IF(.NOT. do_rfft3d(1,Cq_r)) RETURN

    CALL ScalarSum_Transposed(Cq_r,nfft1,nfft2,nfft3)

    IF(.NOT. do_rfft3d(-1,Cq_r)) RETURN

!!$
!!$--- Back to Cartesian CPU partition
!!$

    CALL Transpose_FFTW2Cart

    CALL Grad_Sum(Cq_s,Energy)

    ntot_atom=SIZE(Atoms)

    fp=>FORCE_Pick(i_pa)

    fp(IndBoxP_a_t(:)) % x = fp(IndBoxP_a_t(:)) % x + fx(:)
    fp(IndBoxP_a_t(:)) % y = fp(IndBoxP_a_t(:)) % y + fy(:)
    fp(IndBoxP_a_t(:)) % z = fp(IndBoxP_a_t(:)) % z + fz(:)

!!$
!!$--- Fold forces contributions to atoms inside the cell
!!$
    
!!$    CALL PI__Fold_F(fp,i_p,_INIT_)
!!$    CALL PI__Fold_F(fp,i_p,_FOLD_)

!!$    IF(PI_Nprocs == 1) THEN
!!$       WRITE(60,'(i7,3e17.9)') (n,fp(n) % x, fp(n) % y, fp(n) % z, n=1,natom)
!!$    ELSE
!!$       DO m=1,natom
!!$          n=IndBoxP_t(m)
!!$          IF(groupa(Atoms(n) % Grp_No) % Knwn == 1) THEN
!!$             WRITE(60+PI_Node_Cart,'(i7,3e17.9)') n,fp(m) % x, fp(m) % y&
!!$                  &, fp(m) % z
!!$          END IF
!!$       END DO
!!$    END IF

    eer_i=eer
    CALL EN_Coul_Rec_(eer)
    Energy_i=Energy
    endtime=MPI_WTIME()
    timea=endtime-startime
    WRITE(kprint,*) 'timeo ',timea
  CONTAINS
    SUBROUTINE Fractionals
      INTEGER :: n,m,mm,count0,AtSt,AtEn
      REAL(8) :: w1,w2,w3,x,y,z,v1,v2,v3
      
      DO n=1,natom
         m=IndBoxP_a_t(n)
         chg(n)=Atoms(m) % chg
         x=Atoms(m) % x
         y=Atoms(m) % y
         z=Atoms(m) % z
         w1 = x*recip(1,1)+y*recip(2,1)+z*recip(3,1)
         w2 = x*recip(1,2)+y*recip(2,2)+z*recip(3,2)
         w3 = x*recip(1,3)+y*recip(2,3)+z*recip(3,3)
         fr1(n) = nfft1*(w1-ANINT(w1)+(0.5D0-SIGN(0.5D0,w1-ANINT(w1))))
         fr2(n) = nfft2*(w2-ANINT(w2)+(0.5D0-SIGN(0.5D0,w2-ANINT(w2))))
         fr3(n) = nfft3*(w3-ANINT(w3)+(0.5D0-SIGN(0.5D0,w3-ANINT(w3))))
      END DO
    END SUBROUTINE Fractionals
  END SUBROUTINE PME_
  FUNCTION Initialize_() RESULT(out)
    REAL(8), POINTER :: dummy(:,:,:)
    LOGICAL :: out
    out=.TRUE.

    nfft1=Ewald__Param % nx
    nfft2=Ewald__Param % ny
    nfft3=Ewald__Param % nz

    myfft1=nfft1/PI_npx
    myfft2=nfft2/PI_npy
    myfft3=nfft3/PI_npz
    alltoall_dim=myfft1*myfft2*myfft3/(PI_npx*PI_npy)

    rk_x=PI__Ranks_FFTW (PI_Node_FFTW+1) % nx
    rk_y=PI__Ranks_FFTW (PI_Node_FFTW+1) % ny
    rk_z=PI__Ranks_FFTW (PI_Node_FFTW+1) % nz

    my1_start=rk_x*myfft1+1
    my1_end=(rk_x+1)*myfft1
    my2_start=rk_y*myfft2+1
    my2_end=(rk_y+1)*myfft2
    my3_start=rk_z*myfft3+1
    my3_end=(rk_z+1)*myfft3
    istart=my1_start; iend=my1_end
    jstart=my2_start; jend=my2_end
    kstart=my3_start; kend=my3_end

    IF(PI_Node_FFTW ==0) WRITE(*,100)

    IF(.NOT. do_rfft3d(0,dummy,nfft1,nfft2,nfft3,nfftw3_start,nfftw3_local&
         &,nfftw2_start,nfftw2_local,ndim_fftw1,ndim_fftw2&
         &,ndim_fftw3,PI_Comm_FFTW)) RETURN
    nfftw1=nfft1
    nfftw2=nfft2
    nfftw3=ndim_fftw3

    kwstart=nfftw3_start
    kwend  =nfftw3_start+nfftw3_local-1
    CALL BSP_Moduli
100 FORMAT(/22x,'Finding optimal parameters for FFTWs.'/&
     &     22x,'     This will take a while...'/ /)     
  CONTAINS
    SUBROUTINE BSP_Moduli
      ALLOCATE(bsp_mod1(nfft1),bsp_mod2(nfft2),bsp_mod3(nfft3))
      CALL load_bsp_moduli(bsp_mod1,bsp_mod2,bsp_mod3,nfft1,nfft2&
           &,nfft3,Ewald__Param % order)
    END SUBROUTINE BSP_Moduli
  END FUNCTION Initialize_

  SUBROUTINE Transpose_Cart2FFTW
    INTEGER :: o,p,q,r,s,t,rk2_x,rk2_y,n,mzz

    
    CALL MPI_ALLTOALL(Cq_s,alltoall_dim,MPI_REAL8, Cq_sb,alltoall_dim&
         &,MPI_REAL8,PI_Comm_Z,ierr)


    mzz=myfft3/PI_Nprocs_Z
    DO n=1,PI_Nprocs_Z
       rk2_x=PI__Ranks_Z(n) % nx
       rk2_y=PI__Ranks_Z(n) % ny
       DO t=1,mzz
          q=(n-1)*mzz+t
          DO p=1,myfft2
             s=rk2_y*myfft2+p
             DO o=1,myfft1
                r=rk2_x*myfft1+o
                Cq_r(r,s,t)=Cq_sb(o,p,q)
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE Transpose_Cart2FFTW
  SUBROUTINE Transpose_FFTW2Cart
    INTEGER :: o,p,q,r,s,t,rk2_x,rk2_y,n,mzz

    mzz=myfft3/PI_Nprocs_Z
    DO n=1,PI_Nprocs_Z
       rk2_x=PI__Ranks_Z(n) % nx
       rk2_y=PI__Ranks_Z(n) % ny
       DO t=1,mzz
          q=(n-1)*mzz+t
          DO p=1,myfft2
             s=rk2_y*myfft2+p
             DO o=1,myfft1
                r=rk2_x*myfft1+o
                Cq_sb(o,p,q)=Cq_r(r,s,t)
             END DO
          END DO
       END DO
    END DO

    CALL MPI_ALLTOALL(Cq_sb,alltoall_dim,MPI_REAL8, Cq_s,alltoall_dim&
         &,MPI_REAL8,PI_Comm_Z,ierr)

  END SUBROUTINE Transpose_FFTW2Cart

  INCLUDE 'PME__Sources.f90'
END MODULE PME
