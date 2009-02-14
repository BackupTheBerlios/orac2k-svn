!!$/---------------------------------------------------------------------\
!!$   Copyright  � 2006-2007 Massimo Marchi <Massimo.Marchi at cea.fr>   |
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
MODULE Direct
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Mon Sep 29 2008 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*

#include "parameters.h"
  USE Erfc_Spline, ONLY: Erfc_,tau,c,E_dx=>dx,E_Xbeg=>Xbeg,Erfc_Switch
  USE Print_Defs
  USE PI_ATOM
  USE POTENTIAL
  USE Units
  USE Energies, ONLY: EN_Coul_Dir_=>Coul_Dir_,EN_LJ_=>LJ_
  USE Forces,ONLY: Force, fp_m,fp_l,fp_h, Forces_Init=>Init,Forces_Memory&
       &=>Memory, FORCES_Pick=>Pick, Radii
  USE LennardJones, ONLY: LennardJones__Par
  USE Cell, ONLY: oc,co, Volume
#ifdef HAVE_MPI
  USE mpi
#endif
  USE PI_
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f, errmsg_w
  USE IndBox
  USE PI_Communicate
  
  IMPLICIT none
  PRIVATE 
  PUBLIC Compute,Lists
  
  TYPE(Force), POINTER :: fp(:)
  INTEGER, SAVE :: No_Calls=0
  REAL(8), ALLOCATABLE, SAVE :: ecc6(:),ecc12(:),eccc(:),sig0_i(:),eps0(:),esig0(:)
  INTEGER, ALLOCATABLE, SAVE :: Id_ij(:,:)
  REAL(8), PARAMETER :: a1=0.2548296D0,a2=-0.28449674D0,a3&
       &=1.4214137D0,a4=-1.453152D0,a5=1.0614054D0,qp=0.3275911D0
  REAL(8), SAVE :: alphal,twrtpi,ewald_cut
  
CONTAINS
  SUBROUTINE Compute(i_pa,Initialize)
    INTEGER :: i_pa
    INTEGER, OPTIONAL :: Initialize
    INTEGER :: ierr,nnn,nn,n,m
    REAL(8), DIMENSION(:), ALLOCATABLE, SAVE :: fppx,fppy,fppz,Xg_PBC&
         &,Yg_PBC,Zg_PBC,Xgs_PBC,Ygs_PBC,Zgs_PBC,xcs,ycs,zcs,swrs&
         &,dswrs,cmap2,xmap3,ymap3,zmap3,Xc_PBC&
         &,Yc_PBC,Zc_PBC,Xcs_PBC,Ycs_PBC,Zcs_PBC
    INTEGER, ALLOCATABLE, SAVE :: IndGrp(:),IndGrps(:)&
         &,p_index_j(:),p_index_jj(:)
    INTEGER, ALLOCATABLE, SAVE :: neib(:),neic(:)
    REAL(8) :: startime,endtime,timea,ts1,te1,ts2,te2
    INTEGER, SAVE :: Times_of_Call=0
    Integer, Save :: NoCalls(3)=(/0.0_8,0.0_8,0.0_8/)
    INTEGER :: i_p


    i_p=i_pa-2

    Times_of_Call=Times_of_Call+1
    IF(PRESENT(Initialize)) THEN
       CALL Init
       RETURN
    END IF
    
    IF(ngroup == 0 .AND. natom == 0) THEN
       errmsg_f='Direct lattice forces routine must be called after PI&
            &_Atom_'
       CALL Add_Errors(-1,errmsg_f)
       CALL Print_Errors()
       STOP
    END IF

    fp=>FORCES_Pick(i_pa)

    CALL Memory
    CALL Forces
    No_Calls=No_Calls+1
    NoCalls(i_p)=NoCalls(i_p)+1
  CONTAINS
!!$
!!$--- Initialize things once for all
!!$
    SUBROUTINE Init
      INTEGER :: n,m,ij,ncutoff
      REAL(8) :: ecc_R,aux,ene,rcut,x,bin
      
      twrtpi=2.0d0/SQRT(pi)
      alphal = Ewald__Param % alpha
      ecc_R=-5.0D-1/(efact/1000.0D0)
      n=SIZE(LennardJones__Par % Par_SE)
      ALLOCATE(Id_ij(n,n))
      ALLOCATE(ecc6(n*(n+1)/2),ecc12(n*(n+1)/2),eccc(n*(n+1)/2))
      ALLOCATE(sig0_i(n*(n+1)/2),eps0(n*(n+1)/2),esig0(n*(n+1)/2))
      DO n=1,SIZE(LennardJones__Par % Par_SE)
         DO m=n,SIZE(LennardJones__Par % Par_SE)
            ij=m*(m-1)/2+n
            Id_ij(n,m)=ij
            Id_ij(m,n)=ij

            ecc6(ij)=LennardJones__Par % c6(ij)
            ecc12(ij)=LennardJones__Par % c12(ij)
            ene=0.0D0
            IF(ecc6(ij) /= 0.0D0 .AND. ecc12(ij) /= 0.0D0) THEN
               aux=(ecc12(ij)/ecc6(ij))
               eccc(ij)=3.5_8*aux**(1.0D0/6.0D0)
               sig0_i(ij)=1.0D0/aux**(1.0D0/6.0D0)
               eps0(ij)=ecc6(ij)/aux
               esig0(ij)=sig0_i(ij)*eps0(ij)
            ELSE
               eccc(ij)=0.0D0
            END IF
            eccc(ij)=eccc(ij)**2
         END DO

      END DO

      x=1.0_8
      bin=0.01D0
      DO 
         x=x+bin
         IF(ERFC(alphal*x) < 7.0D-6) EXIT
      END DO
      ewald_cut=x*x

      IF(Erfc_Switch) THEN
         ncutoff=SIZE(Radii)
         rcut=Radii(ncutoff) % out + Radii(ncutoff) % update
         CALL Erfc_(rcut,Ewald__Param % alpha)
      END IF
    END SUBROUTINE Init
!!$
!!$--- Get Memory
!!$
    SUBROUTINE Memory
      IF(ALLOCATED(Xg_PBC)) THEN
         DEALLOCATE(Xg_PBC,Yg_PBC,Zg_PBC,Xgs_PBC&
              &,Ygs_PBC,Zgs_PBC,IndGrp&
              &,IndGrps,neib,xcs,ycs&
              &,zcs,swrs,dswrs,cmap2&
              &,xmap3,ymap3,zmap3,neic,Xc_PBC&
              &,Yc_PBC,Zc_PBC,Xcs_PBC,Ycs_PBC,Zcs_PBC)
         DEALLOCATE(fppx,fppy,fppz,p_index_j&
              &,p_index_jj) 
      END IF

      ALLOCATE(Xg_PBC(ngroup),Yg_PBC(ngroup),Zg_PBC(ngroup),Xgs_PBC(ngroup)&
           &,Ygs_PBC(ngroup),Zgs_PBC(ngroup),IndGrp(ngroup)&
           &,IndGrps(ngroup),neib(ngroup),xcs(ngroup),ycs(ngroup)&
           &,zcs(ngroup),swrs(ngroup),dswrs(ngroup),cmap2(ngroup)&
           &,xmap3(ngroup),ymap3(ngroup),zmap3(ngroup),neic(ngroup)&
           &,Xc_PBC(ngroup),Yc_PBC(ngroup),Zc_PBC(ngroup)&
           &,Xcs_PBC(ngroup),Ycs_PBC(ngroup),Zcs_PBC(ngroup))
      
      ALLOCATE(fppx(natom),fppy(natom),fppz(natom),p_index_j(natom)&
           &,p_index_jj(natom)) 

    END SUBROUTINE Memory
!!$
!!$--- Compute Forces
!!$

#include "DIRECT__Sources.f90"
  END SUBROUTINE Compute

  SUBROUTINE Lists(NShell)
    INTEGER :: NShell
    INTEGER :: ierr,nnn,nn,n,m
    INTEGER, ALLOCATABLE :: neib(:),neic(:)
    REAL(8) :: startime,endtime,timea,ts1,te1,ts2,te2
    INTEGER, SAVE :: Times_of_Call=0
    INTEGER :: i_p


    Times_of_Call=Times_of_Call+1
    IF(ngroup == 0 .AND. natom == 0) THEN
       errmsg_f='Direct lattice Lists routine must be called after PI&
            &_Atom_'
       CALL Add_Errors(-1,errmsg_f)
       CALL Print_Errors()
       STOP
    END IF

    ALLOCATE(neib(ngroup),neic(ngroup))
    CALL Lists_(Nshell)
    No_Calls=No_Calls+1
  CONTAINS
#include "DIRECT__Lists.f90"
  END SUBROUTINE Lists
END MODULE Direct
