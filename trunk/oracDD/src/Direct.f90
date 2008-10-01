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

  USE Potential
  USE Units
  USE Forces, ONLY: Force, rcut_i,rcut_o
  USE LennardJones, ONLY: LennardJones__Par
  USE Cell, ONLY: oc,co, Volume
#ifdef HAVE_MPI
  USE mpi
#endif
  USE PI_
  USE Groups
  USE Atom
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f, errmsg_w
  USE IndBox
  USE Neighbors_S, ONLY: Neighbors__Particles,Neighbors_S__nc, &
       & Neighbors_S__Chain, Neighbors_S__Ind,Chain_xyz&
       &, Head_xyz, nc, clst
  USE PI_Communicate
  
  IMPLICIT none
  PRIVATE 
  PUBLIC Compute, fp
  
  TYPE(Force), ALLOCATABLE, SAVE :: fp(:)
  INTEGER, SAVE :: No_Calls=0
  REAL(8), ALLOCATABLE, SAVE :: ecc6(:),ecc12(:)
  REAL(8), PARAMETER :: a1=0.2548296D0,a2=-0.28449674D0,a3&
       &=1.4214137D0,a4=-1.453152D0,a5=1.0614054D0,qp=0.3275911D0&
       &,twrtpi=2.0d0/SQRT(pi)
  REAL(8), SAVE :: alphal
CONTAINS
  SUBROUTINE Compute(i_p)
    INTEGER :: i_p
    INTEGER :: natom,ngroup,ncx,ncy,ncz

    REAL(8), ALLOCATABLE :: xpg(:),ypg(:),zpg(:),xp0(:),yp0(:),zp0(:)
    REAL(8), DIMENSION(:), ALLOCATABLE :: chg,fppx,fppy,fppz,Xg_PBC&
         &,Yg_PBC,Zg_PBC,Xgs_PBC,Ygs_PBC,Zgs_PBC 
    INTEGER, ALLOCATABLE :: Id(:),Slv(:),IndGrp(:),IndGrps(:)&
         &,p_index_j(:),p_index_jj(:)
    INTEGER, ALLOCATABLE :: grppt(:,:)
    LOGICAL, ALLOCATABLE :: maplg(:)
    TYPE :: Mapnl
       INTEGER, ALLOCATABLE :: ex(:)
    END type Mapnl
    TYPE(Mapnl), ALLOCATABLE :: Maps(:)
    INTEGER, POINTER :: Index_0(:)
    TYPE(Neighbors_S__Ind), POINTER :: ind_xyz(:)
    REAL(8) :: startime,endtime,timea,ts1,te1,ts2,te2


    IF(No_Calls /= 0) THEN
       DEALLOCATE(xpg,ypg,zpg,xp0,yp0,zp0,chg,Id,Slv,grppt,Maps,fppx&
            &,fppy,fppz,maplg,IndGrp,IndGrps,p_index_j,p_index_jj)
       DEALLOCATE(Xg_PBC,Yg_PBC,Zg_PBC,Xgs_PBC,Ygs_PBC,Zgs_PBC)
    ELSE
       CALL Init
    END IF

    CALL Memory

    CALL Gather_Atoms

    IF(.NOT. Neighbors__Particles(i_p,xpg,ypg,zpg)) CALL Print_Errors()

    Ind_xyz=>clst(i_p) % Ind_xyz

    ncx = nc(i_p) % x; ncy = nc(i_p) % y; ncz = nc(i_p) % z

    CALL MPI_BARRIER(PI_Comm,ierr)
    startime=MPI_WTIME()
    CALL Forces
    endtime=MPI_WTIME()
    timea=endtime-startime
    WRITE(*,*) 'time for Force ',timea

    No_Calls=No_Calls+1
  CONTAINS
!!$
!!$--- Initialize once for all
!!$
    SUBROUTINE Init
      INTEGER :: n,m,ij
      
      alphal = Ewald__Param % alpha
      n=SIZE(LennardJones__Par % Par_SE)
      ALLOCATE(ecc6(n*(n+1)/2),ecc12(n*(n+1)/2))
      DO n=1,SIZE(LennardJones__Par % Par_SE)
         DO m=n,SIZE(LennardJones__Par % Par_SE)
            ij=m*(m-1)/2+n
            ecc6(ij)=LennardJones__Par % c6(ij)
            ecc12(ij)=LennardJones__Par % c12(ij)
         END DO
      END DO
    END SUBROUTINE Init
!!$
!!$--- Get Memory
!!$
    SUBROUTINE Memory
      IF(.NOT. IndBox_()) CALL Print_Errors()
      natom=SIZE(IndBox_a_t) ; ngroup=SIZE(IndBox_g_t)
      
      ALLOCATE(xpg(ngroup),ypg(ngroup),zpg(ngroup),grppt(2,ngroup))
      ALLOCATE(Xg_PBC(ngroup),Yg_PBC(ngroup),Zg_PBC(ngroup),Xgs_PBC(ngroup)&
           &,Ygs_PBC(ngroup),Zgs_PBC(ngroup),IndGrp(ngroup),IndGrps(ngroup))
      
      ALLOCATE(xp0(natom),yp0(natom),zp0(natom),chg(natom),Id(natom)&
           &,Slv(natom),Maps(natom),fppx(natom),fppy(natom),fppz(natom)&
           &,maplg(natom),p_index_j(natom),p_index_jj(natom))
      
    END SUBROUTINE Memory
!!$
!!$--- Gather Atoms to the CPU box
!!$
    SUBROUTINE Gather_Atoms
      INTEGER :: n,m,p,q,nn,count0,g1,g2
      DO n=1,ngroup
         m=IndBox_g_t(n)
         xpg(n)=Groupa(m) % xa
         ypg(n)=Groupa(m) % ya
         zpg(n)=Groupa(m) % za
         grppt(1,n)=Groupa(m) % AtSt
         grppt(2,n)=Groupa(m) % AtEn
      END DO
      DO n=1,natom
         m=IndBox_a_t(n)
         xp0(n)=Atoms(m) % xa
         yp0(n)=Atoms(m) % ya
         zp0(n)=Atoms(m) % za
         chg(n)=Atoms(m) % chg
         Id(n)=Atoms(m) % Id_Type
         Slv(n)=Atoms(m) % Id_Slv
      END DO
      
      ALLOCATE(Index_0(SIZE(Atoms)))
      Index_0=-1
      DO n=1,natom
         m=IndBox_a_t(n)
         Index_0(m)=n
      END DO
!!$
!!$--- Change grppt to resect atoms known in the box
!!$
      DO n=1,ngroup
         g1=grppt(1,n)
         g2=grppt(2,n)
         IF(Index_0(g1) /= -1) THEN
            grppt(1,n)=Index_0(g1)
         ELSE
            WRITE(*,*) 'Weired!!!'
            STOP
         END IF
         IF(Index_0(g2) /= -1) THEN
            grppt(2,n)=Index_0(g2)
         ELSE
            WRITE(*,*) 'Weired!!!'
            STOP
         END IF
      END DO
      
!!$
!!$--- Change exclusion atoms 
!!$
      
      DO n=1,natom
         m=IndBox_a_t(n)
         nn=SIZE(Atoms_Tpg(m) % Ex)
         count0=0
         DO q=1,nn
            p=Atoms_Tpg(m) % Ex(q)
            IF(Index_0(p) /= -1) THEN
               count0=count0+1
            END IF
         END DO
         IF(count0 == 0) CYCLE
         ALLOCATE(Maps(n) % ex(count0))
         count0=0
         DO q=1,nn
            p=Atoms_Tpg(m) % Ex(q)
            IF(Index_0(p) /= -1) THEN
               count0=count0+1
               Maps(n) % Ex(count0)=Index_0(Atoms_Tpg(m) % Ex(q))
            END IF
         END DO
      END DO
      
    END SUBROUTINE Gather_Atoms

!!$
!!$--- Compute Forces
!!$

    INCLUDE 'DIRECT__Sources.f90'
  END SUBROUTINE Compute
END MODULE Direct
