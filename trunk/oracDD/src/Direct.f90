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

  USE PI_Atom
  USE Potential
  USE Units
  USE Forces
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
  PUBLIC Compute, fp
  
  TYPE(Force), ALLOCATABLE, SAVE :: fp(:)
  INTEGER, SAVE :: No_Calls=0
  REAL(8), ALLOCATABLE, SAVE :: ecc6(:),ecc12(:)
  INTEGER, ALLOCATABLE, SAVE :: Id_ij(:,:)
  REAL(8), PARAMETER :: a1=0.2548296D0,a2=-0.28449674D0,a3&
       &=1.4214137D0,a4=-1.453152D0,a5=1.0614054D0,qp=0.3275911D0
  REAL(8), SAVE :: alphal,twrtpi
CONTAINS
  SUBROUTINE Compute(i_p,Initialize)
    INTEGER :: i_p
    INTEGER, OPTIONAL :: Initialize
    INTEGER :: ierr,nnn,nn,n
    REAL(8), DIMENSION(:), ALLOCATABLE :: fppx,fppy,fppz,Xg_PBC&
         &,Yg_PBC,Zg_PBC,Xgs_PBC,Ygs_PBC,Zgs_PBC,xcs,ycs,zcs,swrs&
         &,dswrs,cmap2,xmap3,ymap3,zmap3
    INTEGER, ALLOCATABLE :: IndGrp(:),IndGrps(:)&
         &,p_index_j(:),p_index_jj(:)
    INTEGER, ALLOCATABLE :: nei(:)
    REAL(8) :: startime,endtime,timea,ts1,te1,ts2,te2

    IF(ngroup == 0 .AND. natom == 0) THEN
       errmsg_f='Direct lattice forces routine must be called after PI&
            &_Atom_'
       CALL Print_Errors()
       STOP
    END IF
    IF(PRESENT(Initialize)) THEN
       CALL Init
       RETURN
    END IF

    CALL Memory
    CALL Forces
    CALL PI__Fold_F(fp,i_p)

    No_Calls=No_Calls+1
  CONTAINS
!!$
!!$--- Initialize things once for all
!!$
    SUBROUTINE Init
      INTEGER :: n,m,ij
      
      twrtpi=2.0d0/SQRT(pi)
      alphal = Ewald__Param % alpha

      n=SIZE(LennardJones__Par % Par_SE)
      ALLOCATE(Id_ij(n,n))
      ALLOCATE(ecc6(n*(n+1)/2),ecc12(n*(n+1)/2))
      DO n=1,SIZE(LennardJones__Par % Par_SE)
         DO m=n,SIZE(LennardJones__Par % Par_SE)
            ij=m*(m-1)/2+n
            Id_ij(n,m)=ij
            Id_ij(m,n)=ij
            ecc6(ij)=LennardJones__Par % c6(ij)
            ecc12(ij)=LennardJones__Par % c12(ij)
         END DO
      END DO
    END SUBROUTINE Init
!!$
!!$--- Get Memory
!!$
    SUBROUTINE Memory
      ALLOCATE(Xg_PBC(ngroup),Yg_PBC(ngroup),Zg_PBC(ngroup),Xgs_PBC(ngroup)&
           &,Ygs_PBC(ngroup),Zgs_PBC(ngroup),IndGrp(ngroup)&
           &,IndGrps(ngroup),nei(ngroup),xcs(ngroup),ycs(ngroup)&
           &,zcs(ngroup),swrs(ngroup),dswrs(ngroup),cmap2(ngroup)&
           &,xmap3(ngroup),ymap3(ngroup),zmap3(ngroup))
      
      ALLOCATE(fppx(natom),fppy(natom),fppz(natom),p_index_j(natom)&
           &,p_index_jj(natom)) 
      IF(ALLOCATED(fp)) DEALLOCATE(fp)
      ALLOCATE(fp(natom))
      fp(:) % x =0.0_8
      fp(:) % y =0.0_8
      fp(:) % z =0.0_8
    END SUBROUTINE Memory
!!$
!!$--- Compute Forces
!!$

    INCLUDE 'DIRECT__Sources.f90'
  END SUBROUTINE Compute
END MODULE Direct
