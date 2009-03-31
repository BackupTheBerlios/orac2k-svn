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
MODULE IntraAtoms
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Mon Jan 19 2009 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*

#include "parameters.h"
#ifdef HAVE_MPI
  USE mpi
#endif
  USE PI_Cutoffs
  USE PI_
  USE PI_Statistics, ONLY: PI__Write_Stats=>Write_It, PI__Time_It&
       &=>Time_It, PI__Sample_Exchange=>Sample_Exchange, PI__Add_Calls&
       &=>Add_Calls
  USE SystemTpg
  USE SystemPrm
  USE PrmUtilities
  USE Groups
  USE Atom
  USE Cell
  USE Constants, ONLY: max_pars,max_data,max_char
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, error_args, errmsg_f
  USE PointToPoint, ONLY: ResetSecondary_, ZeroSecondary_
  USE PI_Atom, ONLY: PI_Atom_Indexes
  USE Potential, ONLY: Reordering
  USE Forces, ONLY: Radii

  USE Print_Defs
  IMPLICIT none
  PRIVATE
  PUBLIC  IntraAtoms_,IntraParam, Param_Bonds,Param_Angles,Param_Imph&
       &,Param_Dihed,Param_Constr,Indx_Bonds,Indx_Angles,Indx_Imph&
       &,Indx_Dihed,Indx_Int14,Indx_Constr,xp0,yp0,zp0,fpx,fpy,fpz,Id&
       &,Slv,Charge,PI__ShiftIntra, Gather_Forces_,Gather_Coords_&
       &,Scatter_Coords_,Scatter_Forces_, Param_Int14

  type :: MyMaps_
     INTEGER :: Grp_No
     INTEGER :: g0=0
     INTEGER, ALLOCATABLE :: idx(:)
  end type MyMaps_
  TYPE :: Indx
     INTEGER :: NoAtm_s,NoAtm_r
     INTEGER, ALLOCATABLE :: ibuff_r(:)
     INTEGER, ALLOCATABLE :: ibuff_s(:)
  END type Indx
  TYPE :: iBuffer
     TYPE(Indx) :: sh(2)
  END type iBuffer
  TYPE :: IntraParam
     REAL(8), ALLOCATABLE :: Pot(:)
  END type IntraParam

!!$
!!$--- Used by IntraAtoms_
!!$
  REAL(8), DIMENSION(:), ALLOCATABLE, SAVE :: xp0,yp0,zp0,xpa,ypa,zpa,fpx,fpy,fpz
  INTEGER, ALLOCATABLE, SAVE :: Slv(:),Id(:)
  REAL(8), ALLOCATABLE, SAVE :: Charge(:)
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: Intra_t,Intra_p,Intra_conv,knwn
  Type(MyMaps_), ALLOCATABLE, TARGET, SAVE :: Map_n0(:),Map_n1(:)
  INTEGER, SAVE :: natom,natom_p

!!$
!!$--- Used by IndIntraBox_
!!$

  LOGICAL, ALLOCATABLE, SAVE :: oks(:),okt(:)
  TYPE(IntraParam), ALLOCATABLE, SAVE :: Param_Bonds(:),Param_Angles(:)&
       &,Param_Imph(:),Param_Dihed(:),Param_Constr(:),Param_Int14(:)
  INTEGER, ALLOCATABLE, SAVE :: Indx_Bonds(:,:),Indx_Angles(:,:)&
       &,Indx_Imph(:,:),Indx_Dihed(:,:),Indx_Int14(:,:),Indx_Constr(:,:)
!!$
!!$--- Used by shift routines
!!$

  TYPE(iBuffer), SAVE, TARGET :: iShift(6)
  INTEGER, SAVE :: Calls=0
  REAL(8), PARAMETER :: one=1.0D0,two=2.0D0,half=0.5D0
  REAL(8), SAVE :: startime,endtime,startime0,endtime0
  Type(MyMaps_), ALLOCATABLE, SAVE :: MyMaps(:)
  INTEGER, SAVE :: nMyMaps=0
  INTEGER, SAVE :: npx,npy,npz
  LOGICAL, ALLOCATABLE, SAVE :: ok_atm(:),ok_atm2(:)
  REAL(8), SAVE :: MyCutoff(2)=(/11.0_8,11.0_8/) !-- Need to compute it

CONTAINS
  FUNCTION IntraAtoms_() RESULT(out)
    LOGICAL :: out
    INTEGER, DIMENSION(:), ALLOCATABLE :: Extra
    INTEGER :: nextra,n,gg0,nn,Nextra_b

    out=.TRUE.

    CALL ResetSecondary_

!!$
!!$--- Atoms of the primary cell
!!$
    
    natom_p=COUNT(Atoms(:)%knwn == 1)
    IF(ALLOCATED(Intra_P)) CALL IntraAtoms__

!!$
!!$--- Prepare initial shift
!!$

!!$-- Atoms in shift maps can''t be counted twice 
    ALLOCATE(ok_Atm(SIZE(Atoms)),ok_Atm2(SIZE(Atoms)))
    ok_atm=.FALSE.
    ok_atm2=.FALSE.

    ALLOCATE(Intra_P(natom_p))
    Intra_P(:)=PACK((/1:SIZE(Atoms)/),Atoms(:)%knwn == 1)

!!$--- Prepare inclusion maps for atoms of the ghost cell 
!!$--- for shell n0 and n1

    CALL IntraMaps_n0_
    CALL IntraMaps_n1_

    ALLOCATE(knwn(SIZE(Atoms)))
    knwn=Atoms(:) % knwn


!!$--- Shift to half of the ghost cell count all interaction only once

    CALL ResetSecondary_
    CALL PI__ShiftIntra(_N0_,_INIT_,_COUNT_ONCE_)
    CALL PI__ShiftIntra(_N1_,_INIT_,_COUNT_ONCE_)

!!$--- knwn is two for atoms of the ghost cell contributing to the 
!!$--- Intra energy
    
    knwn=Atoms(:) % knwn

!!$--- Shift again to recover the full ghost cell. Full ghost contributes 
!!$--- to primary cell forces. Do not do Fold in this way

    CALL ResetSecondary_
    CALL PI__ShiftIntra(_N0_,_INIT_)
    CALL PI__ShiftIntra(_N1_,_INIT_)

    WHERE(Atoms(:) % knwn == 2 .AND. knwn == 0) knwn=-1

!!$--- natom is the No. of the primary and full ghost cell

    Nextra=COUNT(Atoms(:) % knwn == 2)
    natom=natom_p+Nextra

    ALLOCATE(Intra_T(natom))
    Intra_t(1:natom_p)=Intra_P
    Intra_t(natom_p+1:)=PACK((/1:SIZE(Atoms)/),Atoms(:) % knwn == 2)

!!$
!!$--- Check to see if all atoms needed are there
!!$

    
    Nextra_b=COUNT(ok_Atm2(:))
    IF(Nextra_b <= Nextra) THEN
       ok_Atm=.FALSE.
       nn=natom_p
       DO n=1,Nextra
          nn=nn+1
          ok_Atm(Intra_T(nn))=.TRUE.
       END DO
       DO n=1,SIZE(Atoms)
          IF(.NOT. ok_Atm2(n)) CYCLE
          IF(.NOT. ok_Atm(n)) THEN
             out=.FALSE.
             errmsg_f='Fatal error from IntraAtoms_: Not all needed gh&
                  &ost atoms are included.'
             CALL Add_Errors(-1,errmsg_f)
             RETURN
          END IF
       END DO
    ELSE
       out=.FALSE.
       errmsg_f='Fatal error from IntraAtoms_: Not all needed gh&
            &ost atoms are included.'
       CALL Add_Errors(-1,errmsg_f)
       RETURN
    END IF
!!$    IF(Reordering) CALL Labelling

!!$--- Intra_Conv gathers atoms to the simulation cell from Atoms(:) 

    ALLOCATE(Intra_Conv(SIZE(Atoms)))
    Intra_Conv=-1
    Intra_Conv(Intra_T(:))=(/1:natom/)

    IF(ALLOCATED(xp0)) THEN
       DEALLOCATE(xp0,yp0,zp0,fpx,fpy,fpz,Slv,Id,Charge)
    END IF
    ALLOCATE(xp0(natom),yp0(natom),zp0(natom))
    ALLOCATE(fpx(natom),fpy(natom),fpz(natom))
    ALLOCATE(Slv(natom),Id(natom),Charge(natom))

    Slv(:)=Atoms(Intra_t(:)) % Id_Slv
    Id(:)=Atoms(Intra_t(:)) % Id_Type
    Charge(:)=Atoms(Intra_t(:)) % chg

    IF(.NOT. IndIntraBox_n0_()) CALL Print_Errors()
    IF(.NOT. IndIntraBox_n1_()) CALL Print_Errors()

!!$
!!$--- Do not do any gathering to local atoms 
!!$

    CALL ZeroSecondary_
    out=.TRUE.
  CONTAINS
!!$
!!$--- Privates
!!$
    SUBROUTINE IntraMaps_n0_
      INTEGER :: n,n0,AtSt,AtEn,q,count0,nn,Grp_No
      INTEGER, ALLOCATABLE :: ind_o(:),ind_g(:)
      LOGICAL :: ok
      type(MyMaps_), ALLOCATABLE :: MyGrp(:)


      ALLOCATE(MyGrp(SIZE(Groupa)))

      DO n=1,SIZE(groupa)
         nn=Groupa(n) % AtEn-Groupa(n) % AtSt+1
         ALLOCATE(MyGrp(n) % Idx(nn))
      END DO

      
      n0=SIZE(Intra_P)
      ALLOCATE(ind_o(n0))
      ind_o=0
      count0=0
      DO nn=1,n0
         q=Intra_p(nn)
         ok=.FALSE.
         IF(Choose_Atom(q, Tpg % Bonds, Atoms_Tpg(q) % Bonds)) THEN
            ok=.TRUE.
         END IF
         IF(Choose_Atom(q, Tpg % Angles, Atoms_Tpg(q) % Angles)) THEN
            ok=.TRUE.
         END IF
         IF(Choose_Atom(q, Tpg % Imph, Atoms_Tpg(q) % Imph)) THEN
            ok=.TRUE.
         END IF
         IF(ok .AND. (.NOT. ok_atm(q))) THEN
            ok_atm(q)=.TRUE.
            count0=count0+1
            ind_o(count0)=q

            Grp_No=Atoms(q) % Grp_No
            MyGrp(Grp_No) % g0= MyGrp(Grp_No) % g0 +1
            n0=MyGrp(Grp_No) % g0
            MyGrp(Grp_No) % Idx(MyGrp(Grp_No) % g0)=q
         END IF
      END DO
      
      IF(ALLOCATED(Map_n0)) DEALLOCATE(Map_n0)

      nn=COUNT(MyGrp(:) % g0 /= 0)
      ALLOCATE(Map_n0(nn))

      n=0
      DO nn=1,SIZE(Groupa)
         IF(MyGrp(nn) % g0 /= 0) THEN
            n=n+1
            Map_n0(n) % Grp_no =nn
            Map_n0(n) % g0 =MyGrp(nn) % g0
            ALLOCATE(Map_n0(n) % idx(MyGrp(nn) %g0)) 
            Map_n0(n) % idx=MyGrp(nn) % idx
         END IF
      END DO
    END SUBROUTINE IntraMaps_n0_

    SUBROUTINE IntraMaps_n1_
      INTEGER :: n,n0,AtSt,AtEn,q,count0,nn,Grp_No
      INTEGER, ALLOCATABLE :: ind_o(:)
      LOGICAL :: ok
      type(MyMaps_), ALLOCATABLE :: MyGrp(:)


      ALLOCATE(MyGrp(SIZE(Groupa)))      
      DO n=1,SIZE(groupa)
         nn=Groupa(n) % AtEn-Groupa(n) % AtSt+1
         ALLOCATE(MyGrp(n) % Idx(nn))
      END DO

      n0=SIZE(Intra_P)
      
      ALLOCATE(ind_o(n0))
      ind_o=0
      count0=0
      DO nn=1,n0
         q=Intra_p(nn)
         ok=.FALSE.
         IF(Choose_Atom(q, Tpg % Dihed, Atoms_Tpg(q) % Dihed)) THEN
            ok=.TRUE.
         END IF
         IF(Choose_Atom(q, Tpg % Int14, Atoms_Tpg(q) % Int14)) THEN
            ok=.TRUE.
         END IF
         IF(ok .AND. (.NOT. ok_atm(q))) THEN
            ok_atm(q)=.TRUE.
            count0=count0+1
            ind_o(count0)=q
            Grp_No=Atoms(q) % Grp_No
            MyGrp(Grp_No) % g0= MyGrp(Grp_No) % g0 +1
            n0=MyGrp(Grp_No) % g0
            MyGrp(Grp_No) % Idx(MyGrp(Grp_No) % g0)=q
         END IF
      END DO

      IF(ALLOCATED(Map_n1)) DEALLOCATE(Map_n1)
      nn=COUNT(MyGrp(:) % g0 /= 0)
      ALLOCATE(Map_n1(nn))

      n=0
      DO nn=1,SIZE(Groupa)
         IF(MyGrp(nn) % g0 /= 0) THEN
            n=n+1
            Map_n1(n) % Grp_no =nn
            Map_n1(n) % g0 =MyGrp(nn) % g0
            ALLOCATE(Map_n1(n) % idx(MyGrp(nn) %g0)) 
            Map_n1(n) % idx=MyGrp(nn) % idx
         END IF
      END DO
    END SUBROUTINE IntraMaps_n1_
    FUNCTION Choose_Atom(q, Tpga, Tpgb) RESULT(out)
      LOGICAL :: out
      INTEGER :: q,Tpga(:,:), Tpgb(:)
      
      INTEGER :: n,p
      LOGICAL :: ok
      INTEGER :: r,t,s,ss,ssk
      
      
      
      n=SIZE(Tpgb)
      p=SIZE(Tpga,1)
      ok=.FALSE.
      
      DO r=1,n
         t=Tpgb(r)
         DO s=1,p
            ss=Tpga(s,t)
            IF(ss == q) CYCLE
            ssk=Atoms(ss) % knwn
            IF(ssk == 0) THEN
               ok=.TRUE.
               EXIT
            END IF
         END DO
         IF(ok) EXIT
      END DO
      out=ok

      DO r=1,n
         t=Tpgb(r)
         DO s=1,p
            ss=Tpga(s,t)
            IF(ss == q) CYCLE
            ssk=Atoms(ss) % knwn
            IF(ssk == 0 .AND. (.NOT. ok_atm2(ss))) THEN
               ok_atm2(ss)=.TRUE.
            END IF
         END DO
      END DO
    END FUNCTION Choose_Atom
 
  END FUNCTION IntraAtoms_

  SUBROUTINE IntraAtoms__
    DEALLOCATE(Intra_p,Intra_t,xp0,yp0,zp0,fpx,fpy,fpz,Slv,Id,Charge)
    DEALLOCATE(ok_atm,ok_atm2,knwn,Intra_Conv)
  END SUBROUTINE IntraAtoms__
  SUBROUTINE Gather_Coords_(x0,y0,z0)
    REAL(8) :: x0(:),y0(:),z0(:)
    xp0(:)=x0(Intra_t(:))
    yp0(:)=y0(Intra_t(:))
    zp0(:)=z0(Intra_t(:))
  END SUBROUTINE Gather_Coords_
  SUBROUTINE Gather_Forces_(x0,y0,z0)
    REAL(8) :: x0(:),y0(:),z0(:)
    fpx(:)=x0(Intra_t(:))
    fpy(:)=y0(Intra_t(:))
    fpz(:)=z0(Intra_t(:))
  END SUBROUTINE Gather_Forces_
  SUBROUTINE Scatter_Coords_(x0,y0,z0)
    REAL(8) :: x0(:),y0(:),z0(:)
    x0(Intra_t(:))=xp0(:)
    y0(Intra_t(:))=yp0(:)
    z0(Intra_t(:))=zp0(:)
  END SUBROUTINE Scatter_Coords_
  SUBROUTINE Scatter_Forces_(x0,y0,z0)
    REAL(8) :: x0(:),y0(:),z0(:)
    x0(Intra_t(:))=fpx(:)
    y0(Intra_t(:))=fpy(:)
    z0(Intra_t(:))=fpz(:)
  END SUBROUTINE Scatter_Forces_

#include "INTRAATOMS__Shift.f90"
#include "INTRAATOMS__Params.f90"
END MODULE IntraAtoms
