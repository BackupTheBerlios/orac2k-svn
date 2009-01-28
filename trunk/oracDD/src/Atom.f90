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
MODULE Atom
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Jan 25 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*

#include "config.h"
#define _IND  IndBox_a_p(:)
#include 'forces.h'

#ifdef HAVE_MPI
  USE mpi
#endif
  USE PI_
  USE PI_Collectives
  USE Constants
  USE Units, ONLY: unitc,unite,gascon,Boltz,efact
  USE Errors, ONLY: Add_Errors=>Add, errmsg_f,Print_Errors
  USE SystemTpg
  USE Cell
  USE AtomCnt
  USE AtomBox
  USE SimulationBox
  USE Simulation, ONLY: T, Temp_
  USE PDB
  USE GAUSS
  USE LA_Routines
  USE IndBox, ONLY: IndBox_a_p,IndBox_a_t
  USE Forces, Forces__Pick=>Pick
  USE Energies, ONLY: EN_Kinetic_=>Kinetic_,EN_Temp_=>Temp_
  IMPLICIT none
  PRIVATE
  PUBLIC Atom_,Atom__Tpg_, Atom__PDB, Atom__InitCoords,Atom__,&
       & Atom__Tpg,Atoms_Tpg, Atoms,Atom__vInit_, Atom__Verlet_,&
       & Atom__Correct_,Atom__Write_, Atom__Convert,natom_Slt&
       &,natom_Slv,Atom__KinCompute
  TYPE :: Atom__
     REAL(8) :: x,y,z      ! Coordinates orthogonal frame
     REAL(8) :: xa,ya,za   ! Coordinates reduced frame
     REAL(8) :: vx,vy,vz      ! Velocities orthogonal frame
     REAL(8) :: vxa,vya,vza   ! Velocities reduced frame
     REAL(8) :: chg, mass,pmass
     CHARACTER(len=max_atm) :: Res, beta, betab
     INTEGER :: Res_No, Grp_No, Id_Res, Id_Type, Id_slv, Mol, knwn
  END TYPE Atom__
  TYPE :: Atom__Tpg
     INTEGER, ALLOCATABLE :: Bonds(:)
     INTEGER, ALLOCATABLE :: Angles(:)
     INTEGER, ALLOCATABLE :: dihed(:)
     INTEGER, ALLOCATABLE :: imph(:)
     INTEGER, ALLOCATABLE :: int14(:)
     INTEGER, ALLOCATABLE :: Ex(:)
  END type Atom__Tpg

  TYPE :: Atom__Tpg0
     INTEGER, ALLOCATABLE :: Bnd(:)
  END type Atom__Tpg0

  TYPE(Atom__), ALLOCATABLE, TARGET, SAVE :: Atoms(:)
  TYPE(Atom__Tpg), ALLOCATABLE, SAVE :: Atoms_Tpg(:)
  REAL(8), SAVE :: Temperature
  INTEGER, SAVE :: natom_Slt=0,natom_Slv=0
CONTAINS
  FUNCTION Atom_() RESULT(out)
    LOGICAL :: out
    INTEGER :: n,m,s,nato
    REAL(8), ALLOCATABLE :: massa(:)

    out=.TRUE.
    IF(.NOT. ALLOCATED(AtomCnts)) THEN
       errmsg_f='The No. of Atoms of the system is unknown'
       CALL Add_Errors(-1,errmsg_f)
       out=.FALSE.
       RETURN
    END IF
    nato=SIZE(AtomCnts)
    ALLOCATE(massa(nato))
    massa=0.0D0
    DO n=1,nato
       m=AtomCnts(n) % Grp_No
       massa(m)=massa(m)+AtomCnts(n) % mass
    END DO

    ALLOCATE(Atoms(nato))

    DO n=1,nato
       m=AtomCnts(n) % Grp_No
       Atoms(n) % Grp_No = AtomCnts(n) % Grp_No
       Atoms(n) % Res_No = AtomCnts(n) % Res_No
       Atoms(n) % Id_Type = AtomCnts(n) % Id_Type
       Atoms(n) % Id_Slv = AtomCnts(n) % Id_Slv
       Atoms(n) % Id_Res = AtomCnts(n) % Id_Res
       Atoms(n) % chg = AtomCnts(n) % chg/SQRT(unitc)
       Atoms(n) % mass = AtomCnts(n) % mass
       Atoms(n) % pmass = AtomCnts(n) % mass/massa(m)
       Atoms(n) % Res = AtomCnts(n) % Res
       Atoms(n) % beta = AtomCnts(n) % Beta
       Atoms(n) % betab = AtomCnts(n) % Betab
    END DO
    DO n=1,SIZE(Tpg % Mol_Atm)
       DO m=1,SIZE(Tpg % Mol_Atm(n) % g)
          s=Tpg % Mol_Atm(n) % g (m)
          Atoms(s) % Mol = n
       END DO
    END DO    
    natom_Slv=COUNT(Atoms(:) % Id_Slv == 2)
    natom_Slt=COUNT(Atoms(:) % Id_Slv == 1)
  END FUNCTION Atom_
  FUNCTION Atom__InitCoords() RESULT(out)
    LOGICAL :: out
    INTEGER :: n

    out=.TRUE.
    IF(.NOT. ALLOCATED(Atoms_InBox)) THEN
       errmsg_f='Cannot initialize atomic coordinates: Initial&
            & coordinates were expected, but none were either &
            &built or read in by SimulationBox Module.'
       CALL Add_Errors(-1,errmsg_f)
       out=.FALSE.
       RETURN
    END IF
    
    Atoms(:) % xa = Atoms_InBox(:) % x 
    Atoms(:) % ya = Atoms_InBox(:) % y 
    Atoms(:) % za = Atoms_InBox(:) % z 
    Atoms(:) % x = co(1,1)*Atoms(:) % xa+co(1,2)*Atoms(:) % ya+co(1,3)*Atoms(:) % za    
    Atoms(:) % y = co(2,1)*Atoms(:) % xa+co(2,2)*Atoms(:) % ya+co(2,3)*Atoms(:) % za    
    Atoms(:) % z = co(3,1)*Atoms(:) % xa+co(3,2)*Atoms(:) % ya+co(3,3)*Atoms(:) % za    
  END FUNCTION Atom__InitCoords
  FUNCTION Atom__Convert(Dir,Ind) RESULT(out)
    INTEGER, OPTIONAL :: Ind(:)

    LOGICAL :: out
    INTEGER :: Dir
    out=.TRUE.
    IF(PRESENT(Ind)) THEN
       SELECT CASE(Dir)
       CASE(_XA_TO_X_)
          Atoms(:) % x = co(1,1)*Atoms(:) % xa+co(1,2)*Atoms(:) % ya+co(1,3)*Atoms(:) % za    
          Atoms(:) % y = co(2,1)*Atoms(:) % xa+co(2,2)*Atoms(:) % ya+co(2,3)*Atoms(:) % za    
          Atoms(:) % z = co(3,1)*Atoms(:) % xa+co(3,2)*Atoms(:) % ya+co(3,3)*Atoms(:) % za    
       CASE(_X_TO_XA_)
          Atoms(:) % xa = oc(1,1)*Atoms(:) % x+oc(1,2)*Atoms(:) % y+oc(1,3)*Atoms(:) % z    
          Atoms(:) % ya = oc(2,1)*Atoms(:) % x+oc(2,2)*Atoms(:) % y+oc(2,3)*Atoms(:) % z    
          Atoms(:) % za = oc(3,1)*Atoms(:) % x+oc(3,2)*Atoms(:) % y+oc(3,3)*Atoms(:) % z
       END SELECT
    ELSE
       SELECT CASE(Dir)
       CASE(_XA_TO_X_)
          Atoms(Ind(:)) % x = co(1,1)*Atoms(Ind(:)) % xa+co(1,2)&
               &*Atoms(Ind(:)) % ya+co(1,3)*Atoms(Ind(:)) % za 
          Atoms(Ind(:)) % y = co(2,1)*Atoms(Ind(:)) % xa+co(2,2)&
               &*Atoms(Ind(:)) % ya+co(2,3)*Atoms(Ind(:)) % za 
          Atoms(Ind(:)) % z = co(3,1)*Atoms(Ind(:)) % xa+co(3,2)&
               &*Atoms(Ind(:)) % ya+co(3,3)*Atoms(Ind(:)) % za
       CASE(_X_TO_XA_)
          Atoms(Ind(:)) % xa = oc(1,1)*Atoms(Ind(:)) % x+oc(1,2)&
               &*Atoms(Ind(:)) % y+oc(1,3)*Atoms(Ind(:)) % z
          Atoms(Ind(:)) % ya = oc(2,1)*Atoms(Ind(:)) % x+oc(2,2)&
               &*Atoms(Ind(:)) % y+oc(2,3)*Atoms(Ind(:)) % z
          Atoms(Ind(:)) % za = oc(3,1)*Atoms(Ind(:)) % x+oc(3,2)&
               &*Atoms(Ind(:)) % y+oc(3,3)*Atoms(Ind(:)) % z
       END SELECT
    END IF
  END FUNCTION Atom__Convert
  SUBROUTINE Atom__PDB(unit, nozero_write)
    INTEGER, OPTIONAL :: nozero_write
    INTEGER :: unit
    TYPE(AtomPdb), ALLOCATABLE :: PDB__Coords(:)
    INTEGER :: n,ierr,nn

    ALLOCATE(PDB__Coords(SIZE(Atoms)))
    
    FORALL(n=1:SIZE(Atoms)) PDB__Coords(n) % Serial = n
    DO nn=1,SIZE(Indbox_a_p)
       n=Indbox_a_p(nn)
       PDB__Coords(n) % x = Atoms(n) % x
       PDB__Coords(n) % y = Atoms(n) % y
       PDB__Coords(n) % z = Atoms(n) % z
    END DO
    CALL PI_Gather_(PDB__Coords(:) % x, PDB__Coords(:) % y, PDB__Coords(:) % z)
    IF(Pi_node_cart == 0) THEN
       IF(PRESENT(nozero_write)) THEN
          CALL PDB__Write(unit,PDB__Coords, nozero_write)
       ELSE
          CALL PDB__Write(unit,PDB__Coords)
       END IF
    END IF
#ifdef HAVE_MPI
    CALL MPI_BARRIER(PI_Comm_cart,ierr)
#endif    
  END SUBROUTINE Atom__PDB
  FUNCTION Atom__Tpg_() RESULT(out)
    LOGICAL :: out
    INTEGER :: natom
    TYPE(Atom__Tpg0), POINTER :: Tpg_out(:)
    INTEGER, POINTER :: indx(:),Maps(:)
    INTEGER :: n,nn,m,ia,ib,ic,count_n
    LOGICAL, ALLOCATABLE :: ok(:)
    natom=SIZE(Atoms)
    ALLOCATE(Atoms_Tpg(natom))
    ALLOCATE(Tpg_out(natom))
    
    ALLOCATE(Indx(natom))
    ALLOCATE(ok(natom))
!!$-- Bonds

    CALL Add_Tpg(Tpg % Bonds, tpg_out)
    DO n=1,natom
       IF(ALLOCATED(Tpg_out(n) % Bnd)) THEN
          nn=SIZE(Tpg_out(n) % Bnd)
          ALLOCATE(Atoms_Tpg(n) % Bonds(nn))
          Atoms_Tpg(n) % Bonds=Tpg_out(n) % Bnd
          DEALLOCATE(Tpg_out(n) % Bnd)
       END IF
    END DO

!!$-- Angles

    CALL Add_Tpg(Tpg % Angles, tpg_out)

    DO n=1,natom
       IF(ALLOCATED(Tpg_out(n) % Bnd)) THEN
          nn=SIZE(Tpg_out(n) % Bnd)
          ALLOCATE(Atoms_Tpg(n) % Angles(nn))
          Atoms_Tpg(n) % Angles=Tpg_out(n) % Bnd
          DEALLOCATE(Tpg_out(n) % Bnd)
       END IF
    END DO

!!$-- Diheds

    CALL Add_Tpg(Tpg % Dihed, tpg_out)

    DO n=1,natom
       IF(ALLOCATED(Tpg_out(n) % Bnd)) THEN
          nn=SIZE(Tpg_out(n) % Bnd)
          ALLOCATE(Atoms_Tpg(n) % Dihed(nn))
          Atoms_Tpg(n) % Dihed=Tpg_out(n) % Bnd
          DEALLOCATE(Tpg_out(n) % Bnd)
       END IF
    END DO

!!$-- Imph

    CALL Add_Tpg(Tpg % Imph, tpg_out)

    DO n=1,natom
       IF(ALLOCATED(Tpg_out(n) % Bnd)) THEN
          nn=SIZE(Tpg_out(n) % Bnd)
          ALLOCATE(Atoms_Tpg(n) % Imph(nn))
          Atoms_Tpg(n) % Imph=Tpg_out(n) % Bnd
          DEALLOCATE(Tpg_out(n) % Bnd)
       END IF
    END DO

!!$-- int14

    CALL Add_Tpg(Tpg % Int14, tpg_out)

    DO n=1,natom
       IF(ALLOCATED(Tpg_out(n) % Bnd)) THEN
          nn=SIZE(Tpg_out(n) % Bnd)
          ALLOCATE(Atoms_Tpg(n) % Int14(nn))
          Atoms_Tpg(n) % Int14=Tpg_out(n) % Bnd
          DEALLOCATE(Tpg_out(n) % Bnd)
       END IF
    END DO


    ALLOCATE(Maps(200))
    DO n=1,natom
       Maps=0
       count_n=0
       ok=.TRUE.
       DO nn=1,SIZE(Atoms_Tpg(n) % Bonds)
          m=Atoms_Tpg(n) % Bonds(nn)
          ia=Tpg % Bonds(1,m)
          ib=Tpg % Bonds(2,m)
          IF(n /= ia .AND. ok(ia)) THEN
             count_n=count_n+1
             Maps(count_n)=ia
             ok(ia)=.FALSE.
          END IF
          IF(n /= ib .AND. ok(ib)) THEN
             count_n=count_n+1
             Maps(count_n)=ib
             ok(ib)=.FALSE.
          END IF
       END DO
       DO nn=1,SIZE(Atoms_Tpg(n) % Angles)
          m=Atoms_Tpg(n) % Angles(nn)
          ia=Tpg % Angles(1,m)
          ib=Tpg % Angles(2,m)
          ic=Tpg % Angles(3,m)
          IF(n /= ia .AND. ok(ia)) THEN
             count_n=count_n+1
             Maps(count_n)=ia
             ok(ia)=.FALSE.
          END IF
          IF(n /= ib .AND. ok(ib)) THEN
             count_n=count_n+1
             Maps(count_n)=ib
             ok(ib)=.FALSE.
          END IF
          IF(n /= ic .AND. ok(ic)) THEN
             count_n=count_n+1
             Maps(count_n)=ic
             ok(ic)=.FALSE.
          END IF
       END DO
       DO nn=1,SIZE(Atoms_Tpg(n) % Int14)
          m=Atoms_Tpg(n) % Int14(nn)
          ia=Tpg % Int14(1,m)
          ib=Tpg % Int14(2,m)
          IF(n /= ia .AND. ok(ia)) THEN
             count_n=count_n+1
             Maps(count_n)=ia
             ok(ia)=.FALSE.
          END IF
          IF(n /= ib .AND. ok(ib)) THEN
             count_n=count_n+1
             Maps(count_n)=ib
             ok(ib)=.FALSE.
          END IF
       END DO
       ALLOCATE(Atoms_Tpg(n) % Ex (count_n))
       Atoms_Tpg(n) % Ex (:)=Maps(1:count_n)
    END DO

    out=.TRUE.


  CONTAINS
    SUBROUTINE Add_Tpg(tpg_in, tpg_out)
      INTEGER :: tpg_in(:,:)
      TYPE(Atom__Tpg0) :: Tpg_out(:)

      INTEGER :: d1,d2,n,nn,m

      Indx=0
      d1=SIZE(Tpg_in,1)
      d2=SIZE(Tpg_in,2)

      DO n=1,d2
         DO m=1,d1
            nn=Tpg_in(m,n)
            Indx(nn)=Indx(nn)+1
         END DO
      END DO

      DO n=1,natom
         IF(Indx(n) /= 0) THEN
            ALLOCATE(Tpg_out(n) % Bnd(Indx(n)))
         END IF
      END DO

      Indx=0
      DO n=1,d2
         DO m=1,d1
            nn=Tpg_in(m,n)
            Indx(nn)=Indx(nn)+1
            Tpg_Out(nn) % Bnd(Indx(nn)) = n
         END DO
      END DO
    END SUBROUTINE Add_Tpg
  END FUNCTION Atom__Tpg_

  FUNCTION Atom__Verlet_(dt,level) RESULT(out)
    INTEGER :: level
    REAL(8) :: dt
    LOGICAL :: out
    TYPE(Force), POINTER :: fp(:)
    INTEGER :: n,nn,m
    REAL(8) :: ts2,tfact,mass


    ts2=dt*dt

    out=.TRUE.
    fp=>Forces__Pick(level)

    DO nn=1,SIZE(IndBox_a_p)
       m=IndBox_a_p(nn)
       mass=Atoms(m) % mass
       tfact=0.5_8*ts2/mass

       Atoms(m) % x=Atoms(m) % x+tfact*fp(m) % x+dt*Atoms(m) % vx
       Atoms(m) % y=Atoms(m) % y+tfact*fp(m) % y+dt*Atoms(m) % vy
       Atoms(m) % z=Atoms(m) % z+tfact*fp(m) % z+dt*Atoms(m) % vz
    END DO
  END FUNCTION Atom__Verlet_
  FUNCTION Atom__Correct_(dt,level) RESULT(out)
    INTEGER :: level
    REAL(8) :: dt
    LOGICAL :: out
    TYPE(Force), POINTER :: fp(:)
    INTEGER :: n,nn,m
    REAL(8) :: tfact,mass

    out=.TRUE.
    fp=>Forces__Pick(level)    

    DO nn=1,SIZE(IndBox_a_p)
       m=IndBox_a_p(nn)
       mass=Atoms(m) % mass
       tfact=0.5_8*dt/mass

       Atoms(m) % vx=Atoms(m) % vx+tfact*fp(m) % x
       Atoms(m) % vy=Atoms(m) % vy+tfact*fp(m) % y
       Atoms(m) % vz=Atoms(m) % vz+tfact*fp(m) % z
    END DO
  END FUNCTION Atom__Correct_
  FUNCTION Atom__Write_(Unit,Mode) RESULT(out)
    LOGICAL :: out
    INTEGER :: Mode,Unit
    REAL(8), ALLOCATABLE :: xc(:),yc(:),zc(:)
    INTEGER :: natom,n
    out=.TRUE.

    natom=SIZE(Atoms)
    ALLOCATE(xc(natom),yc(natom),zc(natom))
    xc=Atoms(:) % x
    yc=Atoms(:) % y
    zc=Atoms(:) % z

    CALL PI_Gather_(xc,yc,zc)    

    SELECT CASE(Mode)
    CASE(_SIMPLE_)
       IF(PI_Node_Cart == 0) THEN
          WRITE(90,'(3f12.5,i8)') (xc(n),yc(n),zc(n),n,n=1,natom)
       END IF
    END SELECT
  END FUNCTION Atom__Write_
#define _IN_     IndBox_a_P(:)

!!$---- This module is part of the program oracDD ----*
SUBROUTINE ATOM__KinCompute
  REAL(8), POINTER :: vx(:),vy(:),vz(:),mass(:)
  REAL(8) :: ek_Slt,ek_slv,T_Slv,T_Slt,T_Tot
  INTEGER :: nato_Slt,nato_Slv

  vx=>Atoms(:) % vx
  vy=>Atoms(:) % vy
  vz=>Atoms(:) % vz
  mass=>Atoms(:) % mass

  ek_Slt=SUM(mass(_IN_)*(vx(_IN_)**2+vy(_IN_)**2+vz(_IN_)**2)&
       &,Atoms(_IN_) % Id_Slv == 1) 
  ek_Slv=SUM(mass(_IN_)*(vx(_IN_)**2+vy(_IN_)**2+vz(_IN_)**2)&
       &,Atoms(_IN_) % Id_Slv == 2) 

  ek_Slv=0.5*ek_Slv
  ek_Slv=0.5*ek_Slv
  T_Slv=0.0_8; T_Slt=0.0_8; T_tot=0.0_8
  IF(natom_Slv /= 0) T_Slv=2.0D0*(ek_slv*efact/DBLE(3*natom_Slv))/gascon
  IF(natom_Slt /= 0) T_Slt=2.0D0*(ek_slt*efact/DBLE(3*natom_Slt))/gascon
  T_tot=2.0D0*((ek_slt+ek_Slv)*efact)/(DBLE(3*(natom_Slt+natom_Slv-3))*gascon)

  CALL EN_Kinetic_(ek_Slv,ek_Slt)
  CALL EN_Temp_(T_tot,T_slv,T_Slt)
END SUBROUTINE ATOM__KinCompute

#include "ATOM__vInit.f90"
END MODULE Atom
