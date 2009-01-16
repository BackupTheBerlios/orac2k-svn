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
MODULE PI_Communicate
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Jan 26 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*

#include "config.h"

#ifdef HAVE_MPI
  USE mpi
#endif
  USE PI_Shift, ONLY: SHIFT__iShift_Init=>iShift_Init, SHIFT__Setup&
       &=>Setup, SHIFT__Buff_Shift=>Buff_Shift 
  USE PI_Fold, ONLY: FOLD__iFold_Init=>iFold_Init, FOLD__Setup&
       &=>Setup, FOLD__Buff_Fold=>Buff_Fold 
  USE PI_ShiftIntra, ONLY: SHIFTINTRA__iShift_Init=>iShift_Init, SHIFTINTRA__Setup&
       &=>Setup, SHIFTINTRA__Buff_Shift=>Buff_Shift 
  USE PI_FoldIntra, ONLY: FOLDINTRA__iFold_Init=>iFold_Init, FOLDINTRA__Setup&
       &=>Setup, FOLDINTRA__Buff_Fold=>Buff_Fold 

  USE Integrator, ONLY: Integrator_
  USE Ewald
  USE Errors, ONLY: Add_Errors=>Add, errmsg_f,Print_Errors
  USE PI_Exchange, ONLY: EX_Exchange=>Exchange_
  USE PI_Cutoffs
  USE Geometry
  USE Forces, ONLY: Force, Radii
  USE UNITS
  USE PI_
  USE Groups
  USE Atom
  USE Cell
  USE IndBox
  USE PI_Statistics, ONLY: PI__Write_Stats=>Write_It, PI__Time_It&
       &=>Time_It, PI__Sample_Exchange=>Sample_Exchange, PI__Add_Calls&
       &=>Add_Calls

  IMPLICIT none
  PRIVATE
  PUBLIC PI__Shift,PI__Fold_F, PI__ZeroSecondary,PI__ZeroPrimary,&
       & PI__ShiftIntra, PI__ResetSecondary, PI__ResetSecondaryP&
       &,PI__FoldIntra, PI__Exchange,PI__Initialize_Shell
  TYPE :: Communicate__Chain
     INTEGER :: i,j,k
     INTEGER :: p
  END TYPE Communicate__Chain

  REAL(8), ALLOCATABLE, SAVE :: Buffer(:)
  REAL(8), SAVE :: vp(3),vd(3)
  INTEGER, SAVE :: npx,npy,npz

  TYPE(Communicate__Chain), ALLOCATABLE,SAVE :: Chain_xyz(:)
  INTEGER, ALLOCATABLE,SAVE :: Head_xyz(:)
  REAL(8), SAVE :: Thick(3)
  LOGICAL, SAVE :: do_pme
  REAL(8), PARAMETER :: one=1.0D0,two=2.0D0,half=0.5D0
  INTEGER, SAVE :: Calls=0
  REAL(8), SAVE :: startime,endtime,startime0,endtime0
CONTAINS
  SUBROUTINE PI__Initialize_Shell(i_pa)
    INTEGER :: i_pa
    TYPE(Force), POINTER :: fp(:)=>NULL()

    CALL PI__ResetSecondary
    IF(i_pa >= 3) THEN
       CALL PI__Shift(i_pa,_INIT_)

       do_pme=(i_pa-2 == Integrator_ % Ewald_Shell) .AND. Ewald__Param %&
            & Switch .AND. ( Ewald__Param % nx /= 0 .AND.&
            & Ewald__Param % ny  /= 0 .AND. Ewald__Param % nz /= 0)

       IF(do_pme) THEN
          IF(.NOT. IndBoxP_(Groupa(:) % knwn,Groupa(:) % AtSt,Groupa(:) %&
               & AtEn)) CALL Print_Errors() 
       END IF

       CALL PI__Fold_F(fp,i_pa,_INIT_)

       IF(do_pme) CALL PI__ResetSecondaryP
    ELSE
       CALL PI__ShiftIntra(i_pa,_INIT_)
       CALL PI__FoldIntra(fp,i_pa,_INIT_)
    END IF

  END SUBROUTINE PI__Initialize_Shell
  SUBROUTINE PI__Shift(i_pa,init)
    INTEGER :: init,i_pa
    INTEGER, SAVE :: ShiftTime,source, dest
    INTEGER :: ox,oy,oz,numcell,mpe,mp,m,n
    INTEGER :: nmin,i,j,k,np,AtSt,AtEn,l,q
    INTEGER :: iv(3),Axis,i_p
    LOGICAL :: pme

    i_p=i_pa-2
    do_pme=(i_p == Integrator_ % Ewald_Shell) .AND. Ewald__Param %&
         & Switch .AND. ( Ewald__Param % nx /= 0 .AND.&
         & Ewald__Param % ny  /= 0 .AND. Ewald__Param % nz /= 0)

    npx=PI_npx
    npy=PI_npy
    npz=PI_npz
    CALL Thickness(i_p)
    
!!$
!!$ --- Atoms to send
!!$
    

    iv(1)=npx
    iv(2)=npy
    iv(3)=npz

    IF(PI_Nprocs == 1) RETURN

    CALL SHIFT__Setup(do_pme)

    SELECT CASE(init)
    CASE(_INIT_)
       CALL Shift_it(SHIFT__iShift_init)
    CASE(_EXCHANGE_)
       CALL Shift_it(SHIFT__Buff_Shift)
    END SELECT

!!$
!!$--- Add to Calls counter: One data exchange has occurred
!!$
    CALL MPI_BARRIER(PI_Comm_Cart,ierr)
    CALL PI__Add_Calls
  CONTAINS
    SUBROUTINE Shift_it(Routine)
      INTERFACE
         SUBROUTINE Routine(i_p,Axis,Dir,scnd_half)
           INTEGER, OPTIONAL :: scnd_half
           INTEGER :: Axis,Dir,i_p
         END SUBROUTINE Routine
      END INTERFACE
      IF(do_pme) THEN
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
      ELSE
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
      END IF
    END SUBROUTINE Shift_it
  END SUBROUTINE PI__Shift
!!$
!!$---- Fold forces
!!$
  SUBROUTINE PI__Fold_F(fp,i_pa,init)
    TYPE(Force) :: fp(:)
    INTEGER :: i_pa,init
    INTEGER, SAVE :: ShiftTime,source, dest
    INTEGER :: ox,oy,oz,numcell,mpe,mp,m,n
    INTEGER :: nmin,i,j,k,np,AtSt,AtEn,l,q,nn,grp_no
    INTEGER :: iv(3),Axis,i_p
    INTEGER, SAVE :: MyownCalls=0
    i_p=i_pa-2
    

    CALL Thickness(i_p)

    npx=PI_npx
    npy=PI_npy
    npz=PI_npz
    
!!$
!!$ --- Fold Forces
!!$

    iv(1)=npx
    iv(2)=npy
    iv(3)=npz


    IF(PI_Nprocs == 1) RETURN

    do_pme=(i_p == Integrator_ % Ewald_Shell) .AND. Ewald__Param %&
         & Switch .AND. ( Ewald__Param % nx /= 0 .AND.&
         & Ewald__Param % ny  /= 0 .AND. Ewald__Param % nz /= 0)

    CALL FOLD__Setup(do_pme)

    SELECT CASE(init)
    CASE(_INIT_)
       CALL Fold_it(FOLD__iFold_init)
    CASE(_EXCHANGE_)
       CALL Fold_it(FOLD__Buff_Fold)
    END SELECT

    CALL MPI_BARRIER(PI_Comm_Cart,ierr)
    CALL PI__Add_Calls
    MyOwnCalls=MyOwnCalls+1
  CONTAINS
    SUBROUTINE Fold_it(Routine)
      INTERFACE
         SUBROUTINE Routine(fp,i_p,Axis,Dir)
           USE Forces, ONLY: Force
           TYPE(Force) :: fp(:)
           INTEGER :: Axis,Dir,i_p
         END SUBROUTINE Routine
      END INTERFACE
      IF(do_pme) THEN
         IF(iv(1) == 1 .AND. iv(2) /= 1) THEN
            CALL Routine(fp,i_p,3,_PLUS_)
            CALL Routine(fp,i_p,3,_MINUS_)
            CALL Routine(fp,i_p,2,_PLUS_)
            CALL Routine(fp,i_p,2,_MINUS_)
         ELSE IF(iv(1) == 1 .AND. iv(2) == 1) THEN
            CALL Routine(fp,i_p,3,_PLUS_)
            CALL Routine(fp,i_p,3,_MINUS_)
         ELSE
            CALL Routine(fp,i_p,1,_PLUS_)
            CALL Routine(fp,i_p,1,_MINUS_)
            CALL Routine(fp,i_p,2,_PLUS_)
            CALL Routine(fp,i_p,2,_MINUS_)
            CALL Routine(fp,i_p,3,_MINUS_)
            CALL Routine(fp,i_p,3,_PLUS_)
         END IF
      ELSE
         IF(iv(1) == 1 .AND. iv(2) /= 1) THEN
            CALL Routine(fp,i_p,3,_PLUS_)
            CALL Routine(fp,i_p,3,_MINUS_)
            CALL Routine(fp,i_p,2,_MINUS_)
         ELSE IF(iv(1) == 1 .AND. iv(2) == 1) THEN
            CALL Routine(fp,i_p,3,_MINUS_)
         ELSE
            CALL Routine(fp,i_p,3,_PLUS_)
            CALL Routine(fp,i_p,2,_PLUS_)
            CALL Routine(fp,i_p,3,_MINUS_)
            CALL Routine(fp,i_p,2,_MINUS_)
            CALL Routine(fp,i_p,1,_MINUS_)
         END IF
      END IF
    END SUBROUTINE Fold_it
  END SUBROUTINE PI__Fold_F
  SUBROUTINE PI__ShiftIntra(i_p,init)
    INTEGER :: init,i_p

    INTEGER, SAVE :: ShiftTime,source, dest
    INTEGER :: ox,oy,oz,numcell,mpe,mp,m,n
    INTEGER :: nmin,i,j,k,np,AtSt,AtEn,l,q
    INTEGER :: iv(3),Axis
    INTEGER :: ng
    npx=PI_npx
    npy=PI_npy
    npz=PI_npz

    CALL Thickness(i_p)

!!$
!!$ --- Atoms to send
!!$
    

    iv(1)=npx
    iv(2)=npy
    iv(3)=npz

    IF(PI_Nprocs == 1) RETURN

    CALL SHIFTINTRA__Setup

    SELECT CASE(init)
    CASE(_INIT_EXCHANGE_)
       CALL Shift_it(SHIFTINTRA__iShift_init)
    CASE(_EXCHANGE_ONLY_)
       CALL Shift_it(SHIFTINTRA__Buff_Shift)
    END SELECT


!!$
!!$--- Add to Calls counter: One data exchange has occurred
!!$
    CALL MPI_BARRIER(PI_Comm_Cart,ierr)

    CALL PI__Add_Calls
  CONTAINS
    SUBROUTINE Shift_it(Routine)
      INTERFACE
         SUBROUTINE Routine(i_p,Axis,Dir,scnd_half)
           INTEGER, OPTIONAL :: scnd_half
           INTEGER :: Axis,Dir,i_p
         END SUBROUTINE Routine
      END INTERFACE
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
    END SUBROUTINE Shift_it
  END SUBROUTINE PI__ShiftIntra

!!$
!!$---- Fold forces
!!$
  SUBROUTINE PI__FoldIntra(fp,i_p,init)
    TYPE(Force) :: fp(:)
    INTEGER :: i_p,init
    INTEGER, SAVE :: ShiftTime,source, dest
    INTEGER :: ox,oy,oz,numcell,mpe,mp,m,n
    INTEGER :: nmin,i,j,k,np,AtSt,AtEn,l,q,nn,grp_no
    INTEGER :: iv(3),Axis

    CALL Thickness(i_p)

    npx=PI_npx
    npy=PI_npy
    npz=PI_npz
    
!!$
!!$ --- Fold Forces
!!$

    iv(1)=npx
    iv(2)=npy
    iv(3)=npz


    IF(PI_Nprocs == 1) RETURN

    CALL FOLDINTRA__Setup

    SELECT CASE(init)
    CASE(_INIT_EXCHANGE_)
       CALL Fold_it(FOLDINTRA__iFold_init)
    CASE(_EXCHANGE_ONLY_)
       CALL Fold_it(FOLDINTRA__Buff_Fold)
    END SELECT

    CALL MPI_BARRIER(PI_Comm_Cart,ierr)
    CALL PI__Add_Calls
  CONTAINS
    SUBROUTINE Fold_it(Routine)
      INTERFACE
         SUBROUTINE Routine(fp,i_p,Axis,Dir)
           USE Forces, ONLY: Force
           TYPE(Force) :: fp(:)
           INTEGER :: Axis,Dir,i_p
         END SUBROUTINE Routine
      END INTERFACE
      IF(iv(1) == 1 .AND. iv(2) /= 1) THEN
         CALL Routine(fp,i_p,3,_PLUS_)
         CALL Routine(fp,i_p,3,_MINUS_)
         CALL Routine(fp,i_p,2,_MINUS_)
      ELSE IF(iv(1) == 1 .AND. iv(2) == 1) THEN
         CALL Routine(fp,i_p,3,_MINUS_)
      ELSE
         CALL Routine(fp,i_p,3,_PLUS_)
         CALL Routine(fp,i_p,2,_PLUS_)
         CALL Routine(fp,i_p,3,_MINUS_)
         CALL Routine(fp,i_p,2,_MINUS_)
         CALL Routine(fp,i_p,1,_MINUS_)
      END IF
    END SUBROUTINE Fold_it
  END SUBROUTINE PI__FoldIntra
  SUBROUTINE PI__Exchange
    INTEGER, SAVE :: ShiftTime,source, dest
    INTEGER :: ox,oy,oz,numcell,mpe,mp,m,n
    INTEGER :: nmin,i,j,k,np,AtSt,AtEn,l,q
    INTEGER :: iv(3),Axis,i_p

!!$
!!$ --- Atoms to send
!!$

    npx=PI_npx
    npy=PI_npy
    npz=PI_npz

    iv(1)=npx
    iv(2)=npy
    iv(3)=npz

    IF(PI_Nprocs == 1) RETURN
    CALL Exchange_it

  CONTAINS
    SUBROUTINE Exchange_it
      IF(iv(1) == 1 .AND. iv(2) /= 1) THEN
         CALL EX_Exchange(2,_PLUS_)
         CALL EX_Exchange(3,_PLUS_)
         CALL EX_Exchange(2,_MINUS_)
         CALL EX_Exchange(3,_MINUS_)
      ELSE IF(iv(1) == 1 .AND. iv(2) == 1) THEN
         CALL EX_Exchange(3,_PLUS_)
         CALL EX_Exchange(3,_MINUS_)
      ELSE
         CALL EX_Exchange(1,_MINUS_)
         CALL EX_Exchange(2,_PLUS_)
         CALL EX_Exchange(3,_PLUS_)
         CALL EX_Exchange(1,_PLUS_)
         CALL EX_Exchange(2,_MINUS_)
         CALL EX_Exchange(3,_MINUS_)
      END IF
    END SUBROUTINE Exchange_it
  END SUBROUTINE PI__Exchange

  SUBROUTINE PI__ResetSecondary
    WHERE(Groupa(IndBox_g_t(:)) % knwn /= 1) Groupa(IndBox_g_t(:)) % knwn = 0
    WHERE(Atoms(IndBox_a_t(:)) % knwn /= 1) Atoms(IndBox_a_t(:)) % knwn = 0
  END SUBROUTINE PI__ResetSecondary
  SUBROUTINE PI__ResetSecondaryP
    WHERE(Groupa(IndBoxP_g_t(:)) % knwn /= 1) Groupa(IndBoxP_g_t(:)) % knwn = 0
    WHERE(Atoms(IndBoxP_a_t(:)) % knwn /= 1) Atoms(IndBoxP_a_t(:)) % knwn = 0
  END SUBROUTINE PI__ResetSecondaryP
  SUBROUTINE PI__ZeroSecondary
    INTEGER :: n,AtSt,AtEn,mm

    DO n=1,SIZE(Groupa)
       IF(Groupa(n) % knwn /= 1 ) THEN
          Groupa(n) % knwn = 0
          Groupa(n) % xa=0.0D0
          Groupa(n) % ya=0.0D0
          Groupa(n) % za=0.0D0
          Groupa(n) % x=0.0D0
          Groupa(n) % y=0.0D0
          Groupa(n) % z=0.0D0
          AtSt=Groupa(n) % AtSt
          AtEn=Groupa(n) % AtEn
          DO mm=AtSt,AtEn
             Atoms(mm) % knwn=0
             Atoms(mm) % x = 0.0D0
             Atoms(mm) % y = 0.0D0
             Atoms(mm) % z = 0.0D0
             Atoms(mm) % xa = 0.0D0
             Atoms(mm) % ya = 0.0D0
             Atoms(mm) % za = 0.0D0
          END DO
       END IF          
    END DO
    
  END SUBROUTINE PI__ZeroSecondary
  SUBROUTINE PI__ZeroPrimary
    INTEGER :: n,AtSt,AtEn,mm
    DO n=1,SIZE(Groupa)
       IF(Groupa(n) % knwn /= 2 ) THEN
          Groupa(n) % knwn = 0
          Groupa(n) % xa=0.0D0
          Groupa(n) % ya=0.0D0
          Groupa(n) % za=0.0D0
          Groupa(n) % x=0.0D0
          Groupa(n) % y=0.0D0
          Groupa(n) % z=0.0D0
          AtSt=Groupa(n) % AtSt
          AtEn=Groupa(n) % AtEn
          DO mm=AtSt,AtEn
             Atoms(mm) % x = 0.0D0
             Atoms(mm) % y = 0.0D0
             Atoms(mm) % z = 0.0D0
             Atoms(mm) % xa = 0.0D0
             Atoms(mm) % ya = 0.0D0
             Atoms(mm) % za = 0.0D0
          END DO
       END IF          
    END DO
    
  END SUBROUTINE PI__ZeroPrimary
END MODULE PI_Communicate
