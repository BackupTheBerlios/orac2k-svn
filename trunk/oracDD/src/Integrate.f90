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
MODULE Integrate
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Wed Dec 17 2008 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*

#include "config.h"

#ifdef HAVE_MPI
  USE mpi
#endif
  USE PI_
  USE Pi_Decompose 
  USE PI_Statistics, ONLY: STAT_Write_it=>Write_it
  USE Print_Defs
  USE PI_Atom, PI_Atom_Update_=>Update_
  USE Forces, FORCES_Zero=>Zero, FORCES_Pick=>Pick
  USE Integrator, ONLY: Integrator_
  USE Atom
  USE Groups
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f, errmsg_w
  USE Direct, ONLY: DIR_Forces=>Compute,DIR_Lists=>Lists
  USE Ewald
  USE PME
  USE PI_Communicate
  USE PI_Collectives
  USE IndBox

  USE IntraAtoms, ONLY: IntraAtoms_,PI__ShiftIntra

  USE Rattle, ONLY: RATTLE__Init_=>Init_,RATTLE__Parameters_&
       &=>Parameters_,RATTLE__Verlet_=>Verlet_,RATTLE__Correct_&
       &=>Correct_ 
  USE Intra, ONLY: Intra_n0_,Intra_n1_
  USE Potential, ONLY: Rattle__Param 
  USE Inout, ONLY: Inout__PDB, Inout__
  USE Run
  USE PI_Statistics, ONLY: PI__Time_It=>Time_It, PI__TTime_It&
       &=>TTime_It, PI__Set_time=>Set_Time
  IMPLICIT none
  PRIVATE
  PUBLIC Integrate_
  REAL(8), SAVE :: dt_h, dt_l, dt_m, dt_n1, dt_n0
  INTEGER, SAVE :: n0_,n1_,m_,l_,h_
  INTEGER, SAVE :: Nshell0=0,Nshell=0
  REAL(8), SAVE :: Time_at_Step=0.0_8
  INTEGER, SAVE :: n0,n1,ma,la,ha
  INTEGER, SAVE :: Nstep_N0=0 !!-- How many N0 steps have been run

  TYPE :: Length
     INTEGER :: Nstep=0
     REAL(8) :: Time=0.0_8
  END TYPE Length
  TYPE(Length) :: Iteration
CONTAINS
  FUNCTION Time_Step() RESULT(out)
    REAL(8) :: out
    out=DBLE(Nstep_N0)*dt_N0
  END FUNCTION Time_Step
  SUBROUTINE Init_

    INTEGER :: n
    n0_=Integrator_ % Mult_Intra(1)
    n1_=Integrator_ % Mult_Intra(2)
    m_ =Integrator_ % Mult_Inter(1)
    l_ =Integrator_ % Mult_Inter(2)
    h_ =Integrator_ % Mult_Inter(3)

    dt_h=Integrator_ % t
    dt_l=dt_h/l_
    dt_m=dt_l/m_
    dt_n1=dt_m/n1_
    dt_n0=dt_n1/n0_
    NShell0=SIZE(Radii)
    NShell=Nshell0+2

!!$--- Find out primary and secondary atoms arrays

    IF(.NOT. IndBox_(Groupa(:) % knwn,Groupa(:) % AtSt,Groupa(:) %&
         & AtEn)) CALL Print_Errors()

    IF(.NOT. Atom__vInit_()) CALL Print_Errors()

  END SUBROUTINE Init_
  SUBROUTINE Integrate_
    INTEGER :: iter,n,Iter_End
    REAL(8) :: startime,endtime,timea,starta,enda
    TYPE(Force), POINTER :: fp_d(:)

    
!!$--- Lennard Jones Arrays. _M_ is not used here it could be zero
!!$--- for what it matters 

    CALL DIR_Forces(_M_,_INIT_)

!!$--- Initialize Ewald parameters for the run

    IF(Ewald__Param % Switch) CALL Ewald__Validate

!!$
!!$--- Assign atoms to simulation cells. 
!!$
    CALL PI__AssignAtomsToCells

!!$--- Set to zero coordinates outside primary cell. Should work
!!$---  without it

    CALL Pi__ZeroSecondary
    IF(Nstep_N0 == 0) CALL Init_

!!$
!!$--- Initialize Rattle
!!$
    IF(.NOT. RATTLE__Init_(SIZE(Atoms))) CALL Print_Errors()

!!$--- Define the Intramolecular environment

    IF(.NOT. IntraAtoms_()) CALL Print_Errors()

!!$--- Reset secondary cell. Need it before shifting

    CALL PI__Shift(NShell,_INIT_)
    CALL PI__Shift(NShell,_EXCHANGE_)

!!$--- Setup the primary and secondary atom cells: IndBox_?_? arrays
!!$--- are created

    IF(.NOT. PI_Atom_()) CALL Print_Errors()

!!$
!!$--- Compute outer neighbor lists
!!$

    IF(.NOT. PI_Atom__Neigh_()) CALL Print_Errors()
    CALL DIR_Lists(Nshell)
    IF(Ewald__Param % Switch .AND. Ewald__Param % nx /= 0 .AND.&
         & Ewald__Param % ny  /= 0 .AND. Ewald__Param % nz /= 0) THEN
       CALL PME_(0)
    END IF
    
    CALL Init_TotalShells(NShell)

    CALL MPI_BARRIER(PI_Comm_cart,ierr)
    startime=MPI_WTIME()

    DO n=1,NShell
       CALL FORCES_Zero(n)
       CALL Forces_(n)
    END DO

    IF(.NOT. RATTLE__Parameters_(Atoms(:) % mass,Atoms(:) % knwn))&
         & CALL Print_Errors()

    endtime=MPI_WTIME()
    timea=endtime-startime
    WRITE(kprint,*) 'First time',PI_Node_Cart,timea

    IF(.NOT. STAT_Write_it()) CALL Print_Errors()

    CALL PI__Set_time
    starta=MPI_WTIME()

    Iteration=Get_RunLength()
    
    startime=MPI_WTIME()
    
    DO Iter=1,Iteration % nstep
       CALL Integrate_Shell(NShell)
    END DO

    endtime=MPI_WTIME()
    timea=endtime-startime

    WRITE(*,*) 'Second Time time',PI_Node_Cart,Timea,timea/Iteration % Time
    WRITE(kprint,*) Iteration % Time

    enda=MPI_WTIME()
    CALL PI__TTime_It(starta,enda) 
    IF(.NOT. STAT_Write_it()) CALL Print_Errors()
  END SUBROUTINE Integrate_
  SUBROUTINE Integrate_Shell(n)
    INTEGER :: n
    SELECT CASE(n)
    CASE(_N0_)
       CALL Integrate_N0
    CASE(_N1_)
       CALL Integrate_N1
    CASE(_M_)
       CALL Integrate_M
    CASE(_L_)
       CALL Integrate_L
    CASE(_H_)
       CALL Integrate_H
    END SELECT
  END SUBROUTINE Integrate_Shell
  FUNCTION Rattle_it(dt,Func) RESULT(out)
    LOGICAL :: out
    REAL(8) :: dt
    INTEGER :: n
    INTERFACE
       FUNCTION Func(dt,xp0a,yp0a,zp0a,vpxa,vpya,vpza) RESULT(out)
         LOGICAL :: out
         REAL(8) :: dt,xp0a(:),yp0a(:),zp0a(:),vpxa(:),vpya(:),vpza(:)
       END FUNCTION Func
    END INTERFACE
    out=.TRUE.
    IF(.NOT. Rattle__Param % switch) RETURN  
    out=Func(dt, Atoms(:) % x,  Atoms(:) % y, Atoms(:) % z, Atoms(:) &
         &% vx,  Atoms(:) % vy, Atoms(:) % vz)
  END FUNCTION Rattle_it
  SUBROUTINE Init_TotalShells(Nstart)
    INTEGER :: Nstart
    INTEGER :: n
    DO n=Nstart,3,-1
       CALL PI__Initialize_Shell(n)
    END DO
  CONTAINS
    SUBROUTINE PI__Initialize_Shell(i_pa)
      INTEGER :: i_pa
      TYPE(Force), POINTER :: fp(:)=>NULL()
      LOGICAL :: do_pme

      CALL PI__ResetSecondary

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
    END SUBROUTINE PI__Initialize_Shell

  END SUBROUTINE Init_TotalShells

#include "Integrate__N0.f90"
#include "Integrate__N1.f90"
#include "Integrate__M.f90"
#include "Integrate__L.f90"
#include "Integrate__H.f90"
#include "Integrate__Forces.f90"
#include "Integrate__Utilities.f90"

END MODULE Integrate
