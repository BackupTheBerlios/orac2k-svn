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
  USE Print_Defs
  USE PI_Atom
  USE Forces, ONLY: Radii, FORCES_Zero=>Zero, FORCES_Pick=>Pick, Force
  USE Integrator, ONLY: Integrator_
  USE Atom
  USE Groups
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f, errmsg_w
  USE Direct, ONLY: DIR_Forces=>Compute
  USE Ewald
  USE PME
  USE PI_Communicate
  USE IndBox
  USE IndIntraBox
  USE PI_IntraMaps
  USE Rattle, ONLY: RATTLE__Init_=>Init_,RATTLE__Parameters_&
       &=>Parameters_,RATTLE__Verlet_=>Verlet_,RATTLE__Correct_&
       &=>Correct_ 
  USE Intra, ONLY: Intra_n0_,Intra_n1_
  USE Potential, ONLY: Rattle__Param 
  USE Inout, ONLY: Inout__PDB, Inout__
  USE Run
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

!!$--- Shit it in the smallest possible shell

    CALL PI__Shift(_M_,_INIT_EXCHANGE_)

!!$--- Find out primary and secondary atoms arrays

    IF(.NOT. IndBox_(Groupa(:) % knwn,Groupa(:) % AtSt,Groupa(:) %&
         & AtEn)) CALL Print_Errors()

    IF(.NOT. Atom__vInit_()) CALL Print_Errors()

  END SUBROUTINE Init_
  SUBROUTINE Integrate_
    INTEGER :: iter,n,Iter_End
    REAL(8) :: startime,endtime,timea

    
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

!!$--- Reset secondary cell. Need it before shifting

    CALL PI__ResetSecondary

!!$--- Shift the outer most atoms around the primary cell

    CALL PI__Shift(NShell,_INIT_EXCHANGE_)

!!$--- Setup the primary and secondary atom cells: IndBox_?_? arrays
!!$--- are created
    IF(.NOT. PI_Atom_()) CALL Print_Errors()

!!$
!!$--- Compute outer neighbor lists
!!$

    IF(.NOT. PI_Atom__Neigh_()) CALL Print_Errors()

    IF(Ewald__Param % Switch .AND. Ewald__Param % nx /= 0 .AND. Ewald__Param % ny  /= 0 .AND.&
         & Ewald__Param % nz /= 0) THEN
       CALL PME_(0)
    END IF

    CALL MPI_BARRIER(PI_Comm_cart,ierr)
    startime=MPI_WTIME()

    DO n=NShell,1,-1
       CALL FORCES_Zero(n)
       CALL Forces_(n,_INIT_EXCHANGE_)
    END DO

    IF(.NOT. RATTLE__Parameters_(Atoms(:) % mass,Atoms(:) % knwn))&
         & CALL Print_Errors()

    endtime=MPI_WTIME()
    timea=endtime-startime
    WRITE(*,*) 'First time',PI_Node_Cart,timea

    CALL MPI_BARRIER(PI_Comm_cart,ierr)
    startime=MPI_WTIME()
    
    Iteration=Get_RunLength()
    

    DO Iter=1,Iteration % nstep
       CALL Integrate_Shell(NShell)
    END DO

    endtime=MPI_WTIME()
    timea=endtime-startime
    WRITE(*,*) 'Second Time time',PI_Node_Cart,timea/Iteration % Time

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

#include "Integrate__N0.f90"
#include "Integrate__N1.f90"
#include "Integrate__M.f90"
#include "Integrate__L.f90"
#include "Integrate__H.f90"
#include "Integrate__Forces.f90"
#include "Integrate__Utilities.f90"

END MODULE Integrate
