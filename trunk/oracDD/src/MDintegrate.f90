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
MODULE MDintegrate
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
  USE PI_IntraMaps, ONLY: IntraMaps_n0_, IntraMaps_n1_
  USE Intra, ONLY: Intra_n0_,Intra_n1_
  IMPLICIT none
  PRIVATE
  PUBLIC Integrate
  REAL(8), SAVE :: dt_h, dt_l, dt_m, dt_n1, dt_n0
  INTEGER, SAVE :: n0_,n1_,m_,l_,h_
  INTEGER, SAVE :: nstep=0,Nshell0=0,Nshell=0
  REAL(8), SAVE :: Time_at_Step=0.0_8
  
CONTAINS
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
  SUBROUTINE Integrate
    INTEGER :: n
    REAL(8) :: startime,endtime,timea

    
!!$--- Lennard Jones Arrays

    CALL DIR_Forces(_M_,_INIT_)

!!$--- Initialize Ewald parameters for the run

    CALL Ewald__Validate

!!$
!!$--- Assign atoms to simulation cells. 
!!$
    CALL PI__AssignAtomsToCells

!!$--- Set to zero coordinates outside primary cell. Should work
!!$---  without it

    CALL Pi__ZeroSecondary
    IF(nstep == 0) CALL Init_

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

    IF(Ewald__Param % nx /= 0 .AND. Ewald__Param % ny  /= 0 .AND.&
         & Ewald__Param % nz /= 0) THEN
       CALL PME_(0)
    END IF

    CALL MPI_BARRIER(PI_Comm_cart,ierr)
    startime=MPI_WTIME()

    DO n=NShell,3,-1
       CALL FORCES_Zero(n)
       CALL InterForces_(n,_INIT_EXCHANGE_)
    END DO

!!$
!!$--- Pick all groups within the range of the bonded interactions in
!!$--- the primary cell. Used by ShiftIntra
!!$

    CALL IntraMaps_n0_
    CALL IntraMaps_n1_

!!$--- Shift atoms
    CALL PI__ResetSecondary
    CALL PI__ShiftIntra(_N0_,_INIT_EXCHANGE_)
!!$--- Gets all n0 interactions beloging to the primary cell
    IF(.NOT. IndIntraBox_n0_()) CALL Print_Errors()
    CALL Intra_n0_(_INIT_EXCHANGE_)
    
!!$--- Shift atoms
    CALL PI__ResetSecondary
    CALL PI__ShiftIntra(_N1_,_INIT_EXCHANGE_)
!!$--- Gets all n1 interactions beloging to the primary cell
    IF(.NOT. IndIntraBox_n1_()) CALL Print_Errors()
    CALL Intra_n1_(_INIT_EXCHANGE_)

    endtime=MPI_WTIME()
    timea=endtime-startime
    WRITE(*,*) 'First time',PI_Node_Cart,timea


!!$    CALL FORCES_Memory(SIZE(Atoms))
!!$    CALL PI__AssignAtomsToCells
!!$
!!$    CALL PI__Shift(1,_INIT_EXCHANGE_)
!!$    CALL DIR_Forces(1,_INIT_)
!!$
!!$
!!$    CALL Integrate_m
!!$
!!$    Time_at_Step=Update_Step()
!!$  CONTAINS
  END SUBROUTINE Integrate

!!$  SUBROUTINE Integrate_n0
!!$    DO n0=1,n0_
!!$       IF(.NOT. Atom__Correct_(dt_n0,_N0_)) CALL Print_Errors()
!!$       
!!$       IF(.NOT. Atom__Verlet_(dt_n0,_N0_)) CALL Print_Errors()
!!$       
!!$       CALL GetForces(_N0_)
!!$       
!!$       IF(.NOT. Atom__Correct_(dt_n0,_N0_)) CALL Print_Errors()
!!$    END DO
!!$  END SUBROUTINE Integrate_n0
!!$
!!$  SUBROUTINE Integrate_n1
!!$    DO n1=1,n1_
!!$       IF(.NOT. Atom__Correct_(dt_n1,_N1_)) CALL Print_Errors()
!!$
!!$       CALL Integrate_n0
!!$
!!$       CALL GetForces(_N1_)
!!$       
!!$       IF(.NOT. Atom__Correct_(dt_n1,_N1_)) CALL Print_Errors()
!!$    END DO
!!$  END SUBROUTINE Integrate_n1
!!$
!!$  SUBROUTINE Integrate_m
!!$    DO m=1,m_
!!$       IF(.NOT. Atom__Correct_(dt_m,_M_)) CALL Print_Errors()
!!$
!!$       CALL Initialize_Intra
!!$
!!$       CALL Integrate_n1
!!$
!!$       CALL GetForces(_M_)
!!$       
!!$       IF(.NOT. Atom__Correct_(dt_m,_M_)) CALL Print_Errors()
!!$    END DO
!!$  END SUBROUTINE Integrate_m
!!$
!!$  SUBROUTINE GetForces(Init)
!!$    
!!$    SELECT CASE(Init)
!!$    CASE(_N0_)
!!$       CALL PI__ShiftIntra(_N0_,Track_Intra)
!!$       CALL Intra_n0_(Track_Intra)
!!$    CASE(_N1_)
!!$       CALL PI__ShiftIntra(_N1_,Track_Intra)
!!$       CALL Intra_n1_(Track_Intra)
!!$    CASE(_M_)
!!$
!!$    END SELECT
!!$
!!$  END SUBROUTINE GetForces
!!$  
!!$  SUBROUTINE Initialize_Intra
!!$    
!!$    CALL PI__ZeroSecondary
!!$    CALL PI__ResetSecondary
!!$    
!!$    CALL IntraMaps_n0_
!!$    CALL PI__ShiftIntra(1,_INIT_EXCHANGE_)
!!$    IF(.NOT. IndIntraBox_n0_()) CALL Print_Errors()
!!$    
!!$    CALL PI__ZeroSecondary
!!$    CALL PI__ResetSecondary
!!$    
!!$    CALL IntraMaps_n1_
!!$    CALL PI__ShiftIntra(2,_INIT_EXCHANGE_)
!!$    IF(.NOT. IndIntraBox_n1_()) CALL Print_Errors()
!!$    Track_Intra=0
!!$  END SUBROUTINE Initialize_Intra
!!$  
!!$  FUNCTION Update_Step RESULT(out)
!!$    REAL(8) :: out
!!$    nstep=nstep+1
!!$    out=nstep*Integrator_ % t
!!$  END FUNCTION Update_Step
  SUBROUTINE InterForces_(n,Flag)
    INTEGER :: n,Flag
    LOGICAL :: pme
    TYPE(Force), POINTER :: fp(:)
    pme= (n-2 == Integrator_ % Ewald_Shell)
    CALL PI__ResetSecondary
    IF(pme) THEN
       IF(Ewald__Param % nx /= 0 .AND. Ewald__Param % ny  /= 0 .AND.&
            & Ewald__Param % nz /= 0) THEN
          CALL PI__Shift(n,Flag,_PME_)
       END IF
    ELSE
       CALL PI__Shift(n,Flag)
    END IF
    CALL DIR_Forces(n)
    IF(pme) THEN
       CALL PME_(n)
    END IF
!!$
!!$--- Fold forces contributions to atoms inside the cell
!!$
    fp=>FORCES_Pick(n)
    CALL PI__Fold_F(fp,n,Flag)

  END SUBROUTINE InterForces_
END MODULE MDintegrate
