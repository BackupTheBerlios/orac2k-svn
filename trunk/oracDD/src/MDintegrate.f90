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

  USE PI_
  USE Pi_Decompose 
  USE Print_Defs
  USE PI_Atom
  USE Forces, ONLY: Radii
  USE Integrator, ONLY: Integrator_
  USE Atom
  USE Groups
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f, errmsg_w
  USE Direct, ONLY: DIR_Forces=>Compute
  USE PI_Communicate
  USE IndBox
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
    WRITE(*,*) 'ulla ',COUNT(Atoms(:) % knwn == 2),COUNT(Atoms(:) % knwn == 1)
    WRITE(*,*) 'ulla ',COUNT(Groupa(:) % knwn == 2),COUNT(Groupa(:) % knwn == 1)

  END SUBROUTINE Init_
  SUBROUTINE Integrate
    INTEGER :: n

!!$
!!$--- Assign atoms to simulation cells. 
!!$
    CALL PI__AssignAtomsToCells

!!$--- Set to zero coordinates outside primary cell. Should work
!!$---  without it

    CALL Pi__ZeroSecondary
    IF(nstep == 0) CALL Init_

    CALL PI__ResetSecondary

    CALL PI__Shift(NShell,_INIT_EXCHANGE_)
    WRITE(*,*) 'ulla ',NShell,COUNT(Groupa(:) % knwn == 2),COUNT(Groupa(:) % knwn == 1)

    IF(.NOT. PI_Atom_()) CALL Print_Errors()
    IF(.NOT. PI_Atom__Neigh_()) CALL Print_Errors()

    CALL DIR_Forces(Nshell,_INIT_)

    DO n=NShell,3,-1
       CALL PI__ResetSecondary
       CALL PI__Shift(n,_INIT_EXCHANGE_)
       CALL DIR_Forces(n)
    END DO


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
END MODULE MDintegrate
