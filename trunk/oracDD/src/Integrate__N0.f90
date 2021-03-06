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
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Mon Jan  5 2009 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*
SUBROUTINE Integrate_n0
  INTEGER, SAVE :: counter=0
  INTEGER :: Init,q,n0,nn,m
  Integer :: n
  Real(8) :: tfact

  DO n0=1,n0_

!!$
!!$--- Verlet position step
!!$

     __verlet_vp(dt_n0,fpp_n0)

     IF(.NOT. Rattle_it(dt_n0,RATTLE__Verlet_)) CALL Print_Errors()
     
     Init=Pick_Init(_N0_,counter,NShell)

     IF(Init == _INIT_) THEN
        Call GatherGlobal(atoms(:)%x,atoms(:)%y,atoms(:)%z,p0,IndBox_a_p)
        IF(.NOT.  Atom__Convert(_X_TO_XA_,IndBox_a_p)) CALL Print_Errors()
        IF(.NOT. Groups__Update(IndBox_g_p)) CALL Print_Errors()

        CALL PI__ResetSecondary
        CALL PI__Exchange
        CALL PI__AssignAtomsToCells
        
        CALL PI__Shift(NShell,_INIT_)
        CALL PI__Shift(NShell,_EXCHANGE_)

        IF(.NOT. PI_Atom_()) CALL Print_Errors()
        IF(.NOT. PI_Atom__Neigh_()) CALL Print_Errors()
        CALL DIR_Lists(Nshell)
!!$--- Redefine the Intramolecular environment
        IF(.NOT. IntraAtoms_()) CALL Print_Errors()
        CALL Init_TotalShells(NShell)
     END IF


     Call GatherGlobal(atoms(:)%x,atoms(:)%y,atoms(:)%z,p0,IndBox_a_p)
     CALL FORCES_Zero(_N0_)
     CALL Forces_(_N0_)
     Call GatherLocals(fpp_n0,fp_n0(:)%x,fp_n0(:)%y,fp_n0(:)%z,IndBox_a_p)

     __correct_vp(dt_n0,fpp_n0)

     counter=counter+1
     IF(NShell == _N0_ .AND. Init == 0) THEN
        IF(.NOT. RATTLE__Parameters_(Atoms(:) % mass,Atoms(:) %&
             & knwn)) CALL Print_Errors()
     END IF
     nstep_n0=nstep_n0+1
  END DO
END SUBROUTINE Integrate_n0
