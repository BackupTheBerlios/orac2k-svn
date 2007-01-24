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

PROGRAM OracDD
!!$***********************************************************************
!!$   Time-stamp: <2007-01-13 01:33:22 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Nov 17 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program oracDD ----*


!!$======================== DECLARATIONS ================================*

  USE Errors,ONLY: Print_Errors, Print_Warnings
  USE Tree, ONLY: Tree__Start
  USE Inputs
  USE Grammars
  USE Process
  USE Units
  USE Tops
  USE AtomCnt
  USE SystemTpg
  USE SystemPrm
  USE SimulationBox
  USE Parameters
  USE SecondarySeq
  USE IndSequence
  USE Groups
  
!!$  USE PROCESS_Mod, ONLY:  Inputs, Grammar, Process__Construe=>Construe
!!$  USE TOPOLOGY_Mod, ONLY: Topology__SetupTpg=>SetupTpg, &
!!$       &Topology__SetupPrm=>SetupPrm
!!$  USE SYSTEM_Mod, ONLY: System__Setup=>Setup
  IMPLICIT none
  REAL(8) :: Time_Begin,Time_End
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

  CALL Units_

!!$======================================================================
!!$--  Read Input from command line
!!$======================================================================

  CALL Inputs_

!!$======================================================================
!!$--  Read Command Grammar from the GRAMMAR.dat file, then change
!!$--  Grammar to linked list format: Each environment and command
!!$--  corresponds to a node, its commands and keyword being associated
!!$--  to the next node, respectively. To each node a succint description
!!$--  is enclosed
!!$======================================================================

  CALL Grammars_(Inputs__String)

!!$======================================================================
!!$--  Process the command, store information from the input and reports 
!!$--  of input errors 
!!$======================================================================

  CALL Process_

  CALL BuildSystem

  IF(.NOT. Groups_()) STOP
!!$
!!$  CALL Topology__SetupTpg
!!$
!!$  CALL Topology__SetupPrm
!!$
!!$  CALL System__Setup

!!$  CALL Run_oracS

  STOP
CONTAINS
  SUBROUTINE BuildSystem
    IF(Called_Tpg .AND. Called_Prm) THEN
       CALL Tops_       
       CALL AtomCnt_       
       CALL SystemTpg_       
       CALL SystemPrm_       
       CALL SimulationBox_
       CALL AtomCnt__Update(nunits_Slv)       
       CALL SystemTpg__Update(nunits_Slv)       
       CALL SystemPrm__Update
       CALL SecondarySeq__AddSlv(nunits_Slv)
       CALL IndSequence__Update
       IF(Called_Binary) THEN
          CALL AtomCnt__Write
          CALL SystemTpg__Write
          CALL SystemPrm__Write
          CALL SecondarySeq__Write(kbinary)
          CALL IndSequence__Write
       END IF
    ELSE
       IF(.NOT. AtomCnt__Read()) Call Print_Errors()
       IF(.NOT. SystemTpg__Read()) Call Print_Errors()
       CALL SystemPrm__Read
       IF(.NOT. SecondarySeq__Read(kbinary)) Call Print_Errors()
       IF(.NOT. IndSequence__Read()) Call Print_Errors()
    END IF
    CALL Print_Warnings()
  END SUBROUTINE BuildSystem
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*
END PROGRAM OracDD
