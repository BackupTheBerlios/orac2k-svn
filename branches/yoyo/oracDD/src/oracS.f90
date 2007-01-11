PROGRAM OracS

!!$***********************************************************************
!!$   Time-stamp: <2007-01-11 18:08:56 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Nov 17 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program  ----*


!!$======================== DECLARATIONS ================================*

  USE Tree, ONLY: Tree__Start
  USE Inputs
  USE Grammars
  USE Process
  USE Units
  USE Tops
  USE AtomCnt
  USE SystemTpg
  USE SystemPrm

  
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
  CALL Tops_

  CALL AtomCnt_

  CALL SystemTpg_

  CALL SystemPrm_

  
!!$
!!$  CALL Topology__SetupTpg
!!$
!!$  CALL Topology__SetupPrm
!!$
!!$  CALL System__Setup

!!$  CALL Run_oracS

  STOP
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*
END PROGRAM OracS
