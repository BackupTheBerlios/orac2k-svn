PROGRAM OracS

!!$***********************************************************************
!!$   Time-stamp: <2006-12-15 09:08:41 marchi>                           *
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

  USE Class_Tree, ONLY: Start
  USE ORAC_INPUT_Mod, ONLY: Options, Store, Modify, Single_String, Input_String
  USE GRAMMAR_Mod, ONLY: Create_list=>Create, Grammar_String, Grammar_Collect,&
       & Read_Grammar
  USE ORAC_PROCESS_Mod, ONLY:  Inputs, Grammar, Process_Commands
  USE TOPOLOGY_Mod, ONLY: Setup_Topology, Setup_Parameters
  IMPLICIT none

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

!!$======================================================================
!!$--  Read Input from command line
!!$======================================================================

  CALL Read_Input

!!$======================================================================
!!$--  Change Input to linked list format: Each environment and command
!!$--  correspond to a node, its commands and keyword being associated
!!$--  to the next node, respectively.
!!$======================================================================

  CALL Collect_Input

!!$======================================================================
!!$--  Read Command Grammar from the GRAMMAR.dat file, then change
!!$--  Grammar to linked list format: Each environment and command
!!$--  corresponds to a node, its commands and keyword being associated
!!$--  to the next node, respectively. To each node a succint description
!!$--  is enclosed
!!$======================================================================

  IF(Read_Grammar()) CALL Collect_Grammar

!!$======================================================================
!!$--  Process the command, store information from the input and reports 
!!$--  of input errors 
!!$======================================================================

  CALL Process_Commands

  CALL Setup_Topology

  CALL Setup_Parameters

!!$  CALL Setup_System

!!$  CALL Run_oracS

  STOP
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*
  CONTAINS
   SUBROUTINE Read_Input
     CALL Options
     CALL Store
     CALL Modify
     CALL Single_String
   END SUBROUTINE Read_Input
   SUBROUTINE Collect_Input
     Grammar_Collect=.FALSE.
     CALL start(Inputs)
     CALL Create_List(Input_String,'')
   END SUBROUTINE Collect_Input
   SUBROUTINE Collect_Grammar
     Grammar_Collect=.TRUE.
     CALL start(Grammar)
     CALL Create_List(Grammar_String,'')
   END SUBROUTINE Collect_Grammar
END PROGRAM OracS
