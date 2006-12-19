PROGRAM OracS

!!$***********************************************************************
!!$   Time-stamp: <2006-12-19 14:31:11 marchi>                           *
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
  USE INPUT_Mod, ONLY: Input__Options=>Options, Input__Store=>Store&
       &, Input__Modify=>Modify, Input__Single_String=>Single_String&
       &, Input__Input_String=>Input_String
  USE GRAMMAR_Mod, ONLY: Grammar__Create=>Create, Grammar_String, Grammar_Collect,&
       & Read_Grammar
  USE PROCESS_Mod, ONLY:  Inputs, Grammar, Process_Commands
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
     CALL Input__Options
     CALL Input__Store
     CALL Input__Modify
     CALL Input__Single_String
   END SUBROUTINE Read_Input
   SUBROUTINE Collect_Input
     Grammar_Collect=.FALSE.
     CALL start(Inputs)
     CALL Grammar__Create(Input__Input_String,'')
   END SUBROUTINE Collect_Input
   SUBROUTINE Collect_Grammar
     Grammar_Collect=.TRUE.
     CALL start(Grammar)
     CALL Grammar__Create(Grammar_String,'')
   END SUBROUTINE Collect_Grammar
END PROGRAM OracS
