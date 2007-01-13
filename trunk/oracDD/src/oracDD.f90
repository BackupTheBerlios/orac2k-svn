!!$/---------------------------------------------------------------------\
!!$                                                                      |
!!$  Copyright (C) 2006-2007 Massimo Marchi <Massimo.Marchi@cea.fr>      |
!!$                                                                      |
!!$      This program is free software;  you  can  redistribute  it      |
!!$      and/or modify it under the terms of the GNU General Public      |
!!$      License version 2 as published  by  the  Free  Software         |
!!$      Foundation;                                                     |
!!$                                                                      |
!!$      This program is distributed in the hope that  it  will  be      |
!!$      useful, but WITHOUT ANY WARRANTY; without even the implied      |
!!$      warranty of MERCHANTABILITY or FITNESS  FOR  A  PARTICULAR      |
!!$      PURPOSE.   See  the  GNU  General  Public License for more      |
!!$      details.                                                        |
!!$                                                                      |
!!$      You should have received a copy of the GNU General  Public      |
!!$      License along with this program; if not, write to the Free      |
!!$      Software Foundation, Inc., 59  Temple  Place,  Suite  330,      |
!!$      Boston, MA  02111-1307  USA                                     |
!!$                                                                      |
!!$\---------------------------------------------------------------------/

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
  
  CALL SimulationBox_
  
!!$
!!$  CALL Topology__SetupTpg
!!$
!!$  CALL Topology__SetupPrm
!!$
!!$  CALL System__Setup

!!$  CALL Run_oracS

  STOP
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*
END PROGRAM OracDD