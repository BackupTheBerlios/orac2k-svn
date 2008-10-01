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
MODULE MDRun
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Jan 25 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*

#define _PME_ 1
  USE MPI
  USE PME
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f, errmsg_w
  USE Groups
  USE Atom
  USE Inout
  USE PI_Decompose
  USE PI_
  USE PI_Communicate
  USE NeighCells
  USE Ewald
  USE Parallel, ONLY: PA_npx=>npx
  USE Print_Defs
  USE Direct, ONLY: DIR_Forces=>Compute
  USE Forces, ONLY: rcut_o
  IMPLICIT none
  PRIVATE
  PUBLIC MDRun_
CONTAINS
  SUBROUTINE MDRun_
    LOGICAL :: ok
    REAL(8) :: startime,endtime,timea
    INTEGER :: n
    IF(.NOT. Groups_()) STOP
    IF(.NOT. Atom_()) CALL Print_Errors()
    IF(.NOT. Atom__Tpg_()) CALL Print_Errors()
    IF(.NOT. Atom__InitCoords()) CALL Print_Errors()
    IF(.NOT. Groups__InitCoords()) CALL Print_Errors()
!!$    IF(Inout__PDB % unit /= 0) CALL Atom__PDB(Inout__PDB % unit)
!!$    rcut=(/4.7D0, 7.3D0, 9.7D0/)

    IF(PA_npx == 0) THEN
       IF(PI_nprocs /= 0) THEN
          CALL PI__Decomposition_NB(rcut_o)
       END IF
    END IF

    CALL PI__AssignAtomsToCells


    CALL MPI_BARRIER(PI_Comm,ierr)
    startime=MPI_WTIME()

    CALL PI__Shift(3,_PME_)

    IF(Inout__PDB % unit /= 0 .AND. PI_Node_FFTW == 31) CALL Atom__PDB(Inout__PDB % unit, 1)

    CALL DIR_Forces(3)

    endtime=MPI_WTIME()
    timea=endtime-startime
    WRITE(*,*) 'timeo ',timea
    STOP

    IF(Ewald__Param % nx /= 0 .AND. Ewald__Param % ny  /= 0 .AND.&
         & Ewald__Param % nz /= 0) THEN

       CALL PI__Shift(1,_PME_)
       CALL Ewald__Validate

       CALL MPI_BARRIER(PI_Comm_Cart,ierr)
       startime=MPI_WTIME()

       CALL PME_(3)

       CALL MPI_BARRIER(PI_Comm_Cart,ierr)
       endtime=MPI_WTIME()
       timea=endtime-startime
    END IF
    
  END SUBROUTINE MDRun_
END MODULE MDRun
