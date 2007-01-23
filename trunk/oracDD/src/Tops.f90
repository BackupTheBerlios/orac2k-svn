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
MODULE Tops

!!$***********************************************************************
!!$   Time-stamp: <2007-01-12 12:27:25 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Tue Nov 14 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*
  

  USE Solvent
  USE SecondarySeq
  USE Types
  USE Myparse
  USE Strings, ONLY: MY_Fam, My_Fxm,ST_Concat
  USE Hash_Tops
  USE Resid, ONLY: Resids, Residue=>Topology
  USE Node
  USE STRPAK
  USE Parameters_Globals
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f
  PRIVATE
  PUBLIC Tops_, Tops__Type, App_Char
  TYPE :: Tops__Type
     CHARACTER(len=max_char) :: FField=' ',Type=' ',Residue=' '
     CHARACTER(len=max_char) , DIMENSION(:,:), ALLOCATABLE :: bonds,imph&
          &,acc,acc_,don,don_,dele,ends
     CHARACTER(len=max_char), DIMENSION(:,:), ALLOCATABLE :: mass
     TYPE(list), DIMENSION(:), ALLOCATABLE :: group
  END TYPE Tops__Type
  TYPE(Tops__type), DIMENSION(:), ALLOCATABLE, SAVE, TARGET :: Res_char&
       &,Add_Char
  TYPE(Tops__type), DIMENSION(:), ALLOCATABLE, SAVE :: App_Char
  TYPE(Hash_Tops__Type), SAVE :: iKeys
  CHARACTER(len=max_char), DIMENSION(6), PARAMETER :: Oracs=(/'bonds',&
       &'imph ','don  ','acc  ','terma','backb'/)
  LOGICAL, ALLOCATABLE, SAVE :: ok_Residue(:)
  TYPE(Tops__Type), POINTER :: Resi

  CHARACTER(len=max_char), DIMENSION(:,:), ALLOCATABLE :: share
  TYPE(List), DIMENSION(:), ALLOCATABLE :: shareg
  TYPE(Resids), POINTER :: Ri
CONTAINS
  SUBROUTINE Tops_
    IMPLICIT NONE 
    INTEGER, SAVE :: i_L,n,m,p,o
    REAL(8) :: Time_begin, Time_end
    CALL Hash_Tops_

    IF(Solvent__Param % Build .AND. &
         &Solvent__Param % Added /=  0) &
         &CALL SecondarySeq__AddSlv(Solvent__Param % Added)

    ALLOCATE(Res_Char(SIZE(Residue)))
    ALLOCATE(ok_Residue(SIZE(Residue)))
    ok_residue=.FALSE.
    DO i_l=1,SIZE(Residue)
       Ri=>Residue(i_L)
       IF(TRIM(Ri % Type) == 'mass') THEN
          ok_Residue(i_L)=.TRUE.
          CYCLE
       END IF
       DO n=1,SIZE(Secondary)
          IF(.NOT. ALLOCATED(Secondary(n) % line)) CYCLE
          DO m=1,SIZE(Secondary (n) % line)
             IF(TRIM(Secondary (n) % line (m)) == TRIM(Ri% Type)) THEN
                ok_Residue(i_L)=.TRUE.
                EXIT
             END IF
          END DO
       END DO
       IF(.NOT. ALLOCATED(Patches)) CYCLE
       DO n=1,SIZE(Patches)
          IF(TRIM(patches(n) % pres) == TRIM(Ri% Type) &
               & .OR. TRIM(patches(n) % res) == TRIM(Ri% Type) ) THEN
             ok_Residue(i_L)=.TRUE.
             EXIT
          END IF
       END DO
    END DO
    CALL CPU_TIME(Time_Begin)
    DO i_l=1,SIZE(Residue)
       Ri=>Residue(i_L)
       IF(.NOT. ok_Residue(i_L)) CYCLE
       Res_Char(i_L) % FField=Ri % FField
       Res_Char(i_L) % Type=Ri % Type
       Res_Char(i_L) % Residue=Ri % Residue

       IF(Res_Char(i_L) % FField == 'CHARMM') THEN
          IF(Res_Char (i_L) % Type == 'mass') THEN
             CALL Mass_
          ELSE
             CALL Charmm_
          END IF
       ELSE IF(Residue(i_L) % FField == 'ORAC') THEN
          CALL Orac_
          CALL Charmm_
       ELSE
          STOP
       END IF
    END DO
    CALL Newresidues_
    CALL CPU_TIME(Time_End)
    WRITE(*,*) 'CPU Time = ',Time_End-Time_Begin,' in s '

  CONTAINS
    INCLUDE 'Tops_Orac.f90'
    INCLUDE 'Tops_Charmm.f90'
  END SUBROUTINE Tops_
  INCLUDE 'Newresidues_.f90'
END MODULE Tops
