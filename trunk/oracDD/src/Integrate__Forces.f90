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
SUBROUTINE Forces_(n)
  INTEGER :: n
  INTEGER, PARAMETER :: Flag=1

  IF(n >= 3) THEN
     CALL InterForces_
  ELSE 
     CALL IntraForces_
  END IF
CONTAINS
  SUBROUTINE IntraForces_
    
    CALL PI__ResetSecondary
    SELECT CASE(n)
    CASE(_N0_)
!!$--- Shift atoms
       CALL PI__ShiftIntra(_N0_,Flag)
!!$--- Gets all n0 interactions beloging to the primary cell
       CALL Intra_n0_(Flag)
    CASE(_N1_)
!!$--- Shift atoms
       CALL PI__ShiftIntra(_N1_,Flag)
!!$--- Gets all n1 interactions beloging to the primary cell
       CALL Intra_n1_(Flag)
    END SELECT
    
  END SUBROUTINE IntraForces_
  
  SUBROUTINE InterForces_
    LOGICAL :: pme
    TYPE(Force), POINTER :: fp(:)
    INTEGER :: m
    INTEGER :: callsa=0

    CALL PI__ResetSecondary

    WRITE(*,*) 'Lucillo ',PI_Node_CArt, COUNT(Groupa(:) % knwn == 2)

    pme=(n-2 == Integrator_ % Ewald_Shell) .AND. Ewald__Param % Switch
    IF(pme) THEN
       IF(Ewald__Param % nx /= 0 .AND. Ewald__Param % ny  /= 0 .AND.&
            & Ewald__Param % nz /= 0) THEN
          CALL PI__Shift(n,Flag,_PME_)
          WRITE(*,*) 'Illak'
       END IF
    ELSE
       CALL PI__Shift(n,Flag)
    END IF

    IF(.NOT. PI_Atom_Update_()) CALL Print_Errors()

    IF(n == 3) THEN
       WRITE(*,*) '2029 F',Groupa(Atoms(5756) % Grp_No) % knwn &
            & ,Groupa(2029) % xa,Groupa(2029) % ya
       WRITE(*,*) '2029 F',Atoms(5756) % Grp_No
    END IF


    CALL DIR_Forces(n)
    IF(pme) CALL PME_(n)
     
!!$
!!$--- Fold forces contributions to atoms inside the cell
!!$
    
    fp=>FORCES_Pick(n)
    CALL PI__Fold_F(fp,n,Flag)

!!$--- Reset Secondary region atoms from PME
    
    IF(pme) CALL Pi__ResetSecondaryP
  END SUBROUTINE InterForces_
END SUBROUTINE Forces_
