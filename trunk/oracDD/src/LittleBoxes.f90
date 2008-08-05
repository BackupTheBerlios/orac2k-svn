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
MODULE LittleBoxes
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Tue Aug  5 2008 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*

  USE Atom
  USE Groups
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f, errmsg_w
  IMPLICIT none
  PRIVATE
  PUBLIC LittleBoxes_,LittleBoxes__Update,IndBox_p,IndBox_s,IndBox_t,Atoms_Knwn
  INTEGER :: natom_local
  INTEGER, ALLOCATABLE, SAVE :: indBox_p(:),indBox_s(:),indBox_t(:)
  INTEGER, ALLOCATABLE, SAVE :: Atoms_knwn(:)
CONTAINS
  FUNCTION LittleBoxes_() RESULT(out)
    LOGICAL :: out
    INTEGER :: Known0
    INTEGER :: n,m,AtSt,AtEn,count_p,count_s,count_t
    
    IF(.NOT. ALLOCATED(Atoms_Knwn)) THEN
       ALLOCATE(Atoms_Knwn(SIZE(Atoms)))
    END IF
    Atoms_Knwn(:)=Groupa(Atoms(:) % Grp_No) % knwn
    count_p=0
    count_s=0
    
    DO n=1,SIZE(Atoms_Knwn)
       IF(Atoms_Knwn (n) == 1) THEN
          count_p=count_p+1
       ELSE IF(Atoms_Knwn (n) == 2) THEN
          count_s=count_s+1
       END IF
    END DO
    count_t=count_p+count_s
    IF(ALLOCATED(IndBox_p)) DEALLOCATE(IndBox_p)
    IF(ALLOCATED(IndBox_s)) DEALLOCATE(IndBox_s)
    IF(ALLOCATED(IndBox_t)) DEALLOCATE(IndBox_t)

    ALLOCATE(IndBox_p(count_p)) ; ALLOCATE(IndBox_s(count_s)) 
    ALLOCATE(IndBox_t(count_t)) 
    
    count_p=0 ; count_s=0 

    DO n=1,SIZE(Atoms_Knwn)
       IF(Atoms_Knwn (n) == 1) THEN
          count_p=count_p+1
          IndBox_p(count_p)=n
       ELSE IF(Atoms_Knwn (n) == 2) THEN
          count_s=count_s+(AtEn-AtSt+1)
          IndBox_s(count_s)=n
       END IF
    END DO

    IndBox_t(1:SIZE(IndBox_p))=IndBox_p
    IndBox_t(SIZE(IndBox_p)+1:)=IndBox_s
    out=count_p /= 0
    IF(.NOT. out) THEN
       errmsg_f='No Atoms found in the unit box'
       CALL Add_Errors(-1,errmsg_f)
    END IF
  END FUNCTION LittleBoxes_
  FUNCTION LittleBoxes__Update() RESULT(out)
    LOGICAL :: out
    INTEGER :: Known0
    INTEGER :: n,m,AtSt,AtEn,count_p,count_s,count_t
    
    IF(.NOT. ALLOCATED(Atoms_Knwn)) THEN
       errmsg_f='Cannot do Update of LitteBoxes without initialisation'
       CALL Add_Errors(-1,errmsg_f)
       RETURN
    END IF
    count_p=0
    count_s=0    
    DO n=1,SIZE(Atoms_Knwn)
       IF(Atoms_Knwn (n) == 1) THEN
          count_p=count_p+1
       ELSE IF(Atoms_Knwn (n) == 2) THEN
          count_s=count_s+1
       END IF
    END DO
    count_t=count_p+count_s
    IF(ALLOCATED(IndBox_p)) DEALLOCATE(IndBox_p)
    IF(ALLOCATED(IndBox_s)) DEALLOCATE(IndBox_s)
    IF(ALLOCATED(IndBox_t)) DEALLOCATE(IndBox_t)

    ALLOCATE(IndBox_p(count_p)) ; ALLOCATE(IndBox_s(count_s)) 
    ALLOCATE(IndBox_t(count_t)) 
    
    count_p=0 ; count_s=0 

    DO n=1,SIZE(Atoms_Knwn)
       IF(Atoms_Knwn (n) == 1) THEN
          count_p=count_p+1
          IndBox_p(count_p)=n
       ELSE IF(Atoms_Knwn (n) == 2) THEN
          count_s=count_s+(AtEn-AtSt+1)
          IndBox_s(count_s)=n
       END IF
    END DO

    IndBox_t(1:SIZE(IndBox_p))=IndBox_p
    IndBox_t(SIZE(IndBox_p)+1:)=IndBox_s
    out=count_p /= 0
    IF(.NOT. out) THEN
       errmsg_f='No Atoms found in the unit box'
       CALL Add_Errors(-1,errmsg_f)
    END IF
  END FUNCTION LittleBoxes__Update
END MODULE LittleBoxes
