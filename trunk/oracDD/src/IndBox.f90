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
MODULE IndBox
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
  PUBLIC IndBox_,IndBox_g_p,IndBox_g_t,IndBox_a_p,IndBox_a_t
  INTEGER :: natom_local
  INTEGER, ALLOCATABLE, SAVE :: indBox_a_p(:),indBox_a_t(:)
  INTEGER, ALLOCATABLE, SAVE :: indBox_g_p(:),indBox_g_t(:)
CONTAINS
  FUNCTION IndBox_() RESULT(out)
    LOGICAL :: out
    INTEGER :: n,m,count_a_p,count_a_t,count_G_p,count_G_t

    count_g_p=0
    count_g_t=0
    count_g_p=COUNT(Groupa(:) % knwn == 1)
    count_g_t=COUNT(Groupa(:) % knwn /= 0)
    
    IF(ALLOCATED(IndBox_g_p)) DEALLOCATE(IndBox_g_p)
    IF(ALLOCATED(IndBox_g_t)) DEALLOCATE(IndBox_g_t)

    ALLOCATE(IndBox_g_p(count_g_p))
    ALLOCATE(IndBox_g_t(count_g_t)) 
    
    count_g_p=0 ; count_g_t=0 

    DO n=1,SIZE(Groupa)
       m=Groupa(n) % knwn 
       IF(m /= 0) THEN
          count_g_t=count_g_t+1
          IndBox_g_t(count_g_t)=n
          IF(m == 1) THEN
             count_g_p=count_g_p+1
             IndBox_g_p(count_g_p)=n
          END IF
       END IF
    END DO
    out=count_g_p /= 0
    IF(.NOT. out) THEN
       errmsg_f='No Primary Groups found in the unit box'
       CALL Add_Errors(-1,errmsg_f)
    END IF
    
    count_a_p=0
    count_a_t=0
    count_a_p=COUNT(Groupa(Atoms(:) % Grp_No) % knwn == 1)
    count_a_t=COUNT(Groupa(Atoms(:) % Grp_No) % knwn /= 0)
    
    IF(ALLOCATED(IndBox_a_p)) DEALLOCATE(IndBox_a_p)
    IF(ALLOCATED(IndBox_a_t)) DEALLOCATE(IndBox_a_t)

    ALLOCATE(IndBox_a_p(count_a_p))
    ALLOCATE(IndBox_a_t(count_a_t)) 
    
    count_a_p=0 ; count_a_t=0 

    DO n=1,SIZE(Atoms)
       m=Groupa(Atoms(n) % Grp_No) % knwn 
       IF(m /= 0) THEN
          count_a_t=count_a_t+1
          IndBox_a_t(count_a_t)=n
          IF(m == 1) THEN
             count_a_p=count_a_p+1
             IndBox_a_p(count_a_p)=n
          END IF
       END IF
    END DO
    out=count_a_p /= 0
    IF(.NOT. out) THEN
       errmsg_f='No Primary Atoms found in the unit box'
       CALL Add_Errors(-1,errmsg_f)
    END IF
  END FUNCTION IndBox_
END MODULE IndBox
