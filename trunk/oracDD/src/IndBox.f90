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

  USE PI_
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f, errmsg_w
  IMPLICIT none
  PRIVATE
  PUBLIC IndBox_,IndBox_g_p,IndBox_g_t,IndBox_a_p,IndBox_a_t&
       &,BoxInd_a_p
  INTEGER :: natom_local
  INTEGER, ALLOCATABLE, SAVE :: indBox_a_p(:),indBox_a_t(:),BoxInd_a_p(:)
  INTEGER, ALLOCATABLE, SAVE :: indBox_g_p(:),indBox_g_t(:)
CONTAINS
  FUNCTION IndBox_(g_knwn,g_AtSt,g_AtEn) RESULT(out)
    INTEGER :: g_knwn(:),g_AtSt(:),g_AtEn(:)
    LOGICAL :: out
    INTEGER :: n,m,count_a_p,count_a_t,count_G_p,count_G_t,q,AtSt&
         &,AtEn,natom

    natom=SUM(g_AtEn-g_AtSt)+SIZE(g_knwn)
    count_g_p=0
    count_g_t=0
    count_g_p=COUNT(g_knwn(:) == 1)
    count_g_t=COUNT(g_knwn(:) /= 0)

    IF(ALLOCATED(IndBox_g_p)) DEALLOCATE(IndBox_g_p)
    IF(ALLOCATED(IndBox_g_t)) DEALLOCATE(IndBox_g_t)

    ALLOCATE(IndBox_g_p(count_g_p))
    ALLOCATE(IndBox_g_t(count_g_t)) 
    
    count_g_p=0 ; count_g_t=0 

    DO n=1,SIZE(G_Knwn)
       m=G_knwn(n)
       IF(m /= 0) THEN
          count_g_t=count_g_t+1
          IndBox_g_t(count_g_t)=n
          IF(m == 1) THEN
             count_g_p=count_g_p+1
             IndBox_g_p(count_g_p)=count_g_t
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
    DO n=1,SIZE(G_Knwn)
       AtSt=G_AtSt(n)
       AtEn=G_AtEn(n)
       m=G_knwn(n)
       IF(m == 1) count_a_p=count_a_p+(AtEn-AtSt+1)
       IF(m /= 0) count_a_t=count_a_t+(AtEn-AtSt+1)
    END DO
    IF(ALLOCATED(IndBox_a_p)) DEALLOCATE(IndBox_a_p)
    IF(ALLOCATED(IndBox_a_t)) DEALLOCATE(IndBox_a_t)
    IF(ALLOCATED(BoxInd_a_p)) DEALLOCATE(BoxInd_a_p)

    ALLOCATE(IndBox_a_p(count_a_p))
    ALLOCATE(IndBox_a_t(count_a_t)) 
    ALLOCATE(BoxInd_a_p(natom))
    
    count_a_p=0 ; count_a_t=0 

    BoxInd_a_p=-1
    DO n=1,SIZE(G_Knwn)
       AtSt=G_AtSt(n)
       AtEn=G_AtEn(n)
       m=G_knwn(n)
       DO q=AtSt,AtEn
          IF(m /= 0) THEN
             count_a_t=count_a_t+1
             IndBox_a_t(count_a_t)=q
             IF(m == 1) THEN
                count_a_p=count_a_p+1
                IndBox_a_p(count_a_p)=count_a_t
                BoxInd_a_p(q)=count_a_p
             END IF
          END IF
       END DO
    END DO

    out=count_a_p /= 0
    IF(.NOT. out) THEN
       errmsg_f='No Primary Atoms found in the unit box'
       CALL Add_Errors(-1,errmsg_f)
    END IF
  END FUNCTION IndBox_
END MODULE IndBox
