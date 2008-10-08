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
MODULE IndPatch

!!$***********************************************************************
!!$   Time-stamp: <2007-01-10 10:03:25 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Tue Jan  9 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*

  USE Errors, ONLY: Add_Errors=>Add, errmsg_f, Print_errors
  USE SecondarySeq
  USE Parameters_Globals
  IMPLICIT none
  PRIVATE
  PUBLIC Indpatch_, Indpatch__Type, Patch

  TYPE :: Indpatch__Type
     INTEGER :: one
     INTEGER :: two
     CHARACTER(len=Max_Char) :: l1
     CHARACTER(len=Max_Char) :: l2
  END TYPE Indpatch__Type
  TYPE(Indpatch__type), DIMENSION(:), ALLOCATABLE, SAVE, TARGET :: Ind_Patch
CONTAINS
  FUNCTION Indpatch_() RESULT(out)
    TYPE(Indpatch__type), DIMENSION(:), POINTER :: out
     CHARACTER(len=Max_Char) :: res_i
    INTEGER :: m,n,m1,m2,i

    out=>NULL()
    IF(.NOT. ALLOCATED(Patches)) RETURN
    m=0
    DO n=1,SIZE(patches)
       IF(patches (n) % Type == 'link') THEN
          m=m+1
       END IF
    END DO
    IF(m == 0) RETURN

    ALLOCATE(Ind_Patch(m))
    m=0
    DO n=1,SIZE(patches)
       IF(patches (n) % Type == 'link') THEN
          m=m+1
          Ind_Patch(m) % one=patches(n)%one
          Ind_Patch(m) % Two=patches(n)%Two
          Ind_Patch(m) % l1 ='Link '//TRIM(patches(n)%Res_l(1))
          Ind_Patch(m) % l2 ='Link '//TRIM(patches(n)%Res_l(2))
       END IF
    END DO
    DO i=1,SIZE(Ind_Patch)
       m1=Ind_Patch(i) % one
       m2=Ind_Patch(i) % two
       IF(m1 > SIZE(Secondary(1) % line)&
            & .OR. m2 > SIZE(Secondary(1) % line)) THEN
          errmsg_f='Link patch numbers larger than available solute residues '
          CALL Add_Errors(-1,errmsg_f)
       END IF
    END DO
    CALL Print_Errors()
    n=1
    DO m=1,SIZE(Secondary(n) % line)
       res_i=Secondary(n) % line(m)
       DO i=1,SIZE(Ind_Patch)
          IF(Ind_Patch(i) % one == m) THEN
             res_i=Ind_Patch(i) % l1 
             Secondary(n) % line(m)=TRIM(res_i)
          ELSE IF(Ind_Patch(i) % Two == m) THEN
             res_i=Ind_Patch(i) % l2
             Secondary(n) % line(m)=TRIM(res_i)
          END IF
       END DO
    END DO
    out=>Ind_Patch
  END FUNCTION Indpatch_
  FUNCTION Indpatch__Assign() RESULT(out)
    TYPE(Indpatch__type), DIMENSION(:), POINTER :: out
    out=>NULL()
    IF(.NOT. ALLOCATED(Ind_Patch)) RETURN
    out=>Ind_Patch
  END FUNCTION Indpatch__Assign

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE IndPatch
