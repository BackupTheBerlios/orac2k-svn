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
MODULE Run

!!$***********************************************************************
!!$   Time-stamp: <2007-01-09 10:56:45 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Nov 23 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program ORAC ----*

!!$---- DATA Only Modules -----------------------------------------------*

  USE Constants, ONLY: max_pars,max_data, max_char
  USE Parameters_globals

!!$---- Modules ---------------------------------------------------------*

  USE Tree
  USE Errors, ONLY: Add_Errors=>Add, error_other, error_unr, error_args&
       &, errmsg_f, Print_Errors
  USE Strings, ONLY: MY_Fxm
  USE Myparse 
  USE STRPAK, ONLY: SP_Getnum
!!$---- DATA Statements -------------------------------------------------*

  IMPLICIT none
  PRIVATE
  PUBLIC Run__Scan, Run__Input, Run_
  TYPE :: Run__Input
     REAL(8) :: Time=0.0_8
     REAL(8) :: Reject=0.0_8
     REAL(8) :: Print=0.0_8
     INTEGER :: Control=0
  END type Run__Input
  TYPE(Run__Input), SAVE :: Run_
CONTAINS

!!$---- EXTECUTABLE Statements ------------------------------------------*

  SUBROUTINE Run__Scan
    CHARACTER(len=max_pars) :: line,linea
    INTEGER :: n,nword
    TYPE(Branch), SAVE :: check

    CALL Tree__Check_Tree('&RUN',check)
    IF(.NOT. ASSOCIATED(check%children)) RETURN

    DO n=1,SIZE(check%children)
       line=TRIM(check%children(n))
       nword=MYParse_(line)
       
       linea=strngs(1)

       IF(MY_Fxm('TIME',linea)) THEN
          CALL Time
       ELSE IF(MY_Fxm('CONT',linea)) THEN
          CALL Control
       ELSE IF(MY_Fxm('REJE',linea)) THEN
          CALL Reject
       ELSE IF(MY_Fxm('PRIN',linea)) THEN
          CALL Print
       ELSE
          errmsg_f='Illegal commmands found:'//TRIM(linea)
          CALL Add_Errors(-1,errmsg_f)
       END IF
    END DO
  END SUBROUTINE Run__Scan
  SUBROUTINE Time
    INTEGER ::  nword,iflags,of

    nword=SIZE(strngs)
    SELECT CASE(nword)
    CASE(2)
       CALL SP_Getnum(strngs(2),Run_ % Time,iflags)
       IF(iflags /=0) THEN
          errmsg_f='Cannot convert to Real: '''&
               &//TRIM(strngs(1))//' '//TRIM(strngs(2))&
               &//''' '
          CALL Add_Errors(-1,errmsg_f)
       END IF
    CASE DEFAULT 
       errmsg_f=error_args % g (4)//' 1'
       CALL Add_Errors(-1,errmsg_f)
    END SELECT
    CALL Print_Errors()
  END SUBROUTINE Time
  SUBROUTINE Reject
    INTEGER ::  nword,iflags,of

    nword=SIZE(strngs)
    SELECT CASE(nword)
    CASE(2)
       CALL SP_Getnum(strngs(2),Run_ % Reject,iflags)
       IF(iflags /=0) THEN
          errmsg_f='Cannot convert to Real: '''&
               &//TRIM(strngs(1))//' '//TRIM(strngs(2))&
               &//''' '
          CALL Add_Errors(-1,errmsg_f)
       END IF
    CASE DEFAULT 
       errmsg_f=error_args % g (4)//' 1'
       CALL Add_Errors(-1,errmsg_f)
    END SELECT
    CALL Print_Errors()
  END SUBROUTINE Reject
  SUBROUTINE Print
    INTEGER ::  nword,iflags,of

    nword=SIZE(strngs)
    SELECT CASE(nword)
    CASE(2)
       CALL SP_Getnum(strngs(2),Run_ % Print,iflags)
       IF(iflags /=0) THEN
          errmsg_f='Cannot convert to Real: '''&
               &//TRIM(strngs(1))//' '//TRIM(strngs(2))&
               &//''' '
          CALL Add_Errors(-1,errmsg_f)
       END IF
    CASE DEFAULT 
       errmsg_f=error_args % g (4)//' 1'
       CALL Add_Errors(-1,errmsg_f)
    END SELECT
    CALL Print_Errors()
  END SUBROUTINE Print
  SUBROUTINE Control
    INTEGER ::  nword,iflags,of

    nword=SIZE(strngs)
    SELECT CASE(nword)
    CASE(2)
       CALL SP_Getnum(strngs(2),Run_ % Control,iflags)
       IF(iflags /=0) THEN
          errmsg_f='Cannot convert to integer: '''&
               &//TRIM(strngs(1))//' '//TRIM(strngs(2))&
               &//''' '
          CALL Add_Errors(-1,errmsg_f)
       END IF
    CASE DEFAULT 
       errmsg_f=error_args % g (4)//' 1'
       CALL Add_Errors(-1,errmsg_f)
    END SELECT
    CALL Print_Errors()
  END SUBROUTINE Control

END MODULE Run
