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
MODULE Integrator

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

  USE Node
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
  PUBLIC Integrator__Scan
  TYPE :: Integrator__Input
     REAL(8) :: t=12.0_8
     INTEGER :: Mult_Intra(2)=(/2,2/)
     INTEGER :: Mult_Inter(3)=(/2,3,1/)
     INTEGER :: Ewald_shell
  END type Integrator__Input
  TYPE(Integrator__Input), SAVE :: Integrator_
CONTAINS

!!$---- EXTECUTABLE Statements ------------------------------------------*

  SUBROUTINE Integrator__Scan
    CHARACTER(len=max_pars) :: line,linea
    INTEGER :: n,nword
    TYPE(Branch), SAVE :: check

    CALL Tree__Check_Tree('&INTEGRATOR',check)
    IF(.NOT. ASSOCIATED(check%children)) RETURN

    DO n=1,SIZE(check%children)
       line=TRIM(check%children(n))
       nword=MYParse_(line)
       
       linea=strngs(1)

       IF(MY_Fxm('EWALD',linea)) THEN
          CALL Ewald
       ELSE IF(MY_Fxm('TIME',linea)) THEN
          CALL Timestep
       ELSE IF(MY_Fxm('INTRA',linea)) THEN
          CALL Intra
       ELSE IF(MY_Fxm('INTER',linea)) THEN
          CALL Inter
       ELSE
          errmsg_f='Illegal commmands found:'//TRIM(linea)
          CALL Add_Errors(-1,errmsg_f)
       END IF
    END DO
  CONTAINS
    SUBROUTINE Validate
      LOGICAL :: ex
    END SUBROUTINE Validate
  END SUBROUTINE Integrator__Scan
  SUBROUTINE Ewald
    INTEGER ::  nword,iflags,of

    nword=SIZE(strngs)
    SELECT CASE(nword)
    CASE(2)
       CALL SP_Getnum(strngs(2),Integrator_ % Ewald_Shell,iflags)
       IF(iflags /=0) THEN
          errmsg_f='Cannot convert to integer: '''&
               &//TRIM(strngs(of+1))//' '//TRIM(strngs(of+2))&
               &//''' '
          CALL Add_Errors(-1,errmsg_f)
       END IF
    CASE DEFAULT 
       errmsg_f=error_args % g (4)//' 1'
       CALL Add_Errors(-1,errmsg_f)
    END SELECT
    CALL Print_Errors()
  END SUBROUTINE Ewald

  SUBROUTINE Intra
    INTEGER ::  nword,iflags

    nword=SIZE(strngs)
    SELECT CASE(nword)
    CASE(3)
       CALL SP_Getnum(strngs(2),Integrator_ % Intra(1),iflags)
       IF(iflags /=0) THEN
          errmsg_f='Cannot convert to integer: '''&
               &//TRIM(strngs(3))//''' '
          CALL Add_Errors(-1,errmsg_f)
       END IF
       CALL SP_Getnum(strngs(3),Integrator_ % Intra(2),iflags)
       IF(iflags /=0) THEN
          errmsg_f='Cannot convert to integer: '''&
               &//TRIM(strngs(3))//''' '
          CALL Add_Errors(-1,errmsg_f)
       END IF
    CASE DEFAULT 
       errmsg_f=error_args % g (4)//' 2 '
       CALL Add_Errors(-1,errmsg_f)
    END SELECT
    CALL Print_Errors()
  END SUBROUTINE Intra
  SUBROUTINE Inter
    INTEGER ::  nword,iflags

    nword=SIZE(strngs)
    SELECT CASE(nword)
    CASE(3)
       CALL SP_Getnum(strngs(2),Integrator_ % Inter(1),iflags)
       IF(iflags /=0) THEN
          errmsg_f='Cannot convert to integer: '''&
               &//TRIM(strngs(3))//''' '
          CALL Add_Errors(-1,errmsg_f)
       END IF
       CALL SP_Getnum(strngs(3),Integrator_ % Inter(2),iflags)
       IF(iflags /=0) THEN
          errmsg_f='Cannot convert to integer: '''&
               &//TRIM(strngs(3))//''' '
          CALL Add_Errors(-1,errmsg_f)
       END IF
    CASE DEFAULT 
       errmsg_f=error_args % g (4)//' 2 '
       CALL Add_Errors(-1,errmsg_f)
    END SELECT
    CALL Print_Errors()
  END SUBROUTINE Inter
  SUBROUTINE Timestep
    INTEGER ::  nword,iflags

    nword=SIZE(strngs)
    SELECT CASE(nword)
    CASE(2)
       CALL SP_Getnum(strngs(2),Integrator_ % t,iflags)
       IF(iflags /=0) THEN
          errmsg_f='Cannot convert to Real: '''&
               &//TRIM(strngs(3))//''' '
          CALL Add_Errors(-1,errmsg_f)
       END IF
    CASE DEFAULT 
       errmsg_f=error_args % g (4)//' 1 '
       CALL Add_Errors(-1,errmsg_f)
    END SELECT
    CALL Print_Errors()
  END SUBROUTINE Timestep
END MODULE Integrator
