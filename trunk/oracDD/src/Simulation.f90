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
MODULE Simulation

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

  USE Forces, ONLY: FOR_Init=>Init
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
  PUBLIC Simulation__Scan, T, Pressure, Temp_,Press_

  TYPE :: Temp_
     REAL(8) :: t0=300.0_8,dtemp=10.0_8
     LOGICAL :: Gauss_Init=.FALSE.
     INTEGER :: steps=0
  END type Temp_
  TYPE :: Press_
     REAL(8) :: p0=0.1_8,dpress=0.0_8
     INTEGER :: steps=0
  END type Press_
  TYPE(Temp_), SAVE :: T
  TYPE(Press_), SAVE :: Pressure
CONTAINS

!!$---- EXTECUTABLE Statements ------------------------------------------*

  SUBROUTINE Simulation__Scan
    CHARACTER(len=max_pars) :: line,linea
    INTEGER :: n,nword
    TYPE(Branch), SAVE :: check

    CALL Tree__Check_Tree('&SIMULATION',check)
    IF(.NOT. ASSOCIATED(check%children)) RETURN

    DO n=1,SIZE(check%children)
       line=TRIM(check%children(n))
       nword=MYParse_(line)
       
       linea=strngs(1)

       IF(MY_Fxm('TEMP',linea)) THEN
          CALL Temperature_(TRIM(Line))
       ELSE IF(MY_Fxm('PRESS',linea)) THEN
          CONTINUE
       ELSE
          errmsg_f='Illegal commmands found:'//TRIM(linea)
          CALL Add_Errors(-1,errmsg_f)
       END IF
    END DO
  END SUBROUTINE Simulation__Scan
  SUBROUTINE Temperature_(name)
    CHARACTER(len=*) :: name
    TYPE(Branch), SAVE :: check
    CHARACTER(len=max_pars) :: line,linea
    CHARACTER(len=max_char) :: lab0
    INTEGER :: n,nword,iflags,m
    INTEGER, SAVE :: count0

    CALL Tree__Check_Tree(name,check)
    IF(.NOT. ASSOCIATED(check%children)) RETURN

    count0=0
    DO n=1,SIZE(check%children)
       line=TRIM(check%children(n))
       nword=MYParse_(line)
       linea=strngs(1)

       IF(MY_Fxm('set',linea)) THEN
          SELECT CASE(nword)
          CASE(1)
             errmsg_f=error_args % g (2) //' 2'
             CALL Add_Errors(-1,errmsg_f)
             RETURN
          CASE(2)
             CALL SP_Getnum(strngs(2),T % t0,iflags)
          CASE DEFAULT
             WRITE(lab0,'(i2)') nword
             errmsg_f=error_args % g (3) //' 2 whereas it was '&
                  &//TRIM(lab0)//' : '//TRIM(line)
             CALL Add_Errors(-1,errmsg_f)
             RETURN
          END SELECT
       ELSE IF(MY_Fxm('rescale',linea)) THEN
          SELECT CASE(nword)
          CASE(1)
             errmsg_f=error_args % g (2) //' 3'
             CALL Add_Errors(-1,errmsg_f)
             RETURN
          CASE(2)
             IF(MY_Fxm('gauss',strngs(2))) THEN
                T % Gauss_Init = .TRUE.
             END IF
          CASE(3)
             IF(MY_Fxm('dt',strngs(2))) THEN
                CALL SP_Getnum(strngs(3),T % dtemp,iflags)
             ELSE IF(MY_Fxm('every',strngs(2))) THEN
                CALL SP_Getnum(strngs(3),T % steps,iflags)
             END IF
          CASE DEFAULT
             WRITE(lab0,'(i2)') nword
             errmsg_f=error_args % g (3) //' 3 whereas it was '&
                  &//TRIM(lab0)//' : '//TRIM(line)
             CALL Add_Errors(-1,errmsg_f)
             RETURN
          END SELECT
       ELSE
          errmsg_f='Illegal commmands found:'//TRIM(linea)
          CALL Add_Errors(-1,errmsg_f)
       END IF
    END DO
    
  END SUBROUTINE Temperature_
END MODULE Simulation
