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
MODULE Potential

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
  USE STRPAK

!!$---- DATA Statements -------------------------------------------------*

  IMPLICIT none
  PRIVATE
  PUBLIC Potential__Scan, Ewald__Input, Ewald__Param
  TYPE :: Ewald__Input
     LOGICAL :: Do_not_Change=.FALSE.
     REAL(8) :: alpha
     REAL(8) :: Density
     INTEGER :: nx=0,ny=0,nz=0, order=5
  END TYPE Ewald__Input
  TYPE(Ewald__Input), SAVE :: Ewald__Param
CONTAINS

!!$---- EXTECUTABLE Statements ------------------------------------------*

  SUBROUTINE Potential__Scan
    CHARACTER(len=max_pars) :: line,linea
    INTEGER :: n,nword
    TYPE(Branch), SAVE :: check
    CALL Tree__Check_Tree('&POTENTIAL',check)
    IF(.NOT. ASSOCIATED(check%children)) RETURN

    DO n=1,SIZE(check%children)
       line=TRIM(check%children(n))
       nword=MYParse_(line)
       
       linea=strngs(1)

       IF(MY_Fxm('EWALD',linea)) THEN
          CALL Ewald
       ELSE
          errmsg_f='Illegal commmands found:'//TRIM(linea)
          CALL Add_Errors(-1,errmsg_f)
       END IF
    END DO
    CALL Validate
    CONTAINS
      SUBROUTINE Validate
        LOGICAL :: ex
      END SUBROUTINE Validate
    END SUBROUTINE Potential__Scan
  SUBROUTINE Ewald
    INTEGER ::  nword,iflags,of

    of=1
    nword=SIZE(strngs)
    IF(MY_Fxm('FIX',strngs(2))) THEN
       Ewald__Param % do_not_change=.TRUE.
       of=2
    END IF
    SELECT CASE(nword-of)
    CASE(2)
       CALL SP_Getnum(strngs(of+1),Ewald__Param % alpha,iflags)
       CALL SP_Getnum(strngs(of+2),Ewald__Param % Density,iflags)
       IF(iflags /=0) THEN
          errmsg_f='Cannot convert to integer: '''&
               &//TRIM(strngs(of+1))//' '//TRIM(strngs(of+2))&
               &//''' '
          CALL Add_Errors(-1,errmsg_f)
       END IF
    CASE(3)
       CALL SP_Getnum(strngs(of+1),Ewald__Param % alpha,iflags)
       CALL SP_Getnum(strngs(of+2),Ewald__Param % Density,iflags)
       CALL SP_Getnum(strngs(of+3),Ewald__Param % Order,iflags)
       IF(iflags /=0) THEN
          errmsg_f='Cannot convert to integer: '''&
               &//TRIM(strngs(of+1))//' '//TRIM(strngs(of+2))&
               &//' '//TRIM(strngs(of+3))//''' '
          CALL Add_Errors(-1,errmsg_f)
       END IF
    CASE(4)
       CALL SP_Getnum(strngs(of+1),Ewald__Param % alpha,iflags)
       CALL SP_Getnum(strngs(of+2),Ewald__Param % nx,iflags)
       CALL SP_Getnum(strngs(of+3),Ewald__Param % ny,iflags)
       CALL SP_Getnum(strngs(of+4),Ewald__Param % nz,iflags)
       IF(iflags /=0) THEN
          errmsg_f='Cannot convert to integer: '''&
               &//TRIM(strngs(of+1))//' '//TRIM(strngs(of+2))&
               &//' '//TRIM(strngs(of+3))//' '//TRIM(strngs(of+4))//''' '
          CALL Add_Errors(-1,errmsg_f)
       END IF
    CASE(5)
       CALL SP_Getnum(strngs(of+1),Ewald__Param % alpha,iflags)
       CALL SP_Getnum(strngs(of+2),Ewald__Param % nx,iflags)
       CALL SP_Getnum(strngs(of+3),Ewald__Param % ny,iflags)
       CALL SP_Getnum(strngs(of+4),Ewald__Param % nz,iflags)
       CALL SP_Getnum(strngs(of+5),Ewald__Param % Order,iflags)
       IF(iflags /=0) THEN
          errmsg_f='Cannot convert to integer: '''&
               &//TRIM(strngs(of+1))//' '//TRIM(strngs(of+2))&
               &//' '//TRIM(strngs(of+3))//' '//TRIM(strngs(of+4))&
               &//' '//TRIM(strngs(of+5))//''' '
          CALL Add_Errors(-1,errmsg_f)
       END IF
    CASE DEFAULT 
       errmsg_f=error_args % g (4)//' 2 or 4'
       CALL Add_Errors(-1,errmsg_f)
    END SELECT
    CALL Print_Errors()    
  END SUBROUTINE Ewald
END MODULE Potential
