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
MODULE Inout
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Jan 25 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*

  USE Constants
  USE Tree
  USE Errors, ONLY: Add_Errors=>Add, error_other, error_unr, error_args&
       &, errmsg_f, errmsg_w
  USE Strings, ONLY: MY_Fxm
  USE Myparse 
  USE STRPAK
  IMPLICIT none
  PRIVATE
  PUBLIC Inout__Scan, Inout__, Inout__PDB
  TYPE :: Inout__
     CHARACTER(len=max_char) :: tag=' '
     INTEGER :: unit=0
     REAL(8) :: freq=0.0D0
     CHARACTER(len=max_char) :: file=' '
  END TYPE Inout__
  TYPE(Inout__), SAVE :: Inout__PDB
CONTAINS
  SUBROUTINE Inout__Scan
    CHARACTER(len=max_pars) :: line,linea
    INTEGER :: n,nword,iflags,io
    TYPE(Branch), SAVE :: check

    CALL Tree__Check_Tree('&INOUT',check)
    IF(.NOT. ASSOCIATED(check%children)) RETURN
    DO n=1,SIZE(check%children)
       line=TRIM(check%children(n))
       nword=MYParse_(line)
       linea=strngs(1)
       
       IF(MY_Fxm('PDB',linea)) THEN
          SELECT CASE(nword)
          CASE(3)
             CALL SP_Getnum(strngs(2),Inout__PDB % freq,iflags)
             IF(iflags /= 0) THEN
                errmsg_f='Internal reading error: Module Cell'
                CALL Add_Errors(-1,errmsg_f)
                RETURN
             END IF
             Inout__PDB % File=TRIM(ADJUSTL(strngs(3)))
             CALL CHANNEL(io)
             Inout__PDB % unit = io
             OPEN(unit=io,file=Inout__PDB  % File,form='FORMATTED'&
                  &,status='UNKNOWN', action='WRITE')
          CASE(4)
             IF(My_Fxm('WSC',strngs(2))) THEN
                CALL SP_Getnum(strngs(2),Inout__PDB  % freq,iflags)
                IF(iflags /= 0) THEN
                   errmsg_f='Internal reading error: Module Cell'
                   CALL Add_Errors(-1,errmsg_f)
                   RETURN
                END IF
                Inout__PDB % File=TRIM(ADJUSTL(strngs(3)))            
                Inout__PDB  % tag = 'WSC'
                CALL CHANNEL(io)
                Inout__PDB  % unit = io
                OPEN(unit=io,file=Inout__PDB % File,form='FORMATTED'&
                     &,status='UNKNOWN', action='WRITE')
             ELSE
                errmsg_f=error_unr % g (3)//' : '//TRIM(strngs(2))&
                     &//' is not a defined keword '
                CALL Add_Errors(-1,errmsg_f)                
             END IF
          CASE DEFAULT 
             errmsg_f=error_args % g (4)//' 1, 3 or 6'
             CALL Add_Errors(-1,errmsg_f)
          END SELECT
       ELSE IF(MY_Fxm('DUMP',linea)) THEN
          CONTINUE
       ELSE
          errmsg_f='Illegal commmands found:'//TRIM(linea)
          CALL Add_Errors(-1,errmsg_f)
       END IF
    END DO
  END SUBROUTINE Inout__Scan
END MODULE Inout
