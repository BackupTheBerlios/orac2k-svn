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
MODULE Strings

!!$***********************************************************************
!!$   Time-stamp: <2007-01-12 18:29:49 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Mon Nov 13 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

  USE Errors, ONLY: Add_Errors=> Add, errmsg_f, Print_Errors
  USE CONSTANTS, ONLY: max_pars, max_char, Comms
  USE TYPES
  USE STRPAK
  INTERFACE Myread
     MODULE PROCEDURE MyRead_I4
     MODULE PROCEDURE MyRead_R8
  END INTERFACE
  INTERFACE Myputnum
     MODULE PROCEDURE MyPutNum_i4
     MODULE PROCEDURE MyPutNum_r8
  END INTERFACE
CONTAINS
  FUNCTION MyPutnum_i4(n) RESULT(out)
    INTEGER  :: n
    CHARACTER(len=max_char) :: out
    CHARACTER(len=max_char) :: str1
    INTEGER :: iflag

    CALL STR_Fill(' ',str1)
    CALL PUTI4N(n,1,-1,'I',str1,iflag)
    IF(iflag /= 0) THEN
       errmsg_f='Cannot convert integer to character '
       CALL Add_errors(-1,errmsg_f)
       CALL Print_Errors()
    END IF
    out=TRIM(str1)
  END FUNCTION MyPutnum_i4
  FUNCTION MyPutnum_R8(r) RESULT(out)
    REAL(8) :: r
    CHARACTER(len=max_char) :: out
    CHARACTER(len=max_char) :: str1
    INTEGER :: iflag

    CALL PUTR8N(r,1,-1,'R5',str1,iflag)
    IF(iflag /= 0) THEN
       errmsg_f='Cannot convert real to character '
       CALL Add_errors(-1,errmsg_f)
       CALL Print_Errors()
    END IF
    out=TRIM(str1)
  END FUNCTION MyPutnum_R8
  SUBROUTINE Myread_I4(str1,out)
    CHARACTER(len=*) :: str1
    INTEGER :: out
    
    INTEGER :: iflag
    
    CALL STR_Geti4n(str1,out,iflag)
    IF(iflag /= 0) THEN
       errmsg_f='Cannot convert to integer: '''//TRIM(str1)//''' '
       CALL Add_errors(-1,errmsg_f)
       CALL Print_Errors()
    END IF
  END SUBROUTINE Myread_I4
  SUBROUTINE Myread_R8(str1, out)
    CHARACTER(len=*) :: str1
    REAL(8) :: out
    
    INTEGER :: iflag
    
    CALL STR_Getr8n(str1,out,iflag)
    IF(iflag /= 0) THEN
       errmsg_f='Cannot convert to double: '''//TRIM(str1)//''' '
       CALL Add_errors(-1,errmsg_f)
       CALL Print_Errors()
    END IF
  END SUBROUTINE Myread_R8

  FUNCTION My_FAM(str1,str2)
    IMPLICIT NONE 
    CHARACTER(len=*) :: str1,str2
    LOGICAL :: My_Fam
    INTEGER :: j1,j2,iflag

    My_FAM=.FALSE.
    CALL STR_FAM(str1,str2,j1,j2,iflag)
    IF(iflag == 0) MY_Fam=.TRUE.
  END FUNCTION My_FAM
  FUNCTION My_FXM(str1,str2)
    IMPLICIT NONE 
    CHARACTER(len=*) :: str1,str2
    LOGICAL :: My_Fxm
    INTEGER :: j1,j2,iflag

    My_FXM=.FALSE.
    CALL STR_FXM(str1,str2,j1,j2,iflag)
    IF(iflag == 0) MY_Fxm=.TRUE.
  END FUNCTION My_FXM
  FUNCTION ST_Concat(i,strngs)
    IMPLICIT NONE 
    CHARACTER(len=max_char) :: ST_Concat
    INTEGER :: i
    CHARACTER(len=max_pars), DIMENSION(:) :: strngs

    INTEGER :: n,m,pad
    CHARACTER(len=max_char) :: line,aux 

    CALL STR_FILL(' ',line)
    
    pad=1
    m=0
    DO n=i,SIZE(strngs)
       m=m+LEN_TRIM(strngs(n))+pad
    END DO
    IF(m > LEN(line)) THEN
       WRITE(aux,'(i3)') LEN(line)
       errmsg_f='Cannot contatenate strings. Final string length larger&
            & than hard bound '//TRIM(aux)
       CALL Add_Errors(-1,errmsg_f)
       CALL Print_Errors()
    END IF
    DO n=i,SIZE(strngs)
       line=TRIM(line)//' '//TRIM(strngs(n))
    END DO
    ST_Concat=line
  END FUNCTION ST_Concat
  FUNCTION ST_Clean_Line(line)
    IMPLICIT NONE 
    LOGICAL :: ST_Clean_Line
    CHARACTER(len=max_char) :: line
    INTEGER :: i

    ST_Clean_Line=.FALSE.
    DO i=1,SIZE(Comms)
       CALL STR_TRIM(Comms(i),line)
    END DO
    IF(LEN_TRIM(line) == 0) ST_Clean_Line=.TRUE.
  END FUNCTION ST_Clean_Line
END MODULE Strings
