!!$/---------------------------------------------------------------------\
!!$   Copyright  � 2006-2007 Massimo Marchi <Massimo.Marchi at cea.fr>   |
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
MODULE Inputs
  
!!$***********************************************************************
!!$   Time-stamp: <2007-01-04 15:44:58 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Tue Nov 21 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This MODULE is part of the program ORAC ----*

  USE Errors, ONLY: abort_now
  USE STRPAK
  USE Types
  USE Strings, ONLY: MY_Fam, ST_Clean_Line
  USE Print_Defs
  IMPLICIT none
  PRIVATE
  PUBLIC Inputs_,Inputs__String,Inputs__Filename

  CHARACTER(len=max_data), SAVE :: Inputs__String

  CHARACTER(len=120) :: errmsg
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: ind
  INTEGER, EXTERNAL :: iargc
  LOGICAL, SAVE :: prints=.FALSE.
  INTEGER, SAVE :: kinput
  CHARACTER(len=max_char), SAVE :: file
  CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE, SAVE :: data
CONTAINS
  SUBROUTINE Inputs_
    CALL Options
    CALL Store
    CALL Modify
    CALL Single_String
  END SUBROUTINE Inputs_
  FUNCTION Inputs__filename() RESULT(out)
    CHARACTER(len=max_char) :: out
    out=file
  END FUNCTION Inputs__filename
  SUBROUTINE Store
    CHARACTER(len=max_char) :: line,ui
  
    LOGICAL :: end_of_file=.FALSE.
    INTEGER :: io,c,iopt,iflag
    
    CALL CHANNEL(io)
    kinput=io
    OPEN(unit=kinput, file=file, form='FORMATTED')
    c=0
    DO WHILE(.NOT. end_of_file)
       READ(kinput,'(a)',IOSTAT=iopt) line
       IF(iopt /= 0) EXIT
       IF(ST_Clean_Line(line)) CYCLE
       iflag=0
       DO WHILE(iflag == 0) 
          CALL STR_RAS(CHAR(09),'   ',line,iflag)
       END DO
       IF(prints) WRITE(kprint,*) TRIM(line)
       c=c+1
    END DO
    ALLOCATE(data(c),ind(c))

    REWIND(kinput)
    c=0
    DO WHILE(.NOT. end_of_file)
       READ(kinput,'(a)',IOSTAT=iopt) line
       IF(ST_Clean_Line(line)) CYCLE
       IF(iopt /= 0) EXIT
       c=c+1
       iflag=0
       DO WHILE(iflag == 0) 
          CALL STR_RAS(CHAR(09),'   ',line,iflag)
       END DO
       data(c)=ADJUSTL(TRIM(line))
       ind(c)=LEN_TRIM(line)-LEN_TRIM(ADJUSTL(line))
    END DO
  END SUBROUTINE Store
  SUBROUTINE Modify
    INTEGER :: n,n_d,n_trim,m,o
    LOGICAL :: ok
    CHARACTER(len=max_char) :: dummy
    
    n_d=SIZE(data)

    DO n=1,n_d
       IF(data(n)(1:1) == '&') THEN
          IF(data(n)(1:4) == '&END') THEN
             data(n)(5:5)='}'
          ELSE
             CALL STR_SHFT(1,data(n))
             data(n)(1:1)='{'
          END IF
       END IF
    END DO

    n=0
    DO WHILE(n < n_d)
       n=n+1
       IF(MY_Fam('&',data(n))) CYCLE       
       IF(data(n)(1:3) == 'END') THEN
          data(n)(4:4)='}'
          DO o=n-1,1,-1
             IF(MY_Fam('&',data(o))) CYCLE
             dummy=data(o)
             CALL TRANUC(dummy)
             IF(dummy == data(o) .AND. (ind(n) == ind(o))) THEN
                m=o
                EXIT
             END IF
          END DO
          CALL STR_SHFT(1,data(m))
          data(m)(1:1)='{'
          DO o=m+1,n-1
             CALL STR_SHFT(1,data(o))
             data(o)(1:1)='{'
             n_trim=LEN_TRIM(data(o))
             data(o)(n_trim+1:n_trim+1)='}'
          END DO
       END IF

    END DO

    n=0
    DO WHILE(n < n_d)
       n=n+1
       IF(MY_Fam('&',data(n)) .OR. MY_Fam('{',data(n)) &
            &.OR. MY_Fam('}',data(n))) CYCLE
       CALL STR_SHFT(1,data(n))
       data(n)(1:1)='{'
       n_trim=LEN_TRIM(data(n))
       data(n)(n_trim+1:n_trim+1)='}'
    END DO
  END SUBROUTINE Modify
  SUBROUTINE Single_String
    INTEGER :: n,m
    Inputs__String=' '
    m=0
    DO n=1,SIZE(data)
       Inputs__String=TRIM(Inputs__String)//' '//TRIM(data(n))
       m=m+LEN_TRIM(data(n))
    END DO
    
    DEALLOCATE(data,ind)

  END SUBROUTINE Single_String
  SUBROUTINE Options
    IMPLICIT NONE 
    INTEGER :: n,ntot
    CHARACTER(len=max_char) :: line, linea
    LOGICAL :: ok

    n=0

#ifdef HAVE_F2003_EXT
    ntot=COMMAND_ARGUMENT_COUNT()
#else
    ntot=IARGC()
#endif

    IF(ntot == 0) THEN
       CALL GETARG(0,line)
       errmsg='Usage: '//TRIM(line)//' [-V] input > output '
       CALL abort_now(errmsg)
    END IF
    ok=.FALSE.
    DO WHILE(n < ntot)
       n=n+1
       CALL GETARG(n,line)
       linea=TRIM(line)
       IF(linea(1:1) == '-') THEN
          line=linea(2:)
          IF(line == 'V') THEN
             prints=.TRUE.
          ELSE
             errmsg='Unknown input option. Found '//TRIM(linea)
             CALL abort_now(errmsg)
          END IF
       ELSE
          file=linea
          ok=.TRUE.
          EXIT
       END IF
    END DO
    IF(.NOT. ok) THEN
       CALL GETARG(0,line)
       errmsg='Usage: '//TRIM(line)//' [-V] input > output '
       CALL abort_now(errmsg)
    END IF
  END SUBROUTINE Options

END MODULE Inputs
