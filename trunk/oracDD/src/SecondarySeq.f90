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
MODULE SecondarySeq

!!$***********************************************************************
!!$   Time-stamp: <2007-01-12 12:42:45 marchi>                           *
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

!!$---- This subroutine is part of the program  ----*


!!$======================== DECLARATIONS ================================*

  USE Errors, ONLY: Add_Errors=>Add, errmsg_f, Print_Errors
  USE Parameters_Globals
  USE MyParse
  USE STRPAK
  USE Node
  IMPLICIT none
  PRIVATE
  PUBLIC SecondarySeq_, Secondary, SecondarySeq__type, SecondarySeq__AddSlv&
       &,SecondarySeq__Read, SecondarySeq__Write
  TYPE :: SecondarySeq__Type
     CHARACTER(len=max_char) :: Type
     CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE :: line
  END TYPE SecondarySeq__Type
  TYPE(SecondarySeq__Type), DIMENSION(2), SAVE :: Secondary
CONTAINS
  SUBROUTINE SecondarySeq_(ip,Type,Seq)
    CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE  :: Seq
    CHARACTER(len=max_char) :: Type
    INTEGER :: ip
    INTEGER :: c,m,p,iflag,rept,n,nword
    CHARACTER(len=max_char) :: Line

    IF(.NOT. ALLOCATED(Seq)) RETURN
    Secondary(ip) % type=TRIM(Type)

    IF(.NOT. Node_()) STOP
    DO n=1,SIZE(Seq)
       line=Seq(n)
       nword=MYParse_(line)
       p=0
       DO WHILE(p < nword)
          p=p+1
          IF(strngs(p) /= 'x') THEN
             line=TRIM(strngs(p))
             CALL Node__Push(line)
          ELSE
             IF(p == SIZE(strngs)) THEN
                errmsg_f='Sequence ends. A integer is expected after a ''x'' '
                CALL Add_Errors(-1,errmsg_f)
                RETURN
             END IF
             CALL SP_Getnum(strngs(p+1), rept,iflag)
             IF(iflag /= 0) THEN
                errmsg_f='Reading error. A integer is expected after a ''x'' '
                CALL Add_Errors(-1,errmsg_f)
                RETURN
             END IF
             DO m=1,rept-1
                line=strngs(p-1)
                CALL Node__Push(line)
             END DO
             p=p+1
          END IF
       END DO
    END DO
    c=Node__SIZE()
    ALLOCATE(Secondary (ip) % line(c))
    c=0
    DO WHILE(Node__Pop(line))
       c=c+1
       Secondary(ip) % line (c) = TRIM(line)
    END DO
  END SUBROUTINE SecondarySeq_
  SUBROUTINE SecondarySeq__AddSlv(nunits)
    INTEGER :: nunits
    INTEGER :: n,m,begins,ends

    TYPE(SecondarySeq__Type)  :: Temp
    IF(ALLOCATED(Secondary(2) % line)) THEN
       m=SIZE(Secondary(2) % line)
       ALLOCATE(Temp % line(m*nunits))
       DO n=1,nunits
          begins=(n-1)*m+1
          ends=begins+m-1
          Temp % line (begins:ends)=Secondary(2) % line
       END DO
       DEALLOCATE(Secondary(2) % line)
       ALLOCATE(Secondary(2) % line(m*nunits))
       Secondary(2) % line=Temp % line
    END IF
  END SUBROUTINE Secondaryseq__AddSlv
  SUBROUTINE SecondarySeq__Write(kbinary)
    INTEGER :: kbinary
    INTEGER :: n,s
    DO n=1,2
       s=0
       IF(ALLOCATED(Secondary(n) % line)) s=SIZE(Secondary(n) % line)
       WRITE(kbinary) s
       IF(s /= 0) THEN
          WRITE(kbinary) Secondary(n) % line
       END IF
    END DO
  END SUBROUTINE SecondarySeq__Write
  FUNCTION SecondarySeq__Read(kbinary) RESULT(out)
    INTEGER :: kbinary
    LOGICAL :: out
    INTEGER :: n,s

    DO n=1,2
       READ(kbinary, ERR=100, END=200) s
       IF(s /= 0) THEN
          ALLOCATE(Secondary(n) % line (s))
          READ(kbinary, ERR=100, END=200) Secondary(n) % line
       END IF
    END DO
    out=.TRUE.
    RETURN
100 errmsg_f='Error while reading Lennard-Jones Parameters'
    CALL Add_Errors(-1,errmsg_f)
    out=.FALSE.
    RETURN
200 errmsg_f='End of file found while reading Lennard-Jones Parameters'
    CALL Add_Errors(-1,errmsg_f)
    out=.FALSE.
    RETURN    
  END FUNCTION SecondarySeq__Read
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE SecondarySeq
