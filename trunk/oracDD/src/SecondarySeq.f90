!!$/---------------------------------------------------------------------\
!!$                                                                      |
!!$  Copyright (C) 2006-2007 Massimo Marchi <Massimo.Marchi@cea.fr>      |
!!$                                                                      |
!!$      This program is free software;  you  can  redistribute  it      |
!!$      and/or modify it under the terms of the GNU General Public      |
!!$      License version 2 as published  by  the  Free  Software         |
!!$      Foundation;                                                     |
!!$                                                                      |
!!$      This program is distributed in the hope that  it  will  be      |
!!$      useful, but WITHOUT ANY WARRANTY; without even the implied      |
!!$      warranty of MERCHANTABILITY or FITNESS  FOR  A  PARTICULAR      |
!!$      PURPOSE.   See  the  GNU  General  Public License for more      |
!!$      details.                                                        |
!!$                                                                      |
!!$      You should have received a copy of the GNU General  Public      |
!!$      License along with this program; if not, write to the Free      |
!!$      Software Foundation, Inc., 59  Temple  Place,  Suite  330,      |
!!$      Boston, MA  02111-1307  USA                                     |
!!$                                                                      |
!!$\---------------------------------------------------------------------/
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
  PUBLIC SecondarySeq_, Secondary, SecondarySeq__type, SecondarySeq__AddSlv
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
  SUBROUTINE Secondaryseq__AddSlv(nunits)
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

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE SecondarySeq
