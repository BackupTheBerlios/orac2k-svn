MODULE SecondarySeq

!!$***********************************************************************
!!$   Time-stamp: <2007-01-10 16:32:21 marchi>                           *
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
  PUBLIC SecondarySeq_, Secondary, SecondarySeq__type
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

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE SecondarySeq
