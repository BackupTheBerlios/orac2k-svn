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
MODULE MYPARSE

!!$***********************************************************************
!!$   Time-stamp: <2007-01-04 17:46:50 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Tue Nov 14 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*
  
  USE CONSTANTS, ONLY: max_pars
  USE STRPAK
  IMPLICIT NONE
  PRIVATE
  PUBLIC Myparse_,MyParse__strngs, strngs

  CHARACTER(len=*), PARAMETER :: lst=' ,(){};'
  CHARACTER(len=max_pars), DIMENSION(:), POINTER :: Myparse__strngs=>NULL()
  CHARACTER(len=max_pars), DIMENSION(:), POINTER :: strngs=>NULL()
  
CONTAINS
  FUNCTION Myparse_(line,lst1) RESULT(out)

    INTEGER :: out
    CHARACTER(len=*) :: line
    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: lst1

    CHARACTER(len=max_pars) :: tok
    INTEGER :: count,iflag,nxt
    
    NULLIFY(strngs)

    IF(ASSOCIATED(MyParse__strngs)) DEALLOCATE(MyParse__strngs)

    count=0
    iflag=0
    nxt=1
    DO 
       IF(PRESENT(lst1)) THEN
          CALL Token(0,lst1,line,nxt,tok,iflag)
       ELSE
          CALL Token(0,lst,line,nxt,tok,iflag)
       END IF
       IF(iflag /= 0) EXIT
       count=count+1
    END DO
    out=count
    IF(out == 0) RETURN

    ALLOCATE(MyParse__strngs(count))
    count=0
    iflag=0
    nxt=1
    DO 
       IF(PRESENT(lst1)) THEN
          CALL Token(0,lst1,line,nxt,tok,iflag)
       ELSE
          CALL Token(0,lst,line,nxt,tok,iflag)
       END IF
       IF(iflag /= 0) EXIT
       count=count+1
       MyParse__strngs(count)=TRIM(tok)
    END DO
    strngs=>MyParse__Strngs

  END FUNCTION Myparse_

END MODULE MYPARSE
