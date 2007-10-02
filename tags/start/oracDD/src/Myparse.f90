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
