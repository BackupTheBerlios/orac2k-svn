MODULE MYPARSE_Mod

!!$***********************************************************************
!!$   Time-stamp: <2006-11-24 11:46:35 marchi>                           *
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
  CHARACTER(len=*), PARAMETER :: lst=' ,(){};'

!!$----------------------------- ARGUMENTS ------------------------------*

!!$----------------------- VARIABLES IN COMMON --------------------------*

!!$------------------------- LOCAL VARIABLES ----------------------------*


!!$----------------------- EXECUTABLE STATEMENTS ------------------------*
  
CONTAINS
  SUBROUTINE parse(line,strngs,lst1)
    IMPLICIT NONE
    CHARACTER(len=*) :: line
    CHARACTER(len=max_pars), DIMENSION(:), ALLOCATABLE :: strngs
    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: lst1

    CHARACTER(len=max_pars) :: tok
    INTEGER :: count,iflag,nxt
    
    IF(ALLOCATED(strngs)) DEALLOCATE(strngs)

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
    ALLOCATE(strngs(count))
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
       strngs(count)=TRIM(tok)
    END DO
    
  END SUBROUTINE parse

END MODULE MYPARSE_Mod
