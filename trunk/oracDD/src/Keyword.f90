MODULE Keyword

!!$***********************************************************************
!!$   Time-stamp: <2007-01-11 15:09:45 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Jan  5 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*

  USE Constants
  USE MyParse
  USE Node
  USE Errors, ONLY: Add_Errors=>Add, error_other, errmsg_w, Print_Errors
  USE STRPAK
  IMPLICIT none
  CHARACTER(len=max_char), DIMENSION(3), PARAMETER :: Tops_C=(/'PRESIDUE '&
       &,'RESI     ','PRES     '/)
  CHARACTER(len=max_char), DIMENSION(1), PARAMETER :: Tops_O=(/'RESIDUE  '/)
  CHARACTER(len=max_char), DIMENSION(7), PARAMETER :: Pars=(/'BONDS    ',&
       &'ANGLES   ','DIHEDRALS','IMPROPER ',&
       &'BOND     ','BENDINGS ','TORSION  '/)
  CHARACTER(len=max_char), PARAMETER :: Pars_nbond='NONBONDED'
  PRIVATE
  PUBLIC Keyword_, Keyword__Type, Keyword__Pop, Keyword__Size, Keyword__PopReset

  TYPE :: Keyword__Type
     INTEGER :: Begin, End
     LOGICAL :: Finish
     CHARACTER(len=max_char) :: Type,Residue,FField
  END TYPE Keyword__Type
  INTEGER, DIMENSION(:), ALLOCATABLE :: ind_x
  INTEGER, SAVE :: count0=0,ndata=0
  CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE, SAVE :: data
  LOGICAL, SAVE :: Tpg,charmm
CONTAINS
  SUBROUTINE Keyword_(data_, linea)
    CHARACTER(len=max_char), DIMENSION(:) :: data_
    CHARACTER(len=*) :: linea
    INTEGER :: n,nword,nind_x,count_a
    LOGICAL :: ok,ok_Pars,ok_Tops
    
    charmm=.TRUE.
    Tpg=.TRUE.

    ndata=SIZE(data_)
    IF(ALLOCATED(data)) DEALLOCATE(data)
    ALLOCATE(data(ndata))
    data=data_

    IF(.NOT. Node_()) STOP
    ok=.TRUE.
    ok_Pars=.FALSE.
    ok_Tops=.FALSE.
    DO n=1,ndata
       nword=MyParse_(data(n))
       IF(COUNT(TRIM(strngs(1)) == Pars) /= 0 .AND. nword == 1) THEN
          CALL Node__push(n)
          ok_Pars=.TRUE.

       ELSE IF(COUNT(TRIM(strngs(1)) == Tops_C) /= 0 &
            & .OR. COUNT(TRIM(strngs(1)) == Tops_O) /= 0 ) THEN
          CALL Node__push(n)          
          ok_Tops=.TRUE.
          IF(COUNT(TRIM(strngs(1)) == Tops_O) /= 0) charmm=.FALSE.

       ELSE IF(TRIM(strngs(1)) == TRIM(Pars_nbond) ) THEN
          CALL Node__push(n)
          ok_Pars=.TRUE.

       ELSE IF(TRIM(strngs(1)) == 'MASS' .AND. ok) THEN
          ok=.FALSE.
          CALL Node__push(n)
       END IF
    END DO
    IF(ok_Tops .AND. ok_Pars) THEN
       errmsg_w='In file ' //TRIM(Linea)//' Topology and Parameter keywords&
            & were found. Try to run further.'
       CALL Add_Errors(1,errmsg_w)
    ELSE
       IF(ok_Pars) Tpg=.FALSE.
    END IF

    nind_x=Node__Size()
    IF(ALLOCATED(ind_x)) DEALLOCATE(ind_x)
    ALLOCATE(ind_x(nind_x))
    count_a=0
    DO WHILE(Node__Pop(n))
       count_a=count_a+1
       ind_x(count_a) = n
    END DO
  END SUBROUTINE Keyword_
  FUNCTION Keyword__Size() RESULT(out)
    INTEGER :: out
    out=0
    IF(.NOT. ALLOCATED(ind_x)) RETURN
    out=SIZE(ind_x)
  END FUNCTION Keyword__Size
  SUBROUTINE Keyword__PopReset
    count0=0
  END SUBROUTINE Keyword__PopReset
  FUNCTION Keyword__Pop() RESULT(out)
    TYPE(Keyword__Type) :: out
    INTEGER :: nword

    count0=count0+1

    out % Finish=.FALSE.
    IF(.NOT. ALLOCATED(ind_x)) THEN
       out % Finish=.TRUE.
       RETURN
    END IF
    IF(count0 > SIZE(ind_x)) THEN
       count0=0
       out % Finish=.TRUE.
       RETURN
    ELSE IF(count0 == SIZE(ind_x)) THEN
       out % Begin =ind_x(count0)+1
       out % End = ndata
       nword=MyParse_(data(ind_x(count0)))
    ELSE
       out % Begin =ind_x(count0)+1
       out % End= ind_x(count0+1)-1
       nword=MyParse_(data(ind_x(count0)))
    END IF

    IF(TRIM(strngs(1)) /= 'MASS') THEN
       IF(Tpg) THEN
          out % Type = TRIM(strngs(2))
          out % Residue = TRIM(strngs(1))
       ELSE
          out % Type = ' '
          out % Residue = TRIM(strngs(1))
       END IF
       IF(charmm) THEN
          out % FField = 'CHARMM'
       ELSE
          out % FField = 'ORAC'
       END IF
    ELSE
       out % Begin =ind_x(count0)
       IF(Tpg) THEN
          out % Type = 'mass'
          out % Residue = 'HEADER'
          out % FField = 'CHARMM'
       ELSE
          out % Type = ' '
          out % Residue = 'HEADER'
          out % FField = 'CHARMM'
       END IF
    END IF
    CALL TRANLC(out % Type)
  END FUNCTION Keyword__Pop

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE Keyword
