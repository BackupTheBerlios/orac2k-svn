MODULE CENTER_SOL_Mod

!!$***********************************************************************
!!$   Time-stamp: <2008-03-20 14:19:46 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Sat Oct 15 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

  LOGICAL, SAVE :: Center_object=.FALSE.
  CHARACTER(7), SAVE :: target
  CHARACTER(3), DIMENSION(22), SAVE :: protein=(/ 'ala','arg','asn'&
       &,'asp','cys','gln','glu','gly','hsd','hse','hsp','ile'&
       &,'leu','lys','met','phe','pro','ser','thr','trp','tyr'&
       &,'val'/)
  CHARACTER(3), DIMENSION(3), SAVE :: water=(/'hoh','spc','tip'/)
  LOGICAL, DIMENSION(:), ALLOCATABLE, SAVE :: mask
CONTAINS
  SUBROUTINE Initialize(prsymb,beta,res,ss_index,ntap)

!!$======================== DECLARATIONS ================================*

    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*
    
    INTEGER :: ntap,ss_index(*),res(*)
    CHARACTER(8) :: prsymb(*)
    CHARACTER(7) :: beta(*)

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: i,j
    CHARACTER(3) :: Residue

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    ALLOCATE(mask(ntap))
    mask=.FALSE.
    SELECT CASE (target)
    CASE('WATER')
       DO i=1,ntap
          Residue=prsymb(res(i))(1:3)
          DO j=1,LEN(water)
             IF(Residue == water(j)) mask(i)=.TRUE.
             IF(Residue == upper(water(j))) mask(i)=.TRUE.
          END DO
       END DO
    CASE('PROTEIN')
       DO i=1,ntap
          Residue=prsymb(res(i))(1:3)
          DO j=1,LEN(protein)
             IF(Residue == protein(j)) mask(i)=.TRUE.
             IF(Residue == upper(protein(j))) mask(i)=.TRUE.
          END DO
       END DO
    CASE('SOLUTE')
       DO i=1,ntap
          IF(ss_index(i) == 1) mask(i)=.TRUE.
       END DO
    CASE DEFAULT
       mask=.TRUE.
    END SELECT

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

  END SUBROUTINE Initialize
  SUBROUTINE Compute(xp0,yp0,zp0,mass,ntap)

!!$======================== DECLARATIONS ================================*

    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*
    
    INTEGER :: ntap
    REAL(8) :: xp0(*),yp0(*),zp0(*),co(3,3),mass(*)

!!$------------------------- LOCAL VARIABLES ----------------------------*

    REAL(8) :: xcm,ycm,zcm,totmass
    INTEGER :: i

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

                  
    xcm=0.0D0
    ycm=0.0D0
    zcm=0.0D0
    totmass=0.0D0
    DO i=1,ntap
       IF(mask(i)) THEN
          xcm=xcm+mass(i)*xp0(i)
          ycm=ycm+mass(i)*yp0(i)
          zcm=zcm+mass(i)*zp0(i)
          totmass=totmass+mass(i)
       END IF
    END DO
    xcm=xcm/totmass
    ycm=ycm/totmass
    zcm=zcm/totmass
    DO i=1,ntap
       xp0(i)=xp0(i)-xcm
       yp0(i)=yp0(i)-ycm
       zp0(i)=zp0(i)-zcm
    END DO
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

  END SUBROUTINE Compute
  FUNCTION upper(string) RESULT(out)
    IMPLICIT NONE 
    CHARACTER(3) :: string,out
    INTEGER :: upper_to_lower,len_string,i

    upper_to_lower = IACHAR("a") - IACHAR("A")
    out=string
    len_string=LEN_TRIM(out)
    DO i=1,len_string
       SELECT CASE (out(i:i))
       CASE ('a':'z')
          out(i:i) = ACHAR(IACHAR(out(i:i)) - upper_to_lower)
       END SELECT
    END DO
  END FUNCTION upper
END MODULE CENTER_SOL_Mod
