MODULE IndPatch

!!$***********************************************************************
!!$   Time-stamp: <2007-01-10 10:03:25 marchi>                           *
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

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*

  USE Errors, ONLY: Add_Errors=>Add, errmsg_f, Print_errors
  USE SecondarySeq
  USE Parameters_Globals
  IMPLICIT none
  PRIVATE
  PUBLIC Indpatch_, Indpatch__Type, Patch

  TYPE :: Indpatch__Type
     INTEGER :: one
     INTEGER :: two
     CHARACTER(len=Max_Char) :: l1
     CHARACTER(len=Max_Char) :: l2
  END TYPE Indpatch__Type
  TYPE(Indpatch__type), DIMENSION(:), ALLOCATABLE, TARGET :: Ind_Patch
CONTAINS
  FUNCTION Indpatch_() RESULT(out)
    TYPE(Indpatch__type), DIMENSION(:), POINTER :: out
     CHARACTER(len=Max_Char) :: res_i
    INTEGER :: m,n,m1,m2,i

    out=>NULL()
    IF(.NOT. ALLOCATED(Patches)) RETURN
    m=0
    DO n=1,SIZE(patches)
       IF(patches (n) % Type == 'link') THEN
          m=m+1
       END IF
    END DO
    IF(m == 0) RETURN

    ALLOCATE(Ind_Patch(m))
    m=0
    DO n=1,SIZE(patches)
       IF(patches (n) % Type == 'link') THEN
          m=m+1
          Ind_Patch(m) % one=patches(n)%one
          Ind_Patch(m) % Two=patches(n)%Two
          Ind_Patch(m) % l1 ='Link '//TRIM(patches(n)%Res_l(1))
          Ind_Patch(m) % l2 ='Link '//TRIM(patches(n)%Res_l(2))
       END IF
    END DO
    DO i=1,SIZE(Ind_Patch)
       m1=Ind_Patch(i) % one
       m2=Ind_Patch(i) % two
       IF(m1 > SIZE(Secondary(1) % line)&
            & .OR. m2 > SIZE(Secondary(1) % line)) THEN
          errmsg_f='Link patch numbers larger than available solute residues '
          CALL Add_Errors(-1,errmsg_f)
       END IF
    END DO
    CALL Print_Errors()
    n=1
    DO m=1,SIZE(Secondary(n) % line)
       res_i=Secondary(n) % line(m)
       DO i=1,SIZE(Ind_Patch)
          IF(Ind_Patch(i) % one == m) THEN
             res_i=Ind_Patch(i) % l1 
             Secondary(n) % line(m)=TRIM(res_i)
          ELSE IF(Ind_Patch(i) % Two == m) THEN
             res_i=Ind_Patch(i) % l2
             Secondary(n) % line(m)=TRIM(res_i)
          END IF
       END DO
    END DO
    out=>Ind_Patch
  END FUNCTION Indpatch_
  FUNCTION Indpatch__Assign() RESULT(out)
    TYPE(Indpatch__type), DIMENSION(:), POINTER :: out
    out=>NULL()
    IF(.NOT. ALLOCATED(Ind_Patch)) RETURN
    out=>Ind_Patch
  END FUNCTION Indpatch__Assign

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE IndPatch
