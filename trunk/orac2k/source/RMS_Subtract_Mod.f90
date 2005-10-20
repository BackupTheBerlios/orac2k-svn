MODULE RMS_Subtract_Mod

!!$***********************************************************************
!!$   Time-stamp: <2005-10-20 13:55:30 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Oct 20 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*
  USE INPUT_Mod, ONLY: Read_String, Parser, err_open,err_end,err_unr&
       &,err_fnf,err_args

  TYPE subtract
     CHARACTER(8) :: target
     CHARACTER(7), DIMENSION(:), POINTER :: beta
     INTEGER :: n0
  END TYPE subtract
  TYPE rms_units
     INTEGER :: n0
     INTEGER, DIMENSION(:), POINTER :: ind
  END TYPE rms_units

  LOGICAL, SAVE :: RMS_Subtract
  INTEGER, SAVE :: n_write,krms_sub
  CHARACTER(80), SAVE :: filename
  INTEGER :: n_Res_u,n_Res
  PARAMETER(n_Res_u=25,n_Res=100)
  INTEGER, SAVE :: units
  TYPE(Subtract), DIMENSION(:), ALLOCATABLE, SAVE :: Res_u
  TYPE(rms_units), DIMENSION(:), ALLOCATABLE, SAVE :: atoms

  REAL(8), DIMENSION(:), ALLOCATABLE, SAVE :: xpt,ypt,zpt
  REAL(8), DIMENSION(:), ALLOCATABLE, SAVE :: xp_avg,yp_avg,zp_avg
  REAL(8), DIMENSION(:), ALLOCATABLE, SAVE :: xp_avg2,yp_avg2,zp_avg2
CONTAINS

!!$======================== DECLARATIONS ================================*

  SUBROUTINE Initialize(mres,nbun,res,prsymb,beta,ntap)
    
    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

    INTEGER :: mres(2,*),nbun,res(*),ntap
    CHARACTER(8) :: prsymb(*)
    CHARACTER(7) :: beta(*)

!!$----------------------- VARIABLES IN COMMON --------------------------*

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: i,k1,k2,n,nn
    INTEGER :: count,count2
    INTEGER, DIMENSION (:), ALLOCATABLE :: index
    LOGICAL, DIMENSION(:), ALLOCATABLE :: mask


!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    ALLOCATE(atoms(nbun),index(nbun),mask(ntap))

    DO k1=1,units
       count=0
       index=0
       DO n=1,nbun
          i=mres(1,n)
          IF(prsymb(res(i)) == Res_u(k1)%target) THEN
             count=count+1
             index(count)=n
          END IF
       END DO
       DO nn=1,count
          n=index(nn)             
          mask(mres(1,n):mres(2,n))=.FALSE.
          DO i=mres(1,n),mres(2,n)
             DO k2=1,Res_u(k1)%n0
                IF(Res_u(k1)%beta(k2) == beta(i)) THEN
                   mask(i)=.TRUE.
                END IF
             END DO
          END DO
          count2=0
          DO i=mres(1,n),mres(2,n)
             IF(mask(i)) count2=count2+1
          END DO
          atoms(n)%n0=count2
          ALLOCATE(atoms(n)%ind(count2))
          count2=0
          DO i=mres(1,n),mres(2,n)
             IF(mask(i)) THEN
                count2=count2+1
                atoms(n)%ind(count)=i
             END IF
          END DO
       END DO
    END DO



  END SUBROUTINE Initialize
  SUBROUTINE COMPUTE
    IMPLICIT NONE 
  END SUBROUTINE COMPUTE

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

  INCLUDE 'RMS_Subtract_Read.f90'
END MODULE RMS_Subtract_Mod
