MODULE NUMERICS_Mod

!!$***********************************************************************
!!$   Time-stamp: <2006-12-20 12:08:36 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Wed Dec 20 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC  ----*
  IMPLICIT none
  PRIVATE
  PUBLIC :: MatInv, Determinant
  REAL(8), SAVE :: Determinant
CONTAINS
  SUBROUTINE MatInv(A,B)
    REAL(8), DIMENSION(:,:) :: A
    REAL(8), DIMENSION(:,:) :: B

    INTEGER, DIMENSION(:), ALLOCATABLE :: l
    INTEGER, DIMENSION(:), ALLOCATABLE :: m
    INTEGER :: o,p
    REAL(8) :: d
    
    o=SIZE(A,1)
    p=SIZE(A,2)
    ALLOCATE(l(o),m(o))
    b=a
    IF(o < 20) THEN
       CALL DMINV(b,o,d,l,m)
    END IF
    IF(d == 0.0D0) b=0.0D0
    Determinant=d
    DEALLOCATE(l,m)
  CONTAINS
    INCLUDE 'dminv.f'
  END SUBROUTINE MatInv
END MODULE NUMERICS_Mod
