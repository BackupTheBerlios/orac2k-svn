MODULE Atoms

!!$***********************************************************************
!!$   Time-stamp: <2007-01-04 15:04:41 marchi>                           *
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

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*

  USE CONSTANTS
  USE AtomCnt, ONLY: AtomCnt, atm
  IMPLICIT none
  TYPE :: Atom
     REAL(8) :: x=Huge,y=Huge,z=Huge
     TYPE(AtomCnt), POINTER :: c
  END TYPE Atom
  TYPE(Atom), DIMENSION(:), ALLOCATABLE, SAVE :: AtomsTot, AtomSlv, AtomSlt
  INTEGER, SAVE :: n_Slt,n_Slv

CONTAINS
  SUBROUTINE Init
    INTEGER :: n,m,m1,m2

    m=0
    DO n=1,SIZE(atm)
       IF(atm(n) % Id_Slv == 1) m=m+1
    END DO
    n_Slt=m
    n_Slv=SIZE(atm)-n_Slt
    ALLOCATE(AtomSlt(n_Slt),AtomSlv(n_Slv))
    m1=0
    m2=0
    DO n=1,SIZE(atm)
       IF(atm(n) % Id_Slv == 1) THEN
          m1=m1+1
          AtomSlt(m1) % c => atm(n)
       ELSE 
          m2=m2+1
          AtomSlv(m2) % c => atm(n)
       END IF
    END DO
  END SUBROUTINE Init

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE Atoms
