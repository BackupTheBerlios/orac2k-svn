MODULE Module_Topology
  USE Class_Bending, Bending_Init=> Init
  USE Class_Connect, Connect_Init=> Init, Connect_Destroy=>Destroy,&
       & Connect_Print=> Print 
  REAL, DIMENSION (:,:), ALLOCATABLE, SAVE :: Angle
  TYPE(Bending) :: bndg
CONTAINS
  SUBROUTINE Calc_Topology(concta,m,nato)

!!$***********************************************************************
!!$   Time-stamp: <2002-09-27 10:49:30 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Tue Feb 22 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*

    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

    INTEGER :: m,nato
    INTEGER :: concta(m,*)

!!$------------------------- LOCAL VARIABLES ----------------------------*

    TYPE(Connect), DIMENSION(:), POINTER :: connect
    INTEGER :: langle

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    connect=Connect_Init(concta,m,nato)
    bndg=Bending_Init(connect)
    
    langle=SIZE(bndg % atm,2)
    ALLOCATE(angle(3,langle))
    angle(:,:)=bndg % atm(:,:)
    DEALLOCATE(bndg % atm)

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

  END SUBROUTINE Calc_Topology
END MODULE Module_Topology
