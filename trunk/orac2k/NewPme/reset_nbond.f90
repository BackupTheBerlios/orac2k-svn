SUBROUTINE reset_nbond(nato,beta,atomg,ss_index,Boltz_Group)

!!$***********************************************************************
!!$   Time-stamp: <2005-02-14 17:02:55 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Mon Feb 14 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*

  IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*
  
  INTEGER :: nato,atomg(*),ss_index(*)
  CHARACTER(7) :: beta(*)
  INTEGER :: Boltz_Group(*)

!!$----------------------- VARIABLES IN COMMON --------------------------*

!!$------------------------- LOCAL VARIABLES ----------------------------*

  INTEGER :: i,ig

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

  DO i=1,nato
     ig=atomg(i)
     Boltz_Group(ig)=ss_index(i)
     IF(beta(i) == 'boltz') THEN
        Boltz_Group(ig)=3
     END IF
  END DO

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE reset_nbond
