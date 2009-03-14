SUBROUTINE Make_Pol_Models(label_in,labels_out)

!!$***********************************************************************
!!$   Time-stamp: <2005-02-03 14:44:14 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Tue Jul 27 2004 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*

  IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

  CHARACTER(80) :: Label_in,Labels_out(3)

!!$----------------------- VARIABLES IN COMMON --------------------------*

  INTEGER :: start

!!$------------------------- LOCAL VARIABLES ----------------------------*

  CHARACTER(80) :: snd_part(2)

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

  snd_part(1)='no LJ'
  snd_part(2)='and LJ'
  Label_in=ADJUSTL(Label_in)
  start=INDEX(Label_in,' ')
  Labels_out(1)='Charges Gauss and LJ'
  Labels_out(2)=Label_in(1:start)//snd_part(1)
  Labels_out(3)=Label_in(1:start)//snd_part(2)

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE Make_Pol_Models
