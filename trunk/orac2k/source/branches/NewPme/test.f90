PROGRAM main

!!$***********************************************************************
!!$   Time-stamp: <2009-02-27 12:05:58 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Feb  4 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program  ----*


!!$======================== DECLARATIONS ================================*

  IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

!!$----------------------- VARIABLES IN COMMON --------------------------*

!!$------------------------- LOCAL VARIABLES ----------------------------*

  Type :: My
     Integer :: i,j,k
  End type My

  Type(My) :: hg

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

  hg % i = 1;hg % j = 2;hg % k = 3;
  Call Sub(hg % i, hg % j, hg % k)
  
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*
Contains
  Subroutine sub(g % i, g % j, g % k)
    Type(My) :: g
    Write(*,*) g % i
  End Subroutine sub

END PROGRAM main
