SUBROUTINE DCopy_V3(nstart,nlocal,x1,y1,z1,x,y,z)

!!$***********************************************************************
!!$   Time-stamp: <2004-12-10 09:00:53 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Nov 25 2004 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*

  IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

  INTEGER :: nstart,nlocal
  REAL(8) :: x1(*),y1(*),z1(*),x(*),y(*),z(*)

!!$----------------------- VARIABLES IN COMMON --------------------------*

!!$------------------------- LOCAL VARIABLES ----------------------------*


!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

  CALL dcopy(nlocal,x1(nstart),1,x(nstart),1)
  CALL dcopy(nlocal,y1(nstart),1,y(nstart),1)
  CALL dcopy(nlocal,z1(nstart),1,z(nstart),1)



!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE DCopy_V3
