SUBROUTINE CalcDistConf(x1,y1,z1,x2,y2,z2,w,nato,Dist)

!!$***********************************************************************
!!$   Time-stamp: <2007-10-09 17:10:56 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Sat Apr 14 2001 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*

  IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*
 
  REAL(8) :: x1(*),y1(*),z1(*),x2(*),y2(*),z2(*),w(*),Dist
  INTEGER :: nato

!!$------------------------- LOCAL VARIABLES ----------------------------*

  REAL(8), DIMENSION (:), POINTER :: X,Y,Z
  REAL(8) :: NofW
  INTEGER :: i
  REAL(8) :: x0,y0,z0

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

  NofW=1.0D0/SUM(w(1:nato))
  ALLOCATE(X(nato),Y(nato),Z(nato))
  X=x1(1:nato)-x2(1:nato)
  Y=y1(1:nato)-y2(1:nato)
  Z=z1(1:nato)-z2(1:nato)
  Dist=NofW*SUM((X*X+Y*Y+Z*Z)*w(1:nato))
  DEALLOCATE(X,Y,Z)
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE CalcDistConf
