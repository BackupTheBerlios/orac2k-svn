PROGRAM main

!!$***********************************************************************
!!$   Time-stamp: <2005-02-04 14:19:25 marchi>                           *
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

  USE Neighbors_Mod
  IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

!!$----------------------- VARIABLES IN COMMON --------------------------*

!!$------------------------- LOCAL VARIABLES ----------------------------*


  INTEGER :: o,n,m,count,i,j
  INTEGER, DIMENSION (120000) :: nnlpp
  REAL(8) :: duni
  TYPE(Neighbors), DIMENSION (:), POINTER :: nn_Pol

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*
  
  o=100
  n=20
  count=0
  DO i=1,n
     m=INT(DFLOAT(o)*duni())
     nnlpp(count+1)=m
     DO j=1,m
        nnlpp(count+1+j)=INT(DFLOAT(o)*duni())
     END DO
     IF(i==5) WRITE(*,*) (nnlpp(count+1+j),j=1,4)
     count=count+1+m
  END DO
  CALL Neighbors_init(nnlpp,1,n)
  nn_Pol=>nn

  WRITE(*,*) nn(5)%Neigh(1:4)

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END PROGRAM main
