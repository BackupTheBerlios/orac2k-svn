MODULE MD_Mod

!!$***********************************************************************
!!$   Time-stamp: <2006-02-13 14:41:03 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Feb 10 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************
!!$---- This Module is part of the program ORAC ----*

  
  INTEGER, DIMENSION (:), ALLOCATABLE :: indxi,indxj,indxk
  REAL(8), DIMENSION (:), ALLOCATABLE ::  fpx_n0,fpy_n0,fpz_n0,fcax_n0&
       &,fcay_n0,fcaz_n0,fpx_n1,fpy_n1,fpz_n1,fcax_n1,fcay_n1,fcaz_n1&
       &,fpx_m,fpy_m,fpz_m,fcax_m,fcay_m,fcaz_m,fpx_l,fpy_l,fpz_l&
       &,fcax_l,fpx_h,fpy_h,fpz_h,fcax_h,fcay_h,fcaz_h,fpx,fpy,fpz
  REAL(8), DIMENSION (:), ALLOCATABLE::  vpx,vpy,vpz,vpx1,vpy1,vpz1,vcax&
       &,vcay,vcaz,vcbx,vcby,vcbz 
  REAL(8), DIMENSION (:), ALLOCATABLE ::  xpo,ypo,zpo,xpa,ypa,zpa&
       &,xpga,ypga,zpga
  REAL(8) :: coo(3,3),co2(3,3),oc2(3,3)

CONTAINS

!!$======================== DECLARATIONS ================================*

  SUBROUTINE Init(mapnl)
    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

!!$----------------------- VARIABLES IN COMMON --------------------------*

!!$------------------------- LOCAL VARIABLES ----------------------------*


!!$----------------------- EXECUTABLE STATEMENTS ------------------------*






!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*


!!$======================== DECLARATIONS ================================*

  END SUBROUTINE Init

  SUBROUTINE Simulate(xp0,yp0,zp0)
    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

!!$----------------------- VARIABLES IN COMMON --------------------------*

!!$------------------------- LOCAL VARIABLES ----------------------------*


!!$----------------------- EXECUTABLE STATEMENTS ------------------------*





!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*
  END SUBROUTINE Simulate
END MODULE MD_Mod
