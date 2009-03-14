MODULE CENTER_SOL_Mod

!!$***********************************************************************
!!$   Time-stamp: <2005-10-15 21:52:05 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Sat Oct 15 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

  LOGICAL, SAVE :: Center_solute=.FALSE.

CONTAINS
!!$---- This subroutine is part of the program ORAC ----*
  SUBROUTINE Compute(xp0,yp0,zp0,xpa,ypa,co,oc,mass,ntap,ss_index)

!!$======================== DECLARATIONS ================================*

    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*
    
    INTEGER :: ntap,ss_index(*)
    REAL(8) :: xpa(*),ypa(*),zpa(*),xp0(*),yp0(*),zp0(*),co(3,3),oc(3&
         &,3),mass(*)

!!$------------------------- LOCAL VARIABLES ----------------------------*

    REAL(8) :: xcm,ycm,zcm,totmass

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

                  
    CALL change_frame(co,oc,-1,ntap,xp0,yp0,zp0,xpa,ypa,zpa)

    xcm=0.0D0
    ycm=0.0D0
    zcm=0.0D0
    totmass=0.0D0
    DO i=1,ntap
       IF(ss_index(i) == 1) THEN
          xcm=xcm+mass(i)*xpa(i)
          ycm=ycm+mass(i)*ypa(i)
          zcm=zcm+mass(i)*zpa(i)
          totmass=totmass_mass(i)
       END IF
    END DO
    xcm=xcm/totmass
    ycm=ycm/totmass
    zcm=zcm/totmass
    xpa(1:ntap)=xpa(1:ntap)-xcm
    ypa(1:ntap)=ypa(1:ntap)-ycm
    zpa(1:ntap)=zpa(1:ntap)-zcm
    CALL change_frame(co,oc,1,ntap,xpa,ypa,zpa,xp0,yp0,zp0)

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*
  END SUBROUTINE Compute
END MODULE CENTER_SOL_Mod
