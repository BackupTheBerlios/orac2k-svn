SUBROUTINE ResetCofM(xp0,yp0,zp0,nato,mass)

!!$***********************************************************************
!!$   Time-stamp: <01/04/13 12:45:10 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Apr 12 2001 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*

  IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

  REAL(8) :: xp0(*),yp0(*),zp0(*),mass(*)
  INTEGER :: nato

!!$------------------------- LOCAL VARIABLES ----------------------------*

  INTEGER :: i
  REAL(8) :: xpcm,ypcm,zpcm,TMass

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

  Tmass=0.0D0
  xpcm=0.0D0
  ypcm=0.0D0
  zpcm=0.0D0
  DO i=1,nato
     IF(mass(i) > 2.0D0) THEN
        Tmass=Tmass+mass(i)
        xpcm=xpcm+mass(i)*xp0(i)
        ypcm=ypcm+mass(i)*yp0(i)
        zpcm=zpcm+mass(i)*zp0(i)
     END IF
  END DO
  xpcm=xpcm/Tmass
  ypcm=ypcm/Tmass
  zpcm=zpcm/Tmass
  
  xp0(1:nato)=xp0(1:nato)-xpcm
  yp0(1:nato)=yp0(1:nato)-ypcm
  zp0(1:nato)=zp0(1:nato)-zpcm

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE ResetCofM
