      SUBROUTINE change_frame(co,oc,way,nato,xp0,yp0,zp0,xp1,yp1,zp1)

************************************************************************
*   Time-stamp: <97/02/28 12:02:11 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Fri Feb 28 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nato,way
      REAL*8  xp0(*),yp0(*),zp0(*),xp1(*),yp1(*),zp1(*),co(3,3),oc(3,3)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i
      REAL*8  xpi,ypi,zpi

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      IF(way .EQ. 1) THEN
         DO i=1,nato
            xpi=co(1,1)*xp0(i)+co(1,2)*yp0(i)+co(1,3)*zp0(i)
            ypi=co(2,1)*xp0(i)+co(2,2)*yp0(i)+co(2,3)*zp0(i)
            zpi=co(3,1)*xp0(i)+co(3,2)*yp0(i)+co(3,3)*zp0(i)
            xp1(i)=xpi
            yp1(i)=ypi
            zp1(i)=zpi
         END DO
      ELSE IF(way .EQ. -1) THEN
         DO i=1,nato
            xpi=oc(1,1)*xp0(i)+oc(1,2)*yp0(i)+oc(1,3)*zp0(i)
            ypi=oc(2,1)*xp0(i)+oc(2,2)*yp0(i)+oc(2,3)*zp0(i)
            zpi=oc(3,1)*xp0(i)+oc(3,2)*yp0(i)+oc(3,3)*zp0(i)
            xp1(i)=xpi
            yp1(i)=ypi
            zp1(i)=zpi
         END DO
      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
