      SUBROUTINE RotCoordRow(n,x1,y1,z1,x0,y0,z0,xt0,yt0,zt0,Rot)

************************************************************************
*   Time-stamp: <01/04/01 17:40:03 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Mon Jan  8 2001 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER n
      REAL*8 x1(*),y1(*),z1(*),x0,y0,z0,xt0,yt0,zt0,Rot(3,3)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j
      REAL*8  xd,yd,zd,xc,yc,zc

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO i=1,n
         xd=x1(i)-x0
         yd=y1(i)-y0
         zd=z1(i)-z0
         xc=Rot(1,1)*xd+Rot(1,2)*yd+Rot(1,3)*zd
         yc=Rot(2,1)*xd+Rot(2,2)*yd+Rot(2,3)*zd
         zc=Rot(3,1)*xd+Rot(3,2)*yd+Rot(3,3)*zd
         x1(i)=xc+xt0
         y1(i)=yc+yt0
         z1(i)=zc+zt0
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
