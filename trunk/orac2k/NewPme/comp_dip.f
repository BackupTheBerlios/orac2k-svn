      SUBROUTINE comp_dip(ss_index,co,xpg,ypg,zpg,xp0,yp0,zp0,charge,dip
     &     ,ntap,ngrp,grppt)

************************************************************************
*   Time-stamp: <02/07/06 16:17:00 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Fri Sep  6 1996 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      REAL*8 dip(3,*),xp0(*),yp0(*),zp0(*),charge(*),xpg(*),ypg(*)
     &     ,zpg(*)
      INTEGER ntap,ngrp,grppt(2,*),ss_index(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8 xc,yc,zc,xd,yd,zd,xv,yv,zv,co(3,3)
      INTEGER i,j,k
      INCLUDE 'pbc.h'

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO i=1,2
         dip(1,i)=0.0D0
         dip(2,i)=0.0D0
         dip(3,i)=0.0D0
      END DO
      DO i=1,ngrp
         xd=2.0D0*PBC(xpg(i))
         yd=2.0D0*PBC(ypg(i))
         zd=2.0D0*PBC(zpg(i))
         DO j=grppt(1,i),grppt(2,i)
            xc=xp0(j)-xd
            yc=yp0(j)-yd
            zc=zp0(j)-zd
            xv=co(1,1)*xc+co(1,2)*yc+co(1,3)*zc
            yv=co(2,1)*xc+co(2,2)*yc+co(2,3)*zc
            zv=co(3,1)*xc+co(3,2)*yc+co(3,3)*zc
            k=ss_index(j)
            dip(1,k)=dip(1,k)+charge(j)*xv
            dip(2,k)=dip(2,k)+charge(j)*yv
            dip(3,k)=dip(3,k)+charge(j)*zv
         END DO
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
