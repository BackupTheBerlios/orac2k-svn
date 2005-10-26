      SUBROUTINE trans_center(nprot,nato,ss_point,xp0,yp0,zp0,xpcm,ypcm
     &     ,zpcm)

************************************************************************
*   Time-stamp: <97/07/03 18:02:55 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Jul  3 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nato,nprot,ss_point(*)
      REAL*8  co(3,3),xp0(*),yp0(*),zp0(*),xpcm(*),ypcm(*),zpcm(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER j,m,i
      REAL*8  sumx,sumy,sumz

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      m=ss_point(1)
      sumx=0.0D0
      sumy=0.0D0
      sumz=0.0D0
      DO j=1,m
         i=ss_point(j+1)
         sumx=sumx+xp0(i)
         sumy=sumy+yp0(i)
         sumz=sumz+zp0(i)
      END DO
      sumx=sumx/DBLE(m)
      sumy=sumy/DBLE(m)
      sumz=sumz/DBLE(m)
      DO i=1,nato
         xp0(i)=xp0(i)-sumx
         yp0(i)=yp0(i)-sumy
         zp0(i)=zp0(i)-sumz
      END DO
      DO i=1,nprot
         xpcm(i)=xpcm(i)-sumx
         ypcm(i)=ypcm(i)-sumy
         zpcm(i)=zpcm(i)-sumz
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
