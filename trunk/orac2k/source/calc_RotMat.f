      SUBROUTINE calc_RotMat(wca,xp_ini,yp_ini,zp_ini,xp0,yp0,zp0,xpcc
     &     ,ypcc,zpcc,RotMat,xt_cm,yt_cm,zt_cm,nato_slt)

************************************************************************
*   Time-stamp: <2007-10-09 17:02:14 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Jun 29 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none
      INCLUDE  'parst.h'

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8  wca(*),xp0(*),yp0(*),zp0(*),xp_ini(*),yp_ini(*)
     &     ,zp_ini(*),xpcc,ypcc,zpcc,RotMat(3,3),xt_cm,yt_cm,zt_cm
      INTEGER nato_slt

*----------------------- VARIABLES IN COMMON --------------------------*

      REAL*8  xyz(3,m1),xyz0(3,m1),xyzfit(3,m1),wca2(m1),xpo(m1
     &     ),ypo(m1),zpo(m1),dcm(3)

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8 err,sum,qt(4),q(0:7),d(3,3),xd,yd,zd,xpp,ypp,zpp,xcm,ycm
     &     ,zcm
      INTEGER n,m,i,j,iret
      CHARACTER*80 msg

*----------------------- EXECUTABLE STATEMENTS ------------------------*

*=======================================================================
*----- Calculate XRMS from template ------------------------------------
*=======================================================================

      m=nato_slt
      DO n=1,m
         xyz(1,n)=xp0(n)
         xyz(2,n)=yp0(n)
         xyz(3,n)=zp0(n)
         xyz0(1,n)=xp_ini(n)
         xyz0(2,n)=yp_ini(n)
         xyz0(3,n)=zp_ini(n)
         wca2(n)=wca(n)
      END DO
      sum=0.0D0
      DO i=1,m
         sum=sum+wca2(i)
      END DO
      IF(sum .LT. 4.0D0) THEN
         DO n=1,m
            wca2(n)=1.0D0
         END DO
      END IF
         
      CALL normal(wca2,m)
      CALL rigfit(0,m,xyz0,xyz,wca2,wca2,q,dcm,xyzfit,err,iret,msg)

      d(1,1)=-2.d0*q(2)**2-2.d0*q(3)**2+1.d0
      d(1,2)=2.d0*(-q(0)*q(3)+q(1)*q(2))
      d(1,3)=2.d0*(q(0)*q(2)+q(1)*q(3))
      d(2,1)=2.d0*(q(0)*q(3)+q(1)*q(2))
      d(2,2)=-2.d0*q(1)**2-2.d0*q(3)**2+1.d0
      d(2,3)=2.d0*(-q(0)*q(1)+q(2)*q(3))
      d(3,1)=2.d0*(-q(0)*q(2)+q(1)*q(3))
      d(3,2)=2.d0*(q(0)*q(1)+q(2)*q(3))
      d(3,3)=-2.d0*q(1)**2-2.d0*q(2)**2+1.d0
      DO i=1,3
         DO j=1,3
            RotMat(i,j)=d(i,j)
         END DO
      END DO
      xcm=0.0D0
      ycm=0.0D0
      zcm=0.0D0
      DO i=1,m
         xcm=xcm+wca2(i)*xyz0(1,i)
         ycm=ycm+wca2(i)*xyz0(2,i)
         zcm=zcm+wca2(i)*xyz0(3,i)
      END DO
      xt_cm=xcm
      yt_cm=ycm
      zt_cm=zcm
      xpcc=xt_cm-dcm(1)
      ypcc=yt_cm-dcm(2)
      zpcc=zt_cm-dcm(3)

c$$$      DO i=1,m
c$$$         xpp=xyz(1,i)-xpcc
c$$$         ypp=xyz(2,i)-ypcc
c$$$         zpp=xyz(3,i)-zpcc
c$$$         xd=d(1,1)*xpp+d(1,2)*ypp+d(1,3)*zpp
c$$$         yd=d(2,1)*xpp+d(2,2)*ypp+d(2,3)*zpp
c$$$         zd=d(3,1)*xpp+d(3,2)*ypp+d(3,3)*zpp
c$$$         xd=xd+xt_cm
c$$$         yd=yd+yt_cm
c$$$         zd=zd+zt_cm
c$$$         IF(i .LT. m+1) THEN
c$$$            WRITE(76,'(i8,3f12.5)') i,xd,yd,zd
c$$$            WRITE(76,'(i8,3f12.5)') i,xyzfit(1,i),xyzfit(2,i),xyzfit(3,i
c$$$     &           )
c$$$            WRITE(76,'(i8,3f12.5)') i,xyz0(1,i),xyz0(2,i),xyz0(3,i)
c$$$         END IF
c$$$      END DO
c$$$      STOP
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*
         
      RETURN
      END

