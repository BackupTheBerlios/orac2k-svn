      SUBROUTINE tr_wsc(co,xp0,yp0,zp0,xp1,yp1,zp1,mass,nprot,protl)

************************************************************************
*   Time-stamp: <2005-02-19 16:12:39 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Jul 20 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nprot,protl(*)
      REAL*8  xp0(*),yp0(*),zp0(*),xp1(*),yp1(*),zp1(*),mass(*),co(3,3)

*------------------------- LOCAL VARIABLES ----------------------------*
      
      INTEGER i,j,l,m,n,count,mm,ip_mass
      REAL*8  xg,yg,zg,rs,tmass,xcm,ycm,zcm,xb,yb,zb,xd,yd,zd,xc,yc,zc
     &     ,xcm1,ycm1,zcm1,rs_min,xcmi,ycmi,zcmi,max_mass
      INCLUDE 'pbc.h'

*----------------------- EXECUTABLE STATEMENTS ------------------------*

c$$$-------------------------------------------------------------------*
c$$$-- Find heaviest molecule and use it as the center of the ---------*
c$$$-- conformation                                           ---------*
c$$$-------------------------------------------------------------------*
      
      ip_mass=0
      max_mass=-1.0D0
      count=0
      DO i=1,nprot
         mm=protl(count+1)
         xcm1=0.0D0
         ycm1=0.0D0
         zcm1=0.0D0
         tmass=0.0D0
         DO j=1,mm
            n=protl(count+1+j)
            tmass=tmass+mass(j)
         END DO
         IF(tmass .GT. max_mass) THEN
            max_mass=tmass
            ip_mass=i
         END IF
         count=count+mm+1
      END DO

      count=0
      DO i=1,ip_mass
         mm=protl(count+1)
         IF(i .EQ. ip_mass) THEN
            xcm1=0.0D0
            ycm1=0.0D0
            zcm1=0.0D0
            tmass=0.0D0
            DO j=1,mm
               n=protl(count+1+j)
               xcm1=xcm1+xp0(n)*mass(n)
               ycm1=ycm1+yp0(n)*mass(n)
               zcm1=zcm1+zp0(n)*mass(n)
               tmass=tmass+mass(n)
            END DO
            xcm1=xcm1/tmass
            ycm1=ycm1/tmass
            zcm1=zcm1/tmass
            xcmi=xcm1-2.0D0*PBC(xcm1)
            ycmi=ycm1-2.0D0*PBC(ycm1)
            zcmi=zcm1-2.0D0*PBC(zcm1)
         END IF
         count=count+mm+1
      END DO

      count=0
      DO i=1,nprot
         mm=protl(count+1)
         xcm1=0.0D0
         ycm1=0.0D0
         zcm1=0.0D0
         tmass=0.0D0
         DO j=1,mm
            n=protl(count+1+j)
            xcm1=xcm1+xp0(n)*mass(n)
            ycm1=ycm1+yp0(n)*mass(n)
            zcm1=zcm1+zp0(n)*mass(n)
            tmass=tmass+mass(n)
         END DO
         xcm1=xcm1/tmass
         ycm1=ycm1/tmass
         zcm1=zcm1/tmass
         xcm=xcm1-2.0D0*PBC(xcm1)
         ycm=ycm1-2.0D0*PBC(ycm1)
         zcm=zcm1-2.0D0*PBC(zcm1)
         rs_min=1.0D10
         DO l=-1,1
            DO m=-1,1
               DO n=-1,1
                  xg=xcm-xcmi+DFLOAT(2*l)
                  yg=ycm-ycmi+DFLOAT(2*m)
                  zg=zcm-zcmi+DFLOAT(2*n)
                  xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
                  yc=           co(2,2)*yg+co(2,3)*zg
                  zc=                      co(3,3)*zg
                  rs=DSQRT(xc**2+yc**2+zc**2)
                  IF(rs .LT. rs_min) THEN
                     rs_min=rs
                     xd=xg
                     yd=yg
                     zd=zg
                  END IF
               END DO
            END DO
         END DO
         DO j=1,mm
            n=protl(count+1+j)
            xp1(n)=xp0(n)-xcm1+xd-xcmi
            yp1(n)=yp0(n)-ycm1+yd-ycmi
            zp1(n)=zp0(n)-zcm1+zd-zcmi
         END DO
         count=count+mm+1
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
