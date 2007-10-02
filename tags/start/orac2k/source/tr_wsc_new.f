      SUBROUTINE tr_wsc(co,xp0,yp0,zp0,xp1,yp1,zp1,mass,nprot,protl)

************************************************************************
*   Time-stamp: <00/03/22 13:46:14 sterpone>                             *
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
      
      INTEGER i,j,l,m,n,count,mm
      REAL*8  xg,yg,zg,rs,tmass,xcm,ycm,zcm,xb,yb,zb,xd,yd,zd,xc,yc,zc
     &     ,xcm1,ycm1,zcm1,rs_min,xcmi,ycmi,zcmi
      INCLUDE 'pbc.h'

*----------------------- EXECUTABLE STATEMENTS ------------------------*
      
      count=0
      DO i=1,1
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
         xcmi=xcm1
         ycmi=ycm1
         zcmi=zcm1
      END DO
      count=0
      DO i=1,nprot
         mm=protl(count+1)
         xcm1=0.0D0
         ycm1=0.0D0
         zcm1=0.0D0
         DO j=1,mm
            n=protl(count+1+j)
            xp1(n)=xp0(n)-xcmi
            yp1(n)=yp0(n)-ycmi
            zp1(n)=zp0(n)-zcmi
         END DO
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
         xcm=xcm1
         ycm=ycm1
         zcm=zcm1
         rs_min=1.0D10
         IF(i .EQ. 1) THEN
            xc=co(1,1)*xcm+co(1,2)*ycm+co(1,3)*zcm
            yc=            co(2,2)*ycm+co(2,3)*zcm
            zc=                        co(3,3)*zcm
            WRITE(6,'(3f13.6)') xc,yc,zc
         END IF
         DO l=-1,1
            DO m=-1,1
               DO n=-1,1
                  xg=xcm+DFLOAT(2*l)
                  yg=ycm+DFLOAT(2*m)
                  zg=zcm+DFLOAT(2*n)
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
                  IF(i .EQ. 1) THEN
                     WRITE(6,'(f14.5,3i8)') rs,l,m,n
                  END IF
               END DO
            END DO
         END DO
         IF(i .EQ. 1) THEN
            WRITE(6,'(f14.5,3i8)') rs_min
         END IF
         DO j=1,mm
            n=protl(count+1+j)
            xp1(n)=xp1(n)-xcm1+xd
            yp1(n)=yp1(n)-ycm1+yd
            zp1(n)=zp1(n)-zcm1+zd
         END DO
         count=count+mm+1
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
