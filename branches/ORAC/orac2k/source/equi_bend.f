      SUBROUTINE equi_bend(kequi,gcpu,betb,lbndg,lbend,xp0,yp0,zp0,
     &                     pota,potb,potc,potd)

      IMPLICIT none

      INTEGER lbend,kequi
      INTEGER lbndg(3,*)
      REAL*8  xp0(*),yp0(*),zp0(*),pota(*),potb(*),potc(*),potd(*),gcpu
      CHARACTER*7 betb(*)

      INTEGER i,la,lb,lc
      REAL*8  xr1,xr2,xr3,yr1,yr2,yr3,zr1,zr2,zr3,x12,x32,y12,y32,z12
     &     ,z32,rs12,rs32,xr31,yr31,zr31,rs31,rsp31,x31,y31,z31
      REAL*8  dcc2,cb,bb,pi,aux

*==================== EXECUTABLE STATEMENTS ============================

      pi=4.0D0*DATAN(1.0D0)
      aux=180.0d0/pi

      DO i=1,lbend
         la=lbndg(1,i)
         lb=lbndg(2,i)
         lc=lbndg(3,i)
         xr1=xp0(la)
         yr1=yp0(la)
         zr1=zp0(la)
         xr2=xp0(lb)
         yr2=yp0(lb)
         zr2=zp0(lb)
         xr3=xp0(lc)
         yr3=yp0(lc)
         zr3=zp0(lc)
         x12=xr1-xr2
         y12=yr1-yr2
         z12=zr1-zr2
         x32=xr3-xr2
         y32=yr3-yr2
         z32=zr3-zr2
         rs12=x12**2+y12**2+z12**2
         rs32=x32**2+y32**2+z32**2
         dcc2=DSQRT(rs12*rs32)
         cb=(x12*x32+y12*y32+z12*z32)/dcc2
         bb=DACOS(cb)
*---  Compute the Urey-Bradley term ------------------------------------
         x31=xr3-xr1
         y31=yr3-yr1
         z31=zr3-zr1
         rs31=x31**2+y31**2+z31**2
         rsp31=DSQRT(rs31)

         IF(DABS(gcpu*potc(i)) .LT. 1.0D-7) THEN
             WRITE(kequi,'(
     &         1h ,a4,1x,a4,1x,a4,1x,i5,1x,i5,1x,i5,1x,f7.2,f7.2,f7.2)')
     &         betb(la),betb(lb),betb(lc),la,lb,lc,pota(i)*gcpu,
     &         potb(i)*aux,bb*aux
         else
             WRITE(kequi,'(
     &         1h ,a4,1x,a4,1x,a4,1x,i5,1x,i5,1x,i5,1x,f7.2,f7.2,f7.2,
     &         1x,f7.2,f7.2,f7.2)')
     &         betb(la),betb(lb),betb(lc),la,lb,lc,pota(i)*gcpu,
     &         potb(i)*aux,bb*aux,potc(i)*gcpu,potd(i),rs31
         end if

      END DO
      RETURN
      END
