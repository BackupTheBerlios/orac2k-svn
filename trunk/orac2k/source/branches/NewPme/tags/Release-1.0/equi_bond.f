      SUBROUTINE equi_bond(kequi,gcpu,betb,lbnd,lbond,xp0,yp0,zp0,
     &                     pota,potb)
      IMPLICIT none

      CHARACTER*7 betb(*)
      INTEGER lbond,kequi
      INTEGER lbnd(2,*)
      REAL*8  xp0(*),yp0(*),zp0(*),pota(*),potb(*)

      INTEGER i,la,lb
      REAL*8  xr1,xr2,yr1,yr2,zr1,zr2,x21,y21,z21,rs21, gcpu
    
      DO i=1,lbond
          la=lbnd(1,i)
          lb=lbnd(2,i)
          xr1=xp0(la)
          yr1=yp0(la)
          zr1=zp0(la)
          xr2=xp0(lb)
          yr2=yp0(lb)
          zr2=zp0(lb)
          x21=xr2-xr1
          y21=yr2-yr1
          z21=zr2-zr1
          rs21=DSQRT(x21**2+y21**2+z21**2)
          WRITE(kequi
     &         ,'(1h ,a4,1x,a4,1x,i5,1x,i5,1x,f7.2,f8.3,f8.3)'
     &         )betb(lbnd(1,i)),betb(lbnd(2,i)),lbnd(1,i),lbnd(2,i)
     &         ,pota(i)*gcpu,potb(i), rs21
 
      END DO
      RETURN
      END
