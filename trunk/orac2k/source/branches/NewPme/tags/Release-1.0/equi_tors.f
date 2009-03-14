      SUBROUTINE equi_tors(kequi,gcpu,betb,ltor,ltors,xp0,yp0,zp0,
     &                  pota,potb)

      IMPLICIT none

      CHARACTER*7 betb(*)
      INTEGER ltors,ltor(4,*),kequi
      REAL*8  pota(*),xp0(*),yp0(*),zp0(*),gcpu,potb(*)

      INTEGER i,l1,l2,l3,l4
      REAL*8  xr1,xr2,xr3,xr4,yr1,yr2,yr3,yr4,zr1,zr2,zr3,zr4,x21,x32
     &     ,x43,y21,y32,y43,z21,z32,z43,rsq21,rsq32,rsq43,rsp21,rsp32
     &     ,rsp43
      REAL*8  cb1,cb2,cb3,sb1,sb2,sb3,aux
      REAL*8  bb,coa,pi,angle,aux1,aux2

      pi = 4.0d0*DATAN(1.0d0)

      DO 10 i=1,ltors
          l1=ltor(1,i)
          l2=ltor(2,i)
          l3=ltor(3,i)
          l4=ltor(4,i)

          xr1=xp0(l1)
          yr1=yp0(l1)
          zr1=zp0(l1)
          xr2=xp0(l2)
          yr2=yp0(l2)
          zr2=zp0(l2)
          xr3=xp0(l3)
          yr3=yp0(l3)
          zr3=zp0(l3)
          xr4=xp0(l4)
          yr4=yp0(l4)
          zr4=zp0(l4)
          x21=xr2-xr1
          y21=yr2-yr1
          z21=zr2-zr1
          x32=xr3-xr2
          y32=yr3-yr2
          z32=zr3-zr2
          x43=xr4-xr3
          y43=yr4-yr3
          z43=zr4-zr3
          rsq21=x21**2+y21**2+z21**2
          rsq32=x32**2+y32**2+z32**2
          rsq43=x43**2+y43**2+z43**2
          rsp21=DSQRT(rsq21)
          rsp32=DSQRT(rsq32)
          rsp43=DSQRT(rsq43)
          cb1=(x21*x32+y21*y32+z21*z32)/(rsp21*rsp32)
          cb2=(x43*x32+y43*y32+z43*z32)/(rsp43*rsp32)
          cb3=(x21*x43+y21*y43+z21*z43)/(rsp21*rsp43)
          sb1=DSQRT(1.0d0-cb1**2)
          sb2=DSQRT(1.0d0-cb2**2)
          sb3=DSQRT(1.0d0-cb3**2)
          aux=sb1*sb2
          coa=(cb1*cb2-cb3)/aux
          IF(coa.GT.1.0d0.AND.coa.LT.1.000001) coa=1.0d0
          IF(coa.LT.-1.0d0.AND.coa.GT.-1.000001) coa=-1.0d0
          bb=DACOS(coa)
     
          if ( pota(i).lt.0.0d0 ) then
               angle=180.0
          else
               angle=0.0d0
          endif
c          if ( dabs(gcpu*pota(i)).lt. 1.0D-7) go to 10
          aux1 = dabs(pota(i)) * potb(i)/DABS(potb(i))
          aux2 = dabs(potb(i))
          WRITE(kequi,'(1h ,4(a4,1x),1x,4(i4,1x),
     x                 f9.4,1x,2(f5.1,1x),f8.3)')
     x         betb(ltor(1,i)),betb(ltor(2,i)),betb(ltor(3,i)),
     x         betb(ltor(4,i)),ltor(1,i),ltor(2,i),ltor(3,i),
     x         ltor(4,i),aux1*gcpu,aux2,angle,bb*180/pi 

10    CONTINUE
      RETURN
      END
