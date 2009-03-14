      SUBROUTINE equi_impro(kequi,gcpu,betb,litr,litor,xp0,yp0,zp0,
     &                      pota,potb,potc)

      IMPLICIT none

      CHARACTER*7 betb(*)
      INTEGER litor,litr(4,*),kequi
      REAL*8  pota(*),potb(*),potc(*),xp0(*),yp0(*),zp0(*),gcpu

      INTEGER i,l1,l2,l3,l4,naux
      REAL*8  xr1,xr2,xr3,xr4,yr1,yr2,yr3,yr4,zr1,zr2,zr3,zr4,x21,x32
     &     ,x43,y21,y32,y43,z21,z32,z43,rsq21,rsq32,rsq43,rsp21,rsp32
     &     ,rsp43
      REAL*8  cb1,cb2,cb3,sb1,sb2,sb3,aux,soa,bb
      REAL*8  coa,angle,pi

      pi=4.0d0*DATAN(1.0d0)

      DO 10 i=1,litor
          l1=litr(1,i)
          l2=litr(2,i)
          l3=litr(3,i)
          l4=litr(4,i)
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
          IF(coa.GT.1.AND.coa.LT.1.000001) coa=1.0d0
          IF(coa.LT.-1.0d0.AND.coa.GT.-1.000001) coa=-1.0d0
          bb=DACOS(coa)
          IF(bb.LT.1.0D-06) bb=1.0D-06

          angle=180.0d0

          IF(potc(i) .GT. 0.0D0) THEN

            WRITE(kequi,'(1h ,4(a4,2x),2x,4(i4,1x),3(f9.4,1x))'
     &           )betb(litr(1,i)),betb(litr(2,i)),betb(litr(3,i))
     &           ,betb(litr(4,i)),litr(1,i),litr(2,i),litr(3,i),litr(4,i
     &           ),pota(i)*gcpu,potb(i)*180.d0/pi,bb*180.d0/pi
          ELSE
            naux=DINT(potb(i)+0.5d0)
            WRITE(kequi,'(1h ,4(a4,1x),1x,4(i4,1x),f9.4,i4,f7.1,f9.3)')
     &           betb(litr(1,i)),betb(litr(2,i)),betb(litr(3,i))
     &           ,betb(litr(4,i)),litr(1,i),litr(2,i),litr(3,i),litr(4,i
     &           ),pota(i)*gcpu,naux,angle,bb*180.d0/pi
          END IF
 
10    CONTINUE


*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
