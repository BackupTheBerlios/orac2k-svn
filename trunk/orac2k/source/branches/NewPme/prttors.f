      SUBROUTINE prttors(step,co,xp0,yp0,zp0,tors,torsp)

************************************************************************
*                                                                      *
*     PRTPT will print the proper torsions initial angles and the      *
*     potential multiplicities.                                        *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER torsp,tors(4,*),step
      REAL*8  co(3,3),xp0(*),yp0(*),zp0(*)

*-------------------- VARIABLES IN COMMON ------------------------------

      INCLUDE 'unit.h'
      INCLUDE 'parst.h'
      REAL*4 vtor(m2)
      COMMON /rag1/ vtor

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,l1,l2,l3,l4,ntph
      REAL*8  xr1,xr2,xr3,xr4,yr1,yr2,yr3,yr4,zr1,zr2,zr3,zr4,x21,x32
     &     ,x43,y21,y32,y43,z21,z32,z43,rsq21,rsq32,rsq43,rsp21,rsp32
     &     ,rsp43
      REAL*8  cb1,cb2,cb3,sb1,sb2,sb3,coa,aux

*==================== EXECUTABLE STATEMENTS ============================

      DO i=1,torsp
         l1=tors(1,i)
         l2=tors(2,i)
         l3=tors(3,i)
         l4=tors(4,i)
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
         coa=DACOS(coa)*180.0D0/pi
         vtor(i)=coa
      END DO
      WRITE(ktop) 'T',step,torsp,(vtor(i),i=1,torsp)

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
