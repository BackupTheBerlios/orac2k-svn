      SUBROUTINE comp_abmd_tors(rject,iflag,enout,alpha,litr,litor,xp0
     &     ,yp0,zp0,fpx,fpy,fpz,uconf,angle)

************************************************************************
*   Time-stamp: <97/09/15 08:36:42 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Jul  8 1997 -                                     *
*                                                                      *
************************************************************************


*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER litor,litr(4,*),iflag
      REAL*8  enout,alpha,uconf,angle,xp0(*),yp0(*),zp0(*),fpx(*),fpy(*)
     &     ,fpz(*)
      LOGICAL*4 rject

*-------------------- VARIABLES IN COMMON ------------------------------

      INCLUDE 'unit.h'

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,l1,l2,l3,l4,ntph,tmp,type,nwind
      REAL*8  xr1,xr2,xr3,xr4,yr1,yr2,yr3,yr4,zr1,zr2,zr3,zr4,x21,x32
     &     ,x43,y21,y32,y43,z21,z32,z43,rsq21,rsq32,rsq43,rsp21,rsp32
     &     ,rsp43
      REAL*8  cb1,cb2,cb3,sb1,sb2,sb3,aux,aux1,aux2,auxa,auxb,auxc,dbx11
     &     ,dbx31,dbx12,dbx22,dbx32,dbx13,dbx23,dbx33,dbx24,dbx34,dby11
     &     ,dby31,dby12,dby22,dby32,dby13,dby23,dby33,dby24,dby34,dbz11
     &     ,dbz31,dbz12,dbz22,dbz32,dbz13,dbz23,dbz33,dbz24,dbz34,uux1
     &     ,uux2,uux3,uux4,uuy1,uuy2,uuy3,uuy4,uuz1,uuz2,uuz3,uuz4,aux3
     &     ,soa,bb,enrg,sig,u2334x,u2334y,u2334z,tsign
      REAL*8  coa,qforce,fact
      LOGICAL*4 near0,ok
      save nwind
      DATA nwind/0/

*==================== EXECUTABLE STATEMENTS ============================

      uconf=0.0D0
      fact=180.0D0/pi
      i=litor
      l1=litr(1,i)
      l2=litr(2,i)
      l3=litr(3,i)
      l4=litr(4,i)

*-----------------------------------------------------------------------
*------------ No bonded interaction exists between solvent and solute --
*-----------------------------------------------------------------------
              
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
      aux1=cb1*sb2
      aux2=sb1*cb2
      auxa=-(aux1*coa+aux2)/aux
      auxb=-(aux2*coa+aux1)/aux
      auxc=sb3/aux
      aux1=rsp21*rsp32*sb1
      aux2=rsp32*rsp43*sb2
      aux3=rsp21*rsp43*sb3
      dbx11= x32/aux1-cb1*x21/(rsq21*sb1)
      dbx13=-x21/aux1+cb1*x32/(rsq32*sb1)
      dbx12=-dbx13-dbx11
      dbx22= x43/aux2-cb2*x32/(rsq32*sb2)
      dbx24=-x32/aux2+cb2*x43/(rsq43*sb2)
      dbx23=-dbx22-dbx24
      dbx31=x43/aux3-cb3*x21/(rsq21*sb3)
      dbx33=x21/aux3-cb3*x43/(rsq43*sb3)
      dbx32=-dbx31
      dbx34=-dbx33
      dby11= y32/aux1-cb1*y21/(rsq21*sb1)
      dby13=-y21/aux1+cb1*y32/(rsq32*sb1)
      dby12=-dby13-dby11
      dby22= y43/aux2-cb2*y32/(rsq32*sb2)
      dby24=-y32/aux2+cb2*y43/(rsq43*sb2)
      dby23=-dby22-dby24
      dby31=y43/aux3-cb3*y21/(rsq21*sb3)
      dby33=y21/aux3-cb3*y43/(rsq43*sb3)
      dby32=-dby31
      dby34=-dby33
      dbz11= z32/aux1-cb1*z21/(rsq21*sb1)
      dbz13=-z21/aux1+cb1*z32/(rsq32*sb1)
      dbz12=-dbz13-dbz11
      dbz22= z43/aux2-cb2*z32/(rsq32*sb2)
      dbz24=-z32/aux2+cb2*z43/(rsq43*sb2)
      dbz23=-dbz22-dbz24
      dbz31=z43/aux3-cb3*z21/(rsq21*sb3)
      dbz33=z21/aux3-cb3*z43/(rsq43*sb3)
      dbz32=-dbz31
      dbz34=-dbz33
      bb=DACOS(coa)
      coa=bb
      u2334x=y32*z43-z32*y43
      u2334y=z32*x43-x32*z43
      u2334z=x32*y43-y32*x43
      aux=-(x21*u2334x+y21*u2334y+z21*u2334z)/rsp21/rsp32
     &     /rsp43/sb2
      tsign=DSIGN(1.0D0,aux)
      angle=-(tsign-1.0D0)*pi+tsign*coa+nwind*2.0D0*pi
      
      IF(angle-enout .LT. -pi) THEN
         nwind=nwind+1
         angle=-(tsign-1.0D0)*pi+tsign*coa+nwind*2.0D0*pi
      ELSE IF(angle-enout .GT. pi) THEN
         nwind=nwind-1
         angle=-(tsign-1.0D0)*pi+tsign*coa+nwind*2.0D0*pi
      END IF

      ok=.TRUE.
      IF(near0(enout)) enout=angle
      soa=DSIN(angle)
      
      IF(iflag .LT. 0) THEN
         IF(enout .GT. angle) THEN
            enout=angle
         ELSE
            ok=.FALSE.
         END IF
      ELSE
         IF(enout .LT. angle) THEN
            enout=angle
         ELSE
            ok=.FALSE.
         END IF
      END IF
      IF(.NOT. ok) THEN
         qforce=(2.0d0*alpha*(angle-enout))/soa
         uconf=alpha*(angle-enout)**2
         
         uux1=auxa*dbx11+auxc*dbx31
         uux2=auxa*dbx12+auxb*dbx22+auxc*dbx32
         uux3=auxa*dbx13+auxb*dbx23+auxc*dbx33
         uux4=auxb*dbx24+auxc*dbx34
         uuy1=auxa*dby11+auxc*dby31
         uuy2=auxa*dby12+auxb*dby22+auxc*dby32
         uuy3=auxa*dby13+auxb*dby23+auxc*dby33
         uuy4=auxb*dby24+auxc*dby34
         uuz1=auxa*dbz11+auxc*dbz31
         uuz2=auxa*dbz12+auxb*dbz22+auxc*dbz32
         uuz3=auxa*dbz13+auxb*dbz23+auxc*dbz33
         uuz4=auxb*dbz24+auxc*dbz34
         fpx(l1)=fpx(l1)+qforce*uux1
         fpx(l2)=fpx(l2)+qforce*uux2
         fpx(l3)=fpx(l3)+qforce*uux3
         fpx(l4)=fpx(l4)+qforce*uux4
         fpy(l1)=fpy(l1)+qforce*uuy1
         fpy(l2)=fpy(l2)+qforce*uuy2
         fpy(l3)=fpy(l3)+qforce*uuy3
         fpy(l4)=fpy(l4)+qforce*uuy4
         fpz(l1)=fpz(l1)+qforce*uuz1
         fpz(l2)=fpz(l2)+qforce*uuz2
         fpz(l3)=fpz(l3)+qforce*uuz3
         fpz(l4)=fpz(l4)+qforce*uuz4
      END IF

*================= END OF EXECUTABLE STATEMENTS ========================
         
      RETURN
      END
