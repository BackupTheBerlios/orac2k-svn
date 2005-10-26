      SUBROUTINE fpbend(ss_index,lbndg,lbend,tag,sp,xp0,yp0,zp0,pota
     &     ,potb,potc,potd,ubend_slt,ubend_slv,fpx,fpy,fpz,ma)

************************************************************************
*                                                                      *
*     Compute the contribution from bond angles interaction            *
*     to macromolecular forces and energy.                             *
*                                                                      *
*     CO      :  Transform the positions to orthogo alized        (I)  *
*                crystallographic frame.                               *
*                >> real*8 CO(3,3) <<                                  *
*     LBNDG   :  List of bendings for the macromolecule.          (I)  *
*                >> integer LBNDG(3,N1) <<                             *
*     LBEND   :  Number of bendings.                              (I)  *
*     N1      :  Physical column dimension of LBEND.              (I)  *
*     XP0     :  Coordinates of the macromolecule.                (I)  *
*     YP0        >> real*8 XP0(NATO), YP0(NATO), ZP0(NATO) <<          *
*     ZP0                                                              *
*                                                                      *
*     NATO    :  Number of atoms which the macromolecule is       (I)  *
*                composed of.                                          *
*     POTA    :  List of potential parameters : bending           (I)  *
*                theta angles.                                         *
*                >> real*8 POTA(N1) <<                                 *
*     POTB    :  List of potential parameters : bending           (I)  *
*                force constant.                                       *
*                >> real*8 POTB(N1) <<                                 *
*     UBEND   :  Potential energy contribution from bendings.     (I)  *
*     FPX     :  List of forces acting on the macromolecule       (I)  *
*     FPY        atoms.                                                *
*     FPZ        >> real*8 FPX(NATO), FPY(NATO), FPZ(NATO) <<          *
*                                                                      *
*---- Last update 04/24/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNALS NONE                                                   *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER ma,lbend,ss_index(*),tag(*),sp(*)
      INTEGER lbndg(3,*)
      REAL*8  ubend_slt,ubend_slv
      REAL*8  xp0(*),yp0(*),zp0(*),pota(*),potb(*),potc(*),potd(*)
     &     ,fpx(ma,*),fpy(ma,*),fpz(ma,*)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,la,lb,lc,type,tags,mma,iam,ii
      REAL*8  xr1,xr2,xr3,yr1,yr2,yr3,zr1,zr2,zr3,x12,x32,y12,y32,z12
     &     ,z32,rs12,rs32,uux1,uux2,uux3,uuy1,uuy2,uuy3,uuz1,uuz2,uuz3
     &     ,xr31,yr31,zr31,rs31,rsp31,x31,y31,z31,rsp12,rsp32,k12,r1,r2
      REAL*8  dcc2,cb,sb,bb,qforce,pforce,pi,ubend(2),fact,factor
     &     ,quarter,half,quasi_zero
      logical coupling 
      DATA quasi_zero/1.0D-12/
      DATA quarter,half/0.25D0,0.5D0/
      FACTOR(i)=(1.0D0+DBLE(i/IABS(i)))*half

*==================== EXECUTABLE STATEMENTS ============================

      DO i=1,2
         ubend(i)=0.0D0
      END DO
      pi=4.0D0*DATAN(1.0D0)

*----- bonded interactions: bending

      mma=sp(1)
      DO ii=1,mma
         iam=sp(1+ii)
         i=IABS(iam)
         fact=FACTOR(iam)
         tags=tag(i)
         la=lbndg(1,i)
         lb=lbndg(2,i)
         lc=lbndg(3,i)

*-----------------------------------------------------------------------
*------------ No bonded interaction exists between solvent and solute --
*-----------------------------------------------------------------------
              
         type=ss_index(la)
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
         sb=DSQRT(1.0d0-cb**2)
         bb=DACOS(cb)
c--      switch to singularity free potential for 
c--      linear bending V = 2*K*(cos(theta) + 1)  
         if(abs(pota(i)-pi).lt.0.01) THEN 
           qforce=-2.0d0*potb(i)
         else
           qforce=2.0d0*potb(i)*(bb-pota(i))/sb
         end if
         uux1=x32/dcc2-cb*x12/rs12
         uux3=x12/dcc2-cb*x32/rs32
         uux2=-uux1-uux3
         uuy1=y32/dcc2-cb*y12/rs12
         uuy3=y12/dcc2-cb*y32/rs32
         uuy2=-uuy1-uuy3
         uuz1=z32/dcc2-cb*z12/rs12
         uuz3=z12/dcc2-cb*z32/rs32
         uuz2=-uuz1-uuz3
         fpx(la,tags)=fpx(la,tags)+qforce*uux1
         fpx(lc,tags)=fpx(lc,tags)+qforce*uux3
         fpx(lb,tags)=fpx(lb,tags)+qforce*uux2
         fpy(la,tags)=fpy(la,tags)+qforce*uuy1
         fpy(lc,tags)=fpy(lc,tags)+qforce*uuy3
         fpy(lb,tags)=fpy(lb,tags)+qforce*uuy2
         fpz(la,tags)=fpz(la,tags)+qforce*uuz1
         fpz(lc,tags)=fpz(lc,tags)+qforce*uuz3
         fpz(lb,tags)=fpz(lb,tags)+qforce*uuz2
c--      switch to singularity free potential for 
c--      linear bending V = 2*K(cos(theta) + 1)  
         if(abs(pota(i)-pi).lt.0.01) THEN 
           ubend(type)=ubend(type)+fact*2.d0*potb(i)*(cb + 1.d0) 
         else
           ubend(type)=ubend(type)+fact*potb(i)*(bb-pota(i))**2
         end if

*=======================================================================
*---  Compute the Urey-Bradley term ------------------------------------
*=======================================================================

         x31=xr3-xr1
         y31=yr3-yr1
         z31=zr3-zr1
         rs31=x31**2+y31**2+z31**2
         rsp31=DSQRT(rs31)
         
         qforce=-2.0D0*potd(i)*(potc(i)-rsp31)
         uux1=x31/rsp31
         uuy1=y31/rsp31
         uuz1=z31/rsp31
         uux2=-uux1
         uuy2=-uuy1
         uuz2=-uuz1
         fpx(la,tags)=fpx(la,tags)+qforce*uux1
         fpy(la,tags)=fpy(la,tags)+qforce*uuy1
         fpz(la,tags)=fpz(la,tags)+qforce*uuz1
         fpx(lc,tags)=fpx(lc,tags)+qforce*uux2
         fpy(lc,tags)=fpy(lc,tags)+qforce*uuy2
         fpz(lc,tags)=fpz(lc,tags)+qforce*uuz2
         ubend(type)=ubend(type)+fact*potd(i)*(rsp31-potc(i))**2 
      END DO

      ubend_slt=ubend(1)
      ubend_slv=ubend(2)

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
