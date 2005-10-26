      SUBROUTINE fpbond(ss_index,lbnd,lbond,sp,xp0,yp0,zp0,pota
     &     ,potb,ubond_slt,ubond_slv,fpx,fpy,fpz)

************************************************************************
*                                                                      *
*     Compute the contribution from bond stretching to solute          *
*     forces and energy.                                               *
*                                                                      *
*     CO      :  Transform the positions to orthogo alized        (I)  *
*                crystallographic frame.                               *
*                >> real*8 CO(3,3) <<                                  *
*     LBND    :  List of bends for the solute.                    (I)  *
*                >> integer LBNDG(3,N1) <<                             *
*     LBOND   :  Number of bonds.                                 (I)  *
*     N1      :  Physical column dimension of LBEND.              (I)  *
*     XP0     :  Coordinates of the macromolecule.                (I)  *
*     YP0        >> real*8 XP0(NATO), YP0(NATO), ZP0(NATO) <<          *
*     ZP0                                                              *
*                                                                      *
*     NATO    :  Number of atoms which the macromolecule is       (I)  *
*                composed of.                                          *
*     POTA    :  List of potential parameters : bond              (I)  *
*                dostance.                                             *
*                >> real*8 POTA(N1) <<                                 *
*     POTB    :  List of potential parameters : bond              (I)  *
*                force constant.                                       *
*                >> real*8 POTB(N1) <<                                 *
*     UBOND   :  Potential energy contribution from bond.         (I)  *
*     FPX     :  List of forces acting on the macromolecule       (I)  *
*     FPY        atoms.                                                *
*     FPZ        >> real*8 FPX(NATO), FPY(NATO), FPZ(NATO) <<          *
*                                                                      *
*---- Last update 05/17/94 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi CEA, CE Saclay FRANCE 1994             *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER lbond,ss_index(*),sp(*)
      INTEGER lbnd(2,*)
      REAL*8  ubond_slt,ubond_slv
      REAL*8  xp0(*),yp0(*),zp0(*),pota(*),potb(*),
     x        fpx(*),fpy(*),fpz(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,la,lb,type
      REAL*8  xr1,xr2,yr1,yr2,zr1,zr2,x21,y21,z21,rs21,
     x        uux1,uux2,uuy1,uuy2,uuz1,uuz2,ubond(2)
      REAL*8  qforce
      INTEGER ii,mma,iam
      REAL*8 fact,quarter,half,factor
      DATA quarter,half/0.25D0,0.5D0/
c      FACTOR(i)=(1.0D0+DBLE(i/IABS(i)))*quarter+half
      FACTOR(i)=(1.0D0+DBLE(i/IABS(i)))*half

*==================== EXECUTABLE STATEMENTS ============================
      
      DO i=1,2
         ubond(i)=0.0D0
      END DO

*----- bonded interactions: bond stretching

      mma=sp(1)
      DO ii=1,mma
         iam=sp(1+ii)
         i=IABS(iam)
         fact=FACTOR(iam)
         la=lbnd(1,i)
         lb=lbnd(2,i)
         
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
         
         x21=xr2-xr1
         y21=yr2-yr1
         z21=zr2-zr1
         rs21=DSQRT(x21**2+y21**2+z21**2)
         
         qforce=-2.0D0*potb(i)*(pota(i)-rs21)
         uux1=x21/rs21
         uuy1=y21/rs21
         uuz1=z21/rs21
         uux2=-uux1
         uuy2=-uuy1
         uuz2=-uuz1
         fpx(la)=fpx(la)+qforce*uux1
         fpy(la)=fpy(la)+qforce*uuy1
         fpz(la)=fpz(la)+qforce*uuz1
         
         fpx(lb)=fpx(lb)+qforce*uux2
         fpy(lb)=fpy(lb)+qforce*uuy2
         fpz(lb)=fpz(lb)+qforce*uuz2
         ubond(type)=ubond(type)+fact*potb(i)*(rs21-pota(i))**2
      END DO
      
      ubond_slt=ubond(1)
      ubond_slv=ubond(2)
      
*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
