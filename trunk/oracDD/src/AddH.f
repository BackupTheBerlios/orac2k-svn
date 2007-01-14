!!$/---------------------------------------------------------------------\
!!$                                                                      |
!!$  Copyright (C) 2006-2007 Massimo Marchi <Massimo.Marchi@cea.fr>      |
!!$                                                                      |
!!$      This program is free software;  you  can  redistribute  it      |
!!$      and/or modify it under the terms of the GNU General Public      |
!!$      License version 2 as published  by  the  Free  Software         |
!!$      Foundation;                                                     |
!!$                                                                      |
!!$      This program is distributed in the hope that  it  will  be      |
!!$      useful, but WITHOUT ANY WARRANTY; without even the implied      |
!!$      warranty of MERCHANTABILITY or FITNESS  FOR  A  PARTICULAR      |
!!$      PURPOSE.   See  the  GNU  General  Public License for more      |
!!$      details.                                                        |
!!$                                                                      |
!!$      You should have received a copy of the GNU General  Public      |
!!$      License along with this program; if not, write to the Free      |
!!$      Software Foundation, Inc., 59  Temple  Place,  Suite  330,      |
!!$      Boston, MA  02111-1307  USA                                     |
!!$                                                                      |
!!$\---------------------------------------------------------------------/
      SUBROUTINE addh(ia,nbl,nhb,blist,xh,yh,zh,type,x0,y0,z0,iret,
     x                errmsg)

************************************************************************
*                                                                      *
*     Compute the position of the hydrogen atoms bonded to a given     *
*     atom of the solute.                                              *
*                                                                      *
*     IA      :  Index of the central atom in the solute list.    (I)  *
*     NBL     :  Number of atoms bonded to IA other than hydrogens(I)  *
*     NHB     :  Number of hydrogen bonded to IA.                 (I)  *
*     BLIST   :  List of non-hydrogen atoms bonded to IA.         (I)  *
*                >> integer BLIST(*) <<                                *
*     XH      :  Cartesian coordinates of the solute atoms.       (I)  *
*     YH         >> real*8 XH(*), .... <<                              *
*     ZH                                                               *
*     TYPE    :  Type label of the atom IA.                       (I)  *
*     X0      :  Cartesian coordinates of the hydrogen atoms      (O)  *
*     Y0         >> real*8 X0(*), ... <<                               *
*     Z0                                                               *
*                                                                      *
*---- Last update 06/20/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNAL ADD1H, ADD2H, ADD3H                                     *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER ia,nbl,nhb,iret,blist(*)
      REAL*8  xh(*),yh(*),zh(*),x0(*),y0(*),z0(*)
      CHARACTER*7 type
      CHARACTER*80 errmsg

*-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8 xa(3),ya(3),za(3),ax(3),pax(3)
      REAL*8 length,oh,n2h,n3h,sh,c3h,c2h
      CHARACTER*1 oo,nn,ss,cc,char
      INTEGER i,j,ntot,n
      DATA oh,n2h,n3h,sh,c2h,c3h/0.960d0,1.01d0,1.01d0,1.336d0,
     x     1.08d0,1.09d0/
      DATA oo,nn,ss,cc/'o','n','s','c'/

*==================== EXECUTABLE STATEMENTS ============================

      ntot=nbl+nhb
      IF(ntot.GT.4) THEN
          iret=1
          errmsg=' IN ADDH: Cannot handle atoms with more than 4 bonds*
     x ABORT. '
          RETURN
      END IF

*=======================================================================
*---- Compute group axis -----------------------------------------------
*=======================================================================

      DO 10 i=1,nbl
          j=blist(i)
          xa(i)=xh(ia)-xh(j)
          ya(i)=yh(ia)-yh(j)
          za(i)=zh(ia)-zh(j)
!!$/---------------------------------------------------------------------\
!!$                                                                      |
!!$  Copyright (C) 2006-2007 Massimo Marchi <Massimo.Marchi@cea.fr>      |
!!$                                                                      |
!!$      This program is free software;  you  can  redistribute  it      |
!!$      and/or modify it under the terms of the GNU General Public      |
!!$      License version 2 as published  by  the  Free  Software         |
!!$      Foundation;                                                     |
!!$                                                                      |
!!$      This program is distributed in the hope that  it  will  be      |
!!$      useful, but WITHOUT ANY WARRANTY; without even the implied      |
!!$      warranty of MERCHANTABILITY or FITNESS  FOR  A  PARTICULAR      |
!!$      PURPOSE.   See  the  GNU  General  Public License for more      |
!!$      details.                                                        |
!!$                                                                      |
!!$      You should have received a copy of the GNU General  Public      |
!!$      License along with this program; if not, write to the Free      |
!!$      Software Foundation, Inc., 59  Temple  Place,  Suite  330,      |
!!$      Boston, MA  02111-1307  USA                                     |
!!$                                                                      |
!!$\---------------------------------------------------------------------/

10    CONTINUE
      DO 20 i=1,3
          ax(i)=0.0D0
20    CONTINUE
      DO 30 j=1,nbl
          ax(1)=ax(1)+xa(j)
          ax(2)=ax(2)+ya(j)
          ax(3)=ax(3)+za(j)
30    CONTINUE

*=======================================================================
*---- Calculate bond length from atom type -----------------------------
*=======================================================================

      char=type(1:1)
      IF(char.EQ.oo) THEN
          length=oh
      ELSE IF(char.EQ.nn) THEN
          IF(ntot.EQ.3) THEN
              length=n2h
          ELSE IF(ntot.EQ.4) THEN
              length=n3h
          END IF
      ELSE IF(char.EQ.ss) THEN
          length=sh
      ELSE IF(char.EQ.cc) THEN
          IF(ntot.EQ.3) THEN
              length=c2h
          ELSE IF(ntot.EQ.4) THEN
              length=c3h
          END IF
      END IF

*=======================================================================
*---- Compute position of the hydrogens w.r.t. the atom IA -------------
*=======================================================================

      IF(nhb.EQ.1) THEN
          CALL add1h(ax,length,x0,y0,z0)
      ELSE IF(nhb.EQ.3) THEN
          CALL add3h(ax,length,x0,y0,z0)
      ELSE IF(nhb.EQ.2) THEN
          IF(nbl.EQ.2) THEN
              pax(1)=ya(1)*za(2)-za(1)*ya(2)
              pax(2)=za(1)*xa(2)-xa(1)*za(2)
              pax(3)=xa(1)*ya(2)-ya(1)*xa(2)
              CALL add2h(ax,pax,1,length,x0,y0,z0)
          ELSE IF(nbl .EQ. 0) THEN
              CALL addh2o(x0,y0,z0)
          ELSE
              CALL add2h(ax,pax,0,length,x0,y0,z0)
          END IF
      END IF

*=======================================================================
*---- Compute cartesian coordinates of the hydrogens -------------------
*=======================================================================

      DO 40 i=1,nhb
          x0(i)=x0(i)+xh(ia)
          y0(i)=y0(i)+yh(ia)
          z0(i)=z0(i)+zh(ia)
40    CONTINUE

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
      SUBROUTINE add1h(ax,length,x0,y0,z0)


************************************************************************
*                                                                      *
*     ADD1H will produce coordinates for 1 hydrogen relative           *
*     to a central atom in (0.0, 0.0, 0.0).                            *
*                                                                      *
*     AX      :  Axis of the central atom. Need not be        (I)      *
*                normalized.                                           *
*                >> real*8 AX(3) <<                                    *
*     LENGTH  :  Length of the bond hydrogen-central atom.    (I)      *
*     X0      :  Hydrogen coordinates with respect to the     (O)      *
*     Y0         central atom.                                         *
*     Z0         >> real*8 X0, Y0, Z0 <<                               *
*                                                                      *
*---- Last update 07/08/92 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi CECAM, Orsay France 1992               *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      REAL*8 x0,y0,z0,ax(3),length

*-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8 norm,small
      INTEGER i
      DATA small/1.0D-4/

*==================== EXECUTABLE STATEMENTS ============================

      norm=0.0D0
      DO 10 i=1,3
          norm=norm+ax(i)**2
10    CONTINUE
      norm=DSQRT(norm)

      x0=ax(1)*length/norm+small
      y0=ax(2)*length/norm+small
      z0=ax(3)*length/norm+small

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
      SUBROUTINE add2h(ax,pax,iz,length,x0,y0,z0)

************************************************************************
*                                                                      *
*     ADD2H will produce coordinates for 2 hydrogens relative          *
*     to a central atom in (0.0, 0.0, 0.0).                            *
*                                                                      *
*     AX      :  Axis of the central atom. Need not be        (I)      *
*                normalized.                                           *
*                >> real*8 AX(3) <<                                    *
*     PAX     :  Axis perpendicular to the plane containing   (I)      *
*                the nearest neighbours of the central atom.           *
*                >> real*8 PAX(3) <<                                   *
*     IZ      :  It is 1 if PAX is defined and 0 otherwise.   (I)      *
*     LENGTH  :  Length of the bond hydrogen-central atom.    (I)      *
*     X0      :  Hydrogen coordinates with respect to the     (O)      *
*     Y0         central atom.                                         *
*     Z0         >> real*8 X0(2), Y0(2), Z0(2) <<                      *
*                                                                      *
*---- Last update 07/08/92 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi CECAM, Orsay France 1992.              *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      REAL*8 x0(*),y0(*),z0(*),ax(3),pax(3),length
      INTEGER iz

*-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8 x(2),y(2),z(2),ver(3),vrr(3),vsg(3),d0(3,3)
      REAL*8 norm,alpha,a,b,c,pi,sgns
      INTEGER i,j,m

*==================== EXECUTABLE STATEMENTS ============================

*=======================================================================
*---- Compute the normalized versors -----------------------------------
*=======================================================================

      norm=0.0D0
      DO i=1,3
          norm=norm+ax(i)**2
      END DO
      norm=DSQRT(norm)
      DO i=1,3
          ax(i)=ax(i)/norm
      END DO
      IF(iz.EQ.1) THEN
          norm=0.0D0
          DO i=1,3
              norm=norm+pax(i)**2
          END DO
          norm=DSQRT(norm)
          DO i=1,3
              pax(i)=pax(i)/norm
          END DO
      END IF

*=======================================================================
*---- Calculate the primitive coordinates of the hydrogens -------------
*---- and the versor perpendicular to their plane if required ----------
*=======================================================================

      pi=4.0D0*DATAN(1.0D0)
      alpha=0.50D0*(pi-DACOS(-1.0D0/3.0D0))
      IF(iz.NE.1) alpha=30.0D0*pi/180.0D0
      x(1)=0.0D0
      x(2)=0.0D0
      y(1)=DCOS(alpha)
      y(2)=-y(1)
      z(1)=DSIN(alpha)
      z(2)=z(1)
      a=ax(1)
      b=ax(2)
      c=ax(3)
      alpha=0.0D0
      CALL comd0(a,b,c,alpha,d0)
      IF(iz.EQ.1) THEN
          ver(1)=y(2)*z(1)-z(2)*y(1)
          ver(2)=z(2)*x(1)-x(2)*z(1)
          ver(3)=x(2)*y(1)-y(2)*x(1)
          DO i=1,3
              vrr(i)=0.0D0
              DO j=1,3
                  vrr(i)=vrr(i)+d0(i,j)*ver(j)
              END DO
          END DO
          norm=0.0D0
          DO i=1,3
              norm=norm+vrr(i)**2
          END DO
          norm=DSQRT(norm)
          DO i=1,3
              vrr(i)=vrr(i)/norm
          END DO
          alpha=DACOS(vrr(1)*pax(1)+vrr(2)*pax(2)+vrr(3)*pax(3))

*=======================================================================
*---- Compute the sign of the angle ALPHA ------------------------------
*=======================================================================

          vsg(1)=vrr(2)*pax(3)-pax(2)*vrr(3)
          vsg(2)=vrr(3)*pax(1)-pax(3)*vrr(1)
          vsg(3)=vrr(1)*pax(2)-pax(1)*vrr(2)
          norm=0.0D0
          DO i=1,3
              norm=norm+vsg(i)**2
          END DO       
          norm=DSQRT(norm)
          DO i=1,3
              vsg(i)=vsg(i)/norm
          END DO
          norm=pax(1)*ax(1)+pax(2)*ax(2)+pax(3)*ax(3)
          sgns=vsg(1)*ax(1)+vsg(2)*ax(2)+vsg(3)*ax(3)
          sgns=DSIGN(1.0D0,sgns)
          alpha=(0.50D0*pi-alpha)*sgns
          CALL comd0(a,b,c,alpha,d0)
      END IF

*=======================================================================
*---- Finally rotate the primitive coordinates and find the ------------
*---- required hydrogen coordinates ------------------------------------
*=======================================================================

      DO m=1,2
          x0(m)=d0(1,1)*x(m)+d0(1,2)*y(m)+d0(1,3)*z(m)
          y0(m)=d0(2,1)*x(m)+d0(2,2)*y(m)+d0(2,3)*z(m)
          z0(m)=d0(3,1)*x(m)+d0(3,2)*y(m)+d0(3,3)*z(m)
      END DO
      DO m=1,2
          x0(m)=x0(m)*length
          y0(m)=y0(m)*length
          z0(m)=z0(m)*length
      END DO

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
      SUBROUTINE add3h(ax,length,x0,y0,z0)

************************************************************************
*                                                                      *
*     ADD3H will produce coordinates for 3 hydrogens relative          *
*     to a central atom in (0.0, 0.0, 0.0).                            *
*                                                                      *
*     AX      :  Axis of the central atom. Need not be        (I)      *
*                normalized.                                           *
*                >> real*8 AX(3) <<                                    *
*     LENGTH  :  Length of the bond hydrogen-central atom.    (I)      *
*     X0      :  Hydrogen coordinates with respect to the     (O)      *
*     Y0         central atom.                                         *
*     Z0         >> real*8 X0(3), Y0(3), Z0(3) <<                      *
*                                                                      *
*---- Last update 06/20/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNAL COMD0                                                   *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      REAL*8 x0(*),y0(*),z0(*),ax(3),length

*-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8 x(3),y(3),z(3),d0(3,3)
      REAL*8 norm,alpha,beta,a,b,c,pi
      INTEGER i,m

*==================== EXECUTABLE STATEMENTS ============================

*=======================================================================
*---- Compute the normalized versors -----------------------------------
*=======================================================================

      norm=0.0D0
      DO 10 i=1,3
          norm=norm+ax(i)**2
10    CONTINUE
      norm=DSQRT(norm)
      DO 20 i=1,3
          ax(i)=ax(i)/norm
20    CONTINUE

*=======================================================================
*---- Calculate the primitive coordinates of the hydrogens -------------
*=======================================================================

      pi=4.0D0*DATAN(1.0D0)
      beta=60.0D0*pi/180.0D0
      alpha=DACOS(-1.0D0/3.0D0)-0.50D0*pi
      x(1)=0.0D0
      x(2)= - DCOS(alpha)*DSIN(beta)
      x(3)=DCOS(alpha)*DSIN(beta)
      y(1)=-DCOS(alpha)
      y(2)=DCOS(alpha)*DCOS(beta)
      y(3)=y(2)
      z(1)=DSIN(alpha)
      z(2)=z(1)
      z(3)=z(1)
      a=ax(1)
      b=ax(2)
      c=ax(3)
      alpha=0.0D0
      CALL comd0(a,b,c,alpha,d0)

*=======================================================================
*---- Finally rotate the primitive coordinates and find the ------------
*---- required hydrogen coordinates ------------------------------------
*=======================================================================

      DO 30 m=1,3
          x0(m)=d0(1,1)*x(m)+d0(1,2)*y(m)+d0(1,3)*z(m)
          y0(m)=d0(2,1)*x(m)+d0(2,2)*y(m)+d0(2,3)*z(m)
          z0(m)=d0(3,1)*x(m)+d0(3,2)*y(m)+d0(3,3)*z(m)
30    CONTINUE
      DO 40 m=1,3
          x0(m)=x0(m)*length
          y0(m)=y0(m)*length
          z0(m)=z0(m)*length
40    CONTINUE

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
      SUBROUTINE addh2o(x0,y0,z0)

************************************************************************
*                                                                      *
*     ADDH2O will produce the coordinates of the hydrogens for         *
*     a water molecule.                                                *
*                                                                      *
*     X0      :  Hydrogen coordinates with respect to the     (O)      *
*     Y0         central atom.                                         *
*     Z0         >> real*8 X0(2), Y0(2), Z0(2) <<                      *
*                                                                      *
*---- Last update 06/20/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNAL COMD0                                                   *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      REAL*8 x0(*),y0(*),z0(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8 x(2),y(2),z(2),vx,vy,vz,ux,uy,uz,ranf
      REAL*8 norm,pi,angle,length,alpha,eps,beta
      INTEGER m

!=======================================================================
*---- Parameters for SPC water -----------------------------------------
*=======================================================================

      DATA angle, length, eps/109.47D0, 1.0D0, 1.0D-8/

*==================== EXECUTABLE STATEMENTS ============================

*=======================================================================
*---- Compute a first randomly oriented versor -------------------------
*=======================================================================
      x(1)=2.0D0*(ranf()-0.5D0)
      y(1)=2.0D0*(ranf()-0.5D0)
      z(1)=2.0D0*(ranf()-0.5D0)
      norm=x(1)**2+y(1)**2+z(1)**2
      norm=DSQRT(norm)
      x(1)=x(1)/norm
      y(1)=y(1)/norm
      z(1)=z(1)/norm

*=======================================================================
*---- Compute the second versor at angle "angle" from the first --------
*=======================================================================

      pi=4.0D0*DATAN(1.0D0)
      alpha=angle*pi/180.0D0
      beta=(180-angle)*pi/180.0D0
      vx=2.0D0*(ranf()-0.5D0)
      vy=2.0D0*(ranf()-0.5D0)
      alpha=DCOS(alpha)
      vz=-vx*x(1)-vy*y(1)
      vz=vz/z(1)
      norm=vx**2+vy**2+vz**2
      norm=DSQRT(norm)
      vx=DSIN(beta)*vx/norm
      vy=DSIN(beta)*vy/norm
      vz=DSIN(beta)*vz/norm
      ux=-x(1)*DCOS(beta)
      uy=-y(1)*DCOS(beta)
      uz=-z(1)*DCOS(beta)
      x(2)=vx+ux
      y(2)=vy+uy
      z(2)=vz+uz
      norm=x(2)**2+y(2)**2+z(2)**2
      norm=DSQRT(norm)
      x(2)=x(2)/norm
      y(2)=y(2)/norm
      z(2)=z(2)/norm

*=======================================================================
*---- Create the two new vector ----------------------------------------
*=======================================================================

      DO m=1,2
          x0(m)=x(m)*length
          y0(m)=y(m)*length
          z0(m)=z(m)*length
      END DO
*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
      SUBROUTINE comd0(a,b,c,alpha,d0)

************************************************************************
*                                                                      *
*                                                                      *
*     COMD0 will calculate the rotation matrix D0 which rotate         *
*     a vector onto the direction CHI=(a,b,c) and apply a clockwise    *
*     rotation of alpha radians.                                       *
*                                                                      *
*     A       :  Components of the reference versor.              (I)  *
*     B                                                                *
*     C                                                                *
*     ALPHA   :  Clockwise rotation angle about (a,b,c).          (I)  *
*     D0      :  Rotation matrix.                                 (O)  *
*                                                                      *
*---- Last update 06/20/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNAL NONE                                                    *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT CHARACTER*80(a-z)

*----------------------- ARGUMENTS -------------------------------------

      REAL*8 a,b,c,alpha,d0(3,3)

*-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8 d,cs,sn

*==================== EXECUTABLE STATEMENTS ============================

      d=DSQRT(a**2+b**2)
      IF(DABS(d).GE.1.0D-8) THEN
          cs=DCOS(alpha)
          sn=DSIN(alpha)
          d0(1,1)=(b*cs-a*c*sn)/d
          d0(1,2)=(b*sn+a*c*cs)/d
          d0(1,3)=a
          d0(2,1)=-(a*cs+b*c*sn)/d
          d0(2,2)=(-a*sn+b*c*cs)/d
          d0(2,3)=b
          d0(3,1)=sn*d
          d0(3,2)=-cs*d
          d0(3,3)=c
      ELSE
          cs=DCOS(alpha)
          sn=DSIN(alpha)
          d0(1,1)=cs
          d0(1,2)=sn
          d0(1,3)=0.0D0
          d0(2,1)=-sn
          d0(2,2)=cs
          d0(2,3)=0.0D0
          d0(3,1)=0.0D0
          d0(3,2)=0.0D0
          d0(3,3)=1
      END IF

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
      FUNCTION ranf()
      REAL*8 ranf,duni
      ranf=duni()
      RETURN
      END
