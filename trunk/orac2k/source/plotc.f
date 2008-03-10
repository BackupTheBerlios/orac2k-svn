      SUBROUTINE plotc(co,abmd,gr,gra,fstep,beta,xp0,yp0,zp0,ntap,nres
     &     ,m1,prsymb,charge)

************************************************************************
*                                                                      *
*     Print a snapshot of the solute molecules on unit kplot.          *
*                                                                      *
*     CO      :  Transform the coordinates to the orthogonalized  (I)  *
*                crystallographic frame.                               *
*                >> real*8  CO(3,3) <<                                 *
*     BETA    :  List of labels for each site of the solute.      (I)  *
*                >> character*7 BETA(NTAP) <<                          *
*     XP0     :  Coordinates of the solute molecules, packed by   (I)  *
*     YP0        site.                                                 *
*     ZP0        >> real*8  XP0(NTAP), ... <<                          *
*                                                                      *
*     NTAP    :  Number of sites on the solute molecules.         (I)  *
*                                                                      *
*---- Last update 05/03/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNALS NONE.                                                  *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER ntap,m1
      INTEGER nres(m1,*)
      REAL*8  xp0(*),yp0(*),zp0(*),charge(*),fstep,gr,gra,co(3,3)
      LOGICAL abmd
      CHARACTER*7 beta(*)
      CHARACTER*8 prsymb(*)

*------------------- VARIABLES IN COMMON -------------------------------

      INCLUDE 'unit.h'

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,k,l
      REAL*8  xb,yb,zb,sunitc,a,b,c,alf,bett,gamm
      CHARACTER*5 bet,bet2
      CHARACTER*3 rsd

*==================== EXECUTABLE STATEMENTS ============================


      CALL rotb(a,b,c,alf,bett,gamm,co)
      WRITE(kout,2) fstep
# Les lignes suivantes sont commentees expres pour faire le film
c      WRITE(kout,4) a,b,c,alf,bett,gamm
      WRITE(kout,5) a,b,c,alf,bett,gamm
c      IF(abmd) WRITE(kout,3) gr,gra
      sunitc=DSQRT(unitc)
      k=0
      DO i=1,ntap
         j=nres(i,2)
         xb=xp0(i)
         yb=yp0(i)
         zb=zp0(i)
         bet(1:5)=beta(i)(1:5)
         CALL low_up(bet,5)
         DO l=1,5
            bet2(l:l)=bet(l:l)
         END DO
         rsd(1:3)=prsymb(j)(1:3)
         CALL low_up(rsd,3)
         WRITE(kout,1)'ATOM  ',i,bet2(1:5),rsd,nres(i,1),xb,yb,zb
     &        ,charge(i)*sunitc,DFLOAT(k)
      END DO
      WRITE(kout,'(a)')'END  '
1     FORMAT(a5,i6,1x,a5,a3,1x,i5,4x,3f8.3,2f6.2)
2     FORMAT('COMMENT 1 Configuration at time step ',f11.2,'      ',
     &'                 ')
3     FORMAT('COMMENT 2 ABMD varariables GR ',f12.4,' GRA ',f12.4)
4     FORMAT('COMMENT 1 Cell parameters ',3f10.5/
     &       'COMMENT 1 Cell parameters ',3f10.5)
5     FORMAT('CRYST1',3f9.3,3f7.2,' P  1          1')
C 12345678901234567890123456789012345678901234567890123456789012345678901234567890
C CRYST1   58.397   58.397   42.091  90.00  90.00  90.00 P 43          4  1HRC  83

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
