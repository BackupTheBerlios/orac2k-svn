      SUBROUTINE plot_fragm(fstep,beta,nfragm,fragm,xp0,yp0,zp0,ntap
     &     ,nres,m1,prsymb,charge)

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

      INTEGER ntap,m1,nstart,nend,nfragm
      INTEGER nres(m1,*),fragm(2,*)
      REAL*8  xp0(*),yp0(*),zp0(*),charge(*),fstep
      CHARACTER*7 beta(*)
      CHARACTER*8 prsymb(*)

*------------------- VARIABLES IN COMMON -------------------------------

      INCLUDE 'unit.h'

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,k,l,n
      REAL*8  xb,yb,zb,sunitc
      CHARACTER*5 bet,bet2
      CHARACTER*3 rsd

*==================== EXECUTABLE STATEMENTS ============================

      
      WRITE(kplot_fragm,2) fstep
      sunitc=DSQRT(unitc)
      k=0
      DO n=1,nfragm
         DO i=fragm(1,n),fragm(2,n)
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
            WRITE(kplot_fragm,1)'ATOM  ',i,bet2(1:5),rsd,nres(i,1),xb,yb
     &           ,zb,charge(i)*sunitc,DFLOAT(k)
         END DO
      END DO
      WRITE(kplot_fragm,'(a)')'TER  '
1     FORMAT(a6,i5,1x,a5,a3,2x,i4,4x,3f8.3,2f6.2)
2     FORMAT('COMMENT  1 Configuration at time step ',f11.2,'      ',
     &'                 ')
3     FORMAT('COMMENT  2 ABMD varariables GR ',f12.4,' GRA ',f12.4)

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
