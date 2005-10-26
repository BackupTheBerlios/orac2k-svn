      SUBROUTINE virtual_extract(nbun,atres,beta,xp0,yp0,zp0,nres,m1
     &     ,prsymb,charge,virtual_atoms)

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

      INTEGER ntap,m1,nbun,atres(2,*)
      INTEGER nres(m1,*),virtual_atoms(*),uu
      REAL*8  xp0(*),yp0(*),zp0(*),charge(*)
      CHARACTER*7 beta(*)
      CHARACTER*8 prsymb(*)

*------------------- VARIABLES IN COMMON -------------------------------

      INCLUDE 'unit.h'

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,k,l,ii,kk,ll,m
      REAL*8  xb,yb,zb,sunitc,xc,yc,zc,rsp
      CHARACTER*5 bet,bet2
      CHARACTER*3 rsd
      LOGICAL flag

*==================== EXECUTABLE STATEMENTS ============================

      
      sunitc=DSQRT(unitc)
      m=virtual_atoms(1)
      k=0
      DO ii=1,nbun
         flag=.TRUE.
         DO i=atres(1,ii),atres(2,ii)
            xb=xp0(i)
            yb=yp0(i)
            zb=zp0(i)
            DO ll=1,m
               l=virtual_atoms(1+ll)
               xc=xp0(l)-xb
               yc=yp0(l)-yb
               zc=zp0(l)-zb
               rsp=DSQRT(xc**2+yc**2+zc**2)
               IF(ll .GT. 6 .AND. ll .LT. 15) THEN
                  IF(rsp .LT. 18.0D0) flag=.FALSE.
               END IF
                     
            END DO
         END DO
         IF(flag) THEN
            uu=uu+1
            DO i=atres(1,ii),atres(2,ii)
               kk=kk+1
               j=nres(i,2)
               xb=xp0(i)
               yb=yp0(i)
               zb=zp0(i)
               flag=.TRUE.
               bet(1:5)=beta(i)(1:5)
               CALL low_up(bet,5)
               DO l=1,5
                  bet2(l:l)=bet(l:l)
               END DO
               rsd(1:3)=prsymb(j)(1:3)
               CALL low_up(rsd,3)
               WRITE(88,1)'ATOM  ',kk,bet2(1:5),rsd,uu,xb,yb,zb
     &              ,charge(i)*sunitc,DBLE(k)
            END DO
         END IF
      END DO
      WRITE(88,'(a)')'TER  '
1     FORMAT(a6,i5,1x,a5,a3,1x,i5,4x,3f8.3,2f6.2)
2     FORMAT('COMMENT  1 Configuration at time step ',f11.2,'      ',
     &'                 ')
3     FORMAT('COMMENT  2 ABMD varariables GR ',f12.4,' GRA ',f12.4)

*================= END OF EXECUTABLE STATEMENTS ========================

      STOP
      END
