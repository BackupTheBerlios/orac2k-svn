      SUBROUTINE adjust_bonds(ss_index,ntap,lcnstr,lconstr,xp0,yp0
     &     ,zp0,potbo,ma,iret,errmsg)

************************************************************************
*   Time-stamp: <2005-03-13 23:33:47 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Feb 23 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER ntap,lconstr,lcnstr(2,*),ss_index(*),iret,ma
      REAL*8  xp0(*),yp0(*),zp0(*),potbo(ma,*)
      CHARACTER*80 errmsg

*----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'parst.h'
      REAL*8  fpx(m1),fpy(m1),fpz(m1),fpx1(m1),fpy1(m1),fpz1(m1),xp1(m1)
     &     ,yp1(m1),zp1(m1),gpx(m1),gpy(m1),gpz(m1),hpx(m1),hpy(m1)
     &     ,hpz(m1)
      INTEGER sp(m1)
      COMMON /rag1/ gpx,gpy,gpz,hpx,hpy,hpz,fpx,fpy,fpz,fpx1,fpy1,fpz1
     &     ,xp1,yp1,zp1,sp

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,la,lb,count
      LOGICAL ok
      REAL*8  tol,ubond_slt,ubond_slv,ubond,fret,gg,dgg,gam,xab
     &     ,yab,zab,dpp,tol2
      DATA tol/1.0D-6/

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      sp(1)=lconstr
      DO i=1,lconstr
         sp(1+i)=i
      END DO
      tol2=tol*tol
      DO i=1,ntap
         fpx(i)=0.0D0
         fpy(i)=0.0D0
         fpz(i)=0.0D0
      END DO
      ubond_slt=0.0D0
      ubond_slv=0.0D0
      CALL fpbond(ss_index,lcnstr,lconstr,sp,xp0,yp0,zp0,potbo(1,2)
     &     ,potbo(1,1),ubond_slt,ubond_slv,fpx,fpy,fpz)
      ubond=ubond_slt+ubond_slv

      DO i=1,ntap
         gpx(i)=fpx(i)
         gpy(i)=fpy(i)
         gpz(i)=fpz(i)
         hpx(i)=gpx(i)
         hpy(i)=gpy(i)
         hpz(i)=gpz(i)
      END DO

      ok=.FALSE.
      count=0
      DO WHILE(.NOT. ok)
         IF(count .GT. 200) THEN
            iret=1
            errmsg=
     &' Bond length could not be adjusted to the eq. distance'
     &//' within allowed iterations.'
            RETURN
         END IF
         count=count+1
         fret=0.0D0
         CALL linmin_adjust_bonds(xp0,yp0,zp0,xp1,yp1,zp1,fpx,fpy,fpz
     &     ,fpx1,fpy1,fpz1,ntap,fret)
         DO i=1,ntap
            fpx(i)=0.0D0
            fpy(i)=0.0D0
            fpz(i)=0.0D0
         END DO
         CALL fpbond(ss_index,lcnstr,lconstr,sp,xp0,yp0,zp0,potbo(1,2)
     &        ,potbo(1,1),ubond_slt,ubond_slv,fpx,fpy,fpz)
         ubond=(ubond_slt+ubond_slv)/DBLE(lconstr)
         IF(DABS(ubond) .LE. tol2) ok=.TRUE.
         gg=0.0D0
         dgg=0.0D0
         DO i=1,ntap
            gg=gg+gpx(i)**2+gpy(i)**2+gpz(i)**2
            dgg=dgg+fpx(i)**2+fpy(i)**2+fpz(i)**2
         END DO
         IF(DABS(gg) .EQ. 0.0D0) RETURN
         gam=dgg/gg
         DO i=1,ntap
            gpx(i)=fpx(i)
            gpy(i)=fpy(i)
            gpz(i)=fpz(i)
            hpx(i)=gpx(i)+gam*hpx(i)
            hpy(i)=gpy(i)+gam*hpy(i)
            hpz(i)=gpz(i)+gam*hpz(i)
            fpx(i)=hpx(i)
            fpy(i)=hpy(i)
            fpz(i)=hpz(i)
         END DO
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
