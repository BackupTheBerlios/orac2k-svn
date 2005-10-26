      SUBROUTINE write_rms(krms,xp_avg,yp_avg,zp_avg,xp_avg2,yp_avg2
     &     ,zp_avg2,nstep,fstep,ngrp,grppt,res,beta)

************************************************************************
*   Time-stamp: <97/07/15 19:52:47 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Jul 15 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER ngrp,grppt(2,*),res(*),krms,nstep
      REAL*8  xp_avg(*),yp_avg(*),zp_avg(*),xp_avg2(*),yp_avg2(*)
     &     ,zp_avg2(*),fstep
      CHARACTER*7 beta(*)

*----------------------- VARIABLES IN COMMON --------------------------*


*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER j,i
      REAL*8  fnstep,FLUCT,avg2,avg,step,xd,yd,zd,xc,yc,zc,xf,yf,zf,rsd
     &     ,rsc,rf
      FLUCT(avg2,step,avg)=DABS((avg2-step*avg*avg)/step)

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      WRITE(krms,'('' Tstep = '',f12.2,2x,a2)') fstep
      WRITE(krms,'(''                         SQRT(FLUCT(R))  FLUCT(X)  
     &  FLUCT(Y)    FLUCT(Z)'')') 
      fnstep=DBLE(nstep)
      DO j=1,ngrp
         DO i=grppt(1,j),grppt(2,j)
            xc=xp_avg(i)/fnstep
            yc=yp_avg(i)/fnstep
            zc=zp_avg(i)/fnstep
            xf=FLUCT(xp_avg2(i),fnstep,xc)
            yf=FLUCT(yp_avg2(i),fnstep,yc)
            zf=FLUCT(zp_avg2(i),fnstep,zc)
            rf=xf+yf+zf
            rf=dsqrt(rf)
            IF(DABS(rf) .GT. 1.0D-4) THEN
               WRITE(krms,1000) beta(i)(1:1),i,j,res(i),rf,xf,yf,zf
            END IF
         END DO
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

1000  FORMAT(1x,1a,2x,i6,1x,i6,1x,i5,2x,f12.6,f12.6,f12.6,f12.6)

      RETURN
      END
