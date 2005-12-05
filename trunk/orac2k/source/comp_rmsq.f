      SUBROUTINE comp_rmsq(beta,iatom,time,rms_disp,xpc,ypc,zpc
     &     ,tot_rms_disp,dsq,tnormb)

************************************************************************
*   Time-stamp: <01/01/18 16:55:46 sterpone>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Fri Nov 28 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER time,iatom
      REAL*8  xpc(*),ypc(*),zpc(*),rms_disp(*),tot_rms_disp(*)
     &     ,dsq(0:*),tnormb(*)
      CHARACTER*7 beta

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,no2,j,mwt
      REAL*8  sumsq,fcorr,drx,dry,drz,x(35000),y(35000),sig,siga,sigb
     &     ,chi2,a,b,q

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      mwt=0
      DO i=1,time
         tnormb(i) = 1.0D0/DFLOAT(time-i+1)
      END DO
      sumsq=0.0D0
      DO i=1,time
         dsq(i)=xpc(i)**2+ypc(i)**2+zpc(i)**2
         sumsq=sumsq+2.0D0*dsq(i)
      END DO
      dsq(0)=0.0D0
      dsq(time+1)=0.0D0
      DO i=1,time-1
         sumsq=sumsq-dsq(i-1)-dsq(time-i+2)
         tot_rms_disp(i) = tot_rms_disp(i)+tnormb(i)*(sumsq - 2.0D0
     &        *rms_disp(i))
         IF(beta .EQ. 'o1') THEN
            y(i)=tnormb(i)*(sumsq - 2.0D0*rms_disp(i))
            x(i)=i*0.24D0
         END IF
      END DO
      IF(beta .EQ. 'o1') THEN
         CALL fit(x(100),y(100),8333,sig,mwt,a,b,siga,sigb,chi2,q)
         WRITE(88,'(i8,3e15.7)') iatom,a,b,chi2
      END IF
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
