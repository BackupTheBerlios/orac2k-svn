      SUBROUTINE timcor(ss_index,c,co,vpx,vpy,vpz,nato,dt,nstep,ndelta
     &     ,nprint,nintv)
************************************************************************
*   Time-stamp: <97/06/16 17:38:30 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Feb 20 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8 vpx(*),vpy(*),vpz(*),dt,co(3,3)
      INTEGER nato,nstep,ndelta,nprint,nintv,ss_index(*)
      LOGICAL c

*----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'unit.h'
      INCLUDE 'vacf.h'

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER nbin,n,m,i,nlast,ni,nj,indx,type
      REAL*8  scal,vacft,norm,wf

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      nbin=nstep/ndelta+1
      IF(c) THEN
         nbin=nstep/ndelta
      END IF
      IF(nbin .LE. nintv) THEN
         indx=nbin
         DO m=1,nato
            vpl(1,m,indx)=vpx(m)
            vpl(2,m,indx)=vpy(m)
            vpl(3,m,indx)=vpz(m)
         END DO
      ELSE
         indx=nintv
         DO n=2,nintv
            DO m=1,nato
               DO i=1,3
                  vpl(i,m,n-1)=vpl(i,m,n)
               END DO
            END DO
         END DO
         DO m=1,nato
            vpl(1,m,nintv)=vpx(m)
            vpl(2,m,nintv)=vpy(m)
            vpl(3,m,nintv)=vpz(m)
         END DO
      END IF

      IF(nbin .LT. nintv) THEN
         nlast=nbin
      ELSE
         nlast=nintv
      END IF

      DO ni=1,nlast
         nj=nlast
         i=nj-ni
         DO m=1,nato
            type=ss_index(m)
            scal=vpl(1,m,ni)*vpl(1,m,nj)+vpl(2,m,ni)*vpl(2,m,nj)+
     x           vpl(3,m,ni)*vpl(3,m,nj)
            vafp(i,type)=vafp(i,type)+scal
            weip(i,type)=weip(i,type)+1.0
         END DO
      END DO
      IF(mod(nbin,nprint) .EQ. 0) THEN
         REWIND kvaf
         norm=vafp(0,1)/weip(0,1)
         DO i=0,nlast-1
            vacft=vafp(i,1)/weip(i,1)
            vacft=vacft/norm
            vacf_slt(i)=vacft
         END DO
         norm=vafp(0,2)/weip(0,2)
         DO i=0,nlast-1
            vacft=vafp(i,2)/weip(i,2)
            vacft=vacft/norm
            vacf_slv(i)=vacft
         END DO

         call ffct4(vacf_slt,psp_slt,nlast-1)
         call ffct4(vacf_slv,psp_slt,nlast-1)
         norm=(vafp(0,1)/weip(0,1))+(vafp(0,2)/weip(0,2))
         DO i=0,nlast-1
            vacft=(vafp(i,1)/weip(i,1))+(vafp(i,2)/weip(i,2))
            vacf_tot(i)=vacft/norm
         END DO
         CALL ffct4(vacf_tot,psp_tot,nlast-1)

            
         wf=1.0D5/(2.0D0*2.997925D0*dt)/DBLE(nlast-1)
         WRITE(kvaf,100) 
         DO i=0,nlast-1
            WRITE(kvaf,200) DBLE(i)*dt,vacf_tot(i),vacf_slv(i)
     &           ,vacf_slt(i),DBLE(i)*wf,psp_tot(i),psp_slv(i)
     &           ,psp_slt(i)
         END DO

      END IF

100   FORMAT(/
     &'                 =======           VACF and power spectrum     ',
     &'======='//
     &'                       VACF                                    ',
     &'Intensity'/
     &'     t(fs)    tot       slv       slt           w(cm-1)     tot',
     &'           slv          slt'/)
200   FORMAT(f10.3,3f10.6,5x,f10.3,3e13.5)
      RETURN
      END
