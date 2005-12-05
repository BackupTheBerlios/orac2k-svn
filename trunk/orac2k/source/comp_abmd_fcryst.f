      SUBROUTINE comp_abmd_fcryst(iflag,skkset,nkvect,kvect,alpha,ntap
     &     ,xp0,yp0,zp0,co,fppx,fppy,fppz,uconf,skk)

************************************************************************
*   Time-stamp: <99/04/07 15:05:11 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Jun 13 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER iflag,nkvect,ntap
      REAL*8  skkset,skk,kvect(3,*),alpha,uconf,co(3,3)
      REAL*8  xp0(*),yp0(*),zp0(*),fppx(*),fppy(*),fppz(*)

*---------------------- VARIABLES IN COMMONS --------------------------*
      INCLUDE 'parst.h'
      REAL*8  fpx(m1),fpy(m1),fpz(m1),fskt(10),fski(10),fskr(10)
      COMMON /rag1/ fpx,fpy,fpz

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,m
      REAL*8  kx,ky,kz,fsk,scal,qx,qy,qz,coef,mult,xc,yc,zc
      LOGICAL near0
      
*----------------------- EXECUTABLE STATEMENTS ------------------------*

      uconf=0.0D0
      CALL zeroa(fpx,fpy,fpz,ntap,1)
      fsk=0.0D0
      DO m=1,nkvect
         kx=kvect(1,m)
         ky=kvect(2,m)
         kz=kvect(3,m)
         fskr(m)=0.0D0
         fski(m)=0.0D0
         DO i=1,ntap
            xc=co(1,1)*xp0(i)+co(1,2)*yp0(i)+co(1,3)*zp0(i)
            yc=co(2,1)*xp0(i)+co(2,2)*yp0(i)+co(2,3)*zp0(i)
            zc=co(3,1)*xp0(i)+co(3,2)*yp0(i)+co(3,3)*zp0(i)
            qx=kx*xc
            qy=ky*yc
            qz=kz*zc
            scal=qx+qy+qz
            fskr(m)=fskr(m)+DCOS(scal)
            fski(m)=fski(m)+DSIN(scal)
         END DO
         fskr(m)=fskr(m)/DFLOAT(ntap)
         fski(m)=fski(m)/DFLOAT(ntap)
         fskt(m)=fskr(m)**2+fski(m)**2
         fsk=fsk+fskt(m)
      END DO
      fsk=fsk/DFLOAT(nkvect)

      IF(near0(skkset)) skkset=fsk
      WRITE(6,*) skkset,fsk
      IF(iflag .LT. 0) THEN
         IF(skkset .GE. fsk) THEN
            skkset=fsk
         ELSE
            coef=4.0D0*alpha*(fsk-skkset)/DFLOAT(ntap*nkvect)
            uconf=alpha*(fsk-skkset)**2
            DO m=1,nkvect
               kx=kvect(1,m)
               ky=kvect(2,m)
               kz=kvect(3,m)
               DO i=1,ntap
                  xc=co(1,1)*xp0(i)+co(1,2)*yp0(i)+co(1,3)*zp0(i)
                  yc=co(2,1)*xp0(i)+co(2,2)*yp0(i)+co(2,3)*zp0(i)
                  zc=co(3,1)*xp0(i)+co(3,2)*yp0(i)+co(3,3)*zp0(i)
                  qx=kx*xc
                  qy=ky*yc
                  qz=kz*zc
                  scal=qx+qy+qz
                  mult=coef*(DSIN(scal)*fskr(m)-DCOS(scal)*fski(m))
                  fpx(i)=fpx(i)+mult*kx
                  fpy(i)=fpy(i)+mult*ky
                  fpz(i)=fpz(i)+mult*kz
               END DO
            END DO
         END IF
      ELSE
         IF(skkset .LE. fsk) THEN
            skkset=fsk
         ELSE
            coef=4.0D0*alpha*(fsk-skkset)/DFLOAT(ntap*nkvect)
            uconf=alpha*(fsk-skkset)**2
            DO m=1,nkvect
               kx=kvect(1,m)
               ky=kvect(2,m)
               kz=kvect(3,m)
               DO i=1,ntap
                  xc=co(1,1)*xp0(i)+co(1,2)*yp0(i)+co(1,3)*zp0(i)
                  yc=co(2,1)*xp0(i)+co(2,2)*yp0(i)+co(2,3)*zp0(i)
                  zc=co(3,1)*xp0(i)+co(3,2)*yp0(i)+co(3,3)*zp0(i)
                  qx=kx*xc
                  qy=ky*yc
                  qz=kz*zc
                  scal=qx+qy+qz
                  mult=coef*(DSIN(scal)*fskr(m)-DCOS(scal)*fski(m))
                  fpx(i)=fpx(i)+mult*kx
                  fpy(i)=fpy(i)+mult*ky
                  fpz(i)=fpz(i)+mult*kz
               END DO
            END DO
         END IF
      END IF

      DO i=1,ntap
         fppx(i)=fppx(i)+co(1,1)*fpx(i)+co(1,2)*fpy(i)+co(1,3)*fpz(i)
         fppy(i)=fppy(i)+co(2,1)*fpx(i)+co(2,2)*fpy(i)+co(2,3)*fpz(i)
         fppz(i)=fppz(i)+co(3,1)*fpx(i)+co(3,2)*fpy(i)+co(3,3)*fpz(i)
      END DO

      skk=fsk
      WRITE(6,*) uconf

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
