      SUBROUTINE comp_native(knative,co,xp0,yp0,zp0,beta,nbun,atres
     &     ,native_dist,res)

************************************************************************
*   Time-stamp: <01/08/14 19:09:31 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Fri Mar 12 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER knative,nbun,res(*),atres(2,*)
      REAL*8  xp0(*),yp0(*),zp0(*),co(3,3),native_dist
      CHARACTER*7 beta(*)
      
*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,map,ii,jj
      REAL*8  cut2,xpi,ypi,zpi,xc,yc,zc,xd,yd,zd,rs,resi,resj
      LOGICAL ok
      INCLUDE 'pbc.h'

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      WRITE(knative,'(''# Native contact list with Reside-Residue''
     &     ,'' cutoff of '',f10.3)') native_dist
      map=0
      cut2=native_dist**2
      DO ii=1,nbun-1
         DO jj=ii+4,nbun
            ok=.FALSE.
            DO i=atres(1,ii),atres(2,ii)
               IF(beta(i)(1:1) .NE. 'h') THEN
                  resi=res(i)
                  xpi=xp0(i)
                  ypi=yp0(i)
                  zpi=zp0(i)
                  DO j=atres(1,jj),atres(2,jj)
                     IF(beta(j)(1:1) .NE. 'h') THEN
                        resj=res(j)
                        xc=xpi-xp0(j)
                        yc=ypi-yp0(j)
                        zc=zpi-zp0(j)
                        xc=xc-2.0D0*PBC(xc)
                        yc=yc-2.0D0*PBC(yc)
                        zc=zc-2.0D0*PBC(zc)
                        xd=co(1,1)*xc+co(1,2)*yc+co(1,3)*zc
                        yd=co(2,1)*xc+co(2,2)*yc+co(2,3)*zc
                        zd=co(3,1)*xc+co(3,2)*yc+co(3,3)*zc
                        rs=xd*xd+yd*yd+zd*zd
                        IF(rs .LT. cut2) THEN
c$$$                           IF(ii .EQ. 22 .AND. jj .EQ. 137) THEN
c$$$                              WRITE(*,*) i,j,DSQRT(rs)
c$$$                              xd=co(1,1)*xpi+co(1,2)*ypi+co(1,3)*zpi
c$$$                              yd=co(2,1)*xpi+co(2,2)*ypi+co(2,3)*zpi
c$$$                              zd=co(3,1)*xpi+co(3,2)*ypi+co(3,3)*zpi
c$$$                              WRITE(*,'(i8,3f15.6)') i,xd,yd,zd
c$$$                              xd=co(1,1)*xp0(j)+co(1,2)*yp0(j)+co(1,3)
c$$$     &                             *zp0(j)
c$$$                              yd=co(2,1)*xp0(j)+co(2,2)*yp0(j)+co(2,3)
c$$$     &                             *zp0(j)
c$$$                              zd=co(3,1)*xp0(j)+co(3,2)*ypi+co(3,3)*zpi
c$$$                              WRITE(*,'(i8,3f15.6)') i,xd,yd,zd
c$$$                           END IF
                           
                           ok=.TRUE.
                           GOTO 100
                        END IF
                     END IF
                  END DO
               END IF
            END DO
100         CONTINUE
            IF(ok) THEN
               map=map+1
               WRITE(knative,'(3i9)') map,ii,jj
            END IF
         END DO
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
