      SUBROUTINE comp_contacts(co,xp0,yp0,zp0,beta,alpha,atres,listp
     &     ,list,native_dist,listp_o,list_o)

************************************************************************
*   Time-stamp: <99/03/14 17:29:09 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Mar 14 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER kprint,atres(2,*),listp,list(2,*)
      REAL*8  listp_o,list_o(*)
      REAL*8  xp0(*),yp0(*),zp0(*),co(3,3),native_dist,alpha
      CHARACTER*7 beta(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,ii,jj,ia
      REAL*8  cut2,xpi,ypi,zpi,xc,yc,zc,xd,yd,zd,rs
      LOGICAL ok
      REAL*8 r,sum,nat,aux
      INCLUDE 'pbc.h'
      NAT(r)=0.5D0*(1.0D0-(DEXP(-r)-DEXP(r))/(DEXP(-r)+DEXP(r)))

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      cut2=native_dist**2
      DO ia=1,listp
         ii=list(1,ia)
         jj=list(2,ia)
         ok=.FALSE.
         sum=1.0D0
         DO i=atres(1,ii),atres(2,ii)
            IF(beta(i)(1:1) .NE. 'h') THEN
               xpi=xp0(i)
               ypi=yp0(i)
               zpi=zp0(i)
               DO j=atres(1,jj),atres(2,jj)
                  IF(beta(j)(1:1) .NE. 'h') THEN
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
                     r=(DSQRT(rs)-native_dist)*alpha
                     aux=NAT(r)
                     sum=sum*aux
                  END IF
               END DO
            END IF
         END DO
         list_o(ia)=1.0D0-sum
      END DO
      listp_o=0.0D0
      DO i=1,listp
         listp_o=listp_o+DNINT(list_o(i))
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
