      SUBROUTINE P_get_neighbors(nstart,nend,ngrp,xpg,ypg,zpg,co
     &     ,nnlpp,cut)

************************************************************************
*   Time-stamp: <99/03/09 09:22:41 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Mon Jan 18 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nstart,nend,ngrp
      INTEGER nnlpp(*)
      REAL*8  xpg(*),ypg(*),zpg(*),co(3,3),cut

*----------------------- VARIABLES IN COMMON --------------------------*


*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,map
      REAL*8  xpgi,ypgi,zpgi,xgg,ygg,zgg,xc,yc,zc,drj,cut2
      INCLUDE 'pbc.h'

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO i=1,ngrp
         nnlpp(i)=0
      END DO
      cut2=cut*cut
      DO i=nstart,nend
         xpgi=xpg(i)
         ypgi=ypg(i)
         zpgi=zpg(i)
         map=0
         DO j=i+1,ngrp
            xgg=xpgi-xpg(j)-2.0D0*PBC(xpgi-xpg(j))
            ygg=ypgi-ypg(j)-2.0D0*PBC(ypgi-ypg(j))
            zgg=zpgi-zpg(j)-2.0D0*PBC(zpgi-zpg(j))
            xc=co(1,1)*xgg+co(1,2)*ygg+co(1,3)*zgg
            yc=            co(2,2)*ygg+co(2,3)*zgg
            zc=                        co(3,3)*zgg
            drj=xc**2+yc**2+zc**2
            IF(drj .LT. cut2) THEN
               map=map+1
            END IF
         END DO
         nnlpp(i)=map
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
