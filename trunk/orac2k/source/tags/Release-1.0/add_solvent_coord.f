      SUBROUTINE add_solvent_coord(xp0,yp0,zp0,xpa,ypa,zpa,ntap,nmol
     &     ,nato_slv,n1,iret,errmsg)

************************************************************************
*   Time-stamp: <2005-02-26 11:34:10 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed Feb 12 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8 xp0(*),yp0(*),zp0(*),xpa(*),ypa(*),zpa(*)
      INTEGER ntap,nato_slv,nmol,iret,n1
      CHARACTER*80 errmsg

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER nts,n

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      nts=nato_slv*nmol
      IF(n1 .LT. nts+ntap) THEN
         iret=1
         errmsg=
     &'After solvent generation: No. of atoms exceeds physical'
     &//' dimensions. Change config.h and recmp.'
         WRITE(*,*) nts,nato_slv,nmol,ntap
         RETURN

      END IF
      DO n=1,nts
         xp0(ntap+n)=xpa(n)
         yp0(ntap+n)=ypa(n)
         zp0(ntap+n)=zpa(n)
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
