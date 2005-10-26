      SUBROUTINE comp_thermos_energy(n,ntot,t,qmass,eta,etap,uceh,hpot
     &     ,temph)

************************************************************************
*   Time-stamp: <99/02/18 15:43:21 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Apr  5 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER n,ntot(*)
      REAL*8  t,qmass(*),eta(*),etap(*),hpot,uceh,temph

*----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'unit.h'

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,m

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      m=0
      hpot=0.0D0
      uceh=0.0D0
      DO i=1,n
         hpot=hpot+(DBLE(ntot(i))+1.0D0)*boltz*t*eta(i)
         uceh=uceh+0.5D0*qmass(i)*etap(i)**2
         IF(DABS(qmass(i)) .GT. 1.0D-5) m=m+1
      END DO

      hpot=hpot/unite

      uceh=uceh*efact
      temph=2.0D0*uceh/(DBLE(m)*gascon)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
