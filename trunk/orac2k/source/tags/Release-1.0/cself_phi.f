      SUBROUTINE cself_phi(nstart,nend,alphal,charge,volume,phi)
      
************************************************************************
*                                                                      *
*                                                                      *
*     Compute the Ewald self term.                                     *
*                                                                      *
*---- Last update 06/12/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNALS  NONE                                                  *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER nstart,nend
      REAL*8  charge(*),phi(*),alphal,volume

*------------------ VARIABLES IN COMMON --------------------------------

      INCLUDE 'unit.h'

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER n
      REAL(8) :: TotCharge

*==================== EXECUTABLE STATEMENTS ============================


*=======================================================================
*---- Compute the self term --------------------------------------------
*=======================================================================
      
      TotCharge=0.0D0
      DO n=nstart,nend
         phi(n)=phi(n)-2.0D0*charge(n)*alphal/DSQRT(pi)
         TotCharge=TotCharge+charge(n)
      END DO
      phi(nstart:nend)=phi(nstart:nend)-TotCharge*pi/alphal**2/volume

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
