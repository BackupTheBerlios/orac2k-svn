      SUBROUTINE covfod(spring,abmd_tors)

************************************************************************
*                                                                      *
*     Convert folding parameters to program units.                     *
*                                                                      *
*---- Last update 04/01/91 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi UC Berkeley 1991                       *
*                                                                      *
*     EXTERNALS NONE                                                   *
*                                                                      *
*                                                                      *
************************************************************************

!======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      REAL*8  spring
      LOGICAL abmd_tors

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i

*-------------------- VARIABLES IN COMMON ------------------------------

      INCLUDE 'unit.h'


!==================== EXECUTABLE STATEMENTS ============================

      IF(abmd_tors) THEN
         spring=1000.0D0*spring/(unite*avogad)
      ELSE
         spring=1000.0D0*spring/(unite*avogad)
      END IF

!================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
