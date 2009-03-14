      SUBROUTINE close_files

************************************************************************
*   Time-stamp: <00/10/30 15:41:52 sterpone>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Aug  6 1996 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'unit.h'

*------------------------- LOCAL VARIABLES ----------------------------*


*----------------------- EXECUTABLE STATEMENTS ------------------------*

      
      IF(kgdata .NE. 0) THEN
         CLOSE(kgdata)
      END IF
      IF(kcnfig_i .NE. 0) THEN
         CLOSE(kcnfig_i)
      END IF
      IF(kcnfig_o .NE. 0) THEN
         CLOSE(kcnfig_o)
      END IF
      IF(knlist .NE. 0) THEN
         CLOSE(knlist)
      END IF
      IF(kgofr .NE. 0) THEN
         CLOSE(kgofr)
      END IF
      IF(kplot .NE. 0) THEN
         CLOSE(kplot)
      END IF
      IF(kconf .NE. 0) THEN
         CLOSE(kconf)
      END IF
      IF(kpot .NE. 0) THEN
         CLOSE(kpot)
      END IF
      IF(kout .NE. 0) THEN
         CLOSE(kout)
      END IF
      IF(kfield .NE. 0) THEN
         CLOSE(kfield)
      END IF
      IF(kvaf .NE. 0) THEN
         CLOSE(kvaf)
      END IF
      IF(krms .NE. 0) THEN
         CLOSE(krms)
      END IF
      IF(kxrms .NE. 0) THEN
         CLOSE(kxrms)
      END IF
      IF(kxrms_atm .NE. 0) THEN
         CLOSE(kxrms_atm)
      END IF
      IF(ksol1 .NE. 0) THEN
         CLOSE(ksol1)
      END IF
      IF(ksol2 .NE. 0) THEN
         CLOSE(ksol2)
      END IF
      IF(kgroup .NE. 0) THEN
         CLOSE(kgroup)
      END IF
      IF(ktop .NE. 0) THEN
         CLOSE(ktop)
      END IF
      IF(ktemplate .NE. 0) THEN
         CLOSE(ktemplate)
      END IF
      IF(ktest .NE. 0) THEN
         CLOSE(ktest)
      END IF
      IF(iurest .NE. 0) THEN
         CLOSE(iurest)
      END IF
      IF(kspec .NE. 0) THEN
         CLOSE(kspec)
      END IF
      IF(kgr .NE. 0) THEN
         CLOSE(kgr)
      END IF
      IF(kplot_fragm .NE. 0) THEN
         CLOSE(kplot_fragm)
      END IF
      IF(kgofr_sk .NE. 0) THEN
         CLOSE(kgofr_sk)
      END IF
      
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
