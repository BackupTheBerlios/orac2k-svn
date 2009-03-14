      SUBROUTINE correc_exp_scale(nstart,nend,nstart_cm,nend_cm
     &     ,cpress,ma,etap,tm4,vcax,vcay,vcaz,vpx,vpy,vpz,vco)

************************************************************************
*   Time-stamp: <02/09/07 15:00:07 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Feb 10 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER ma(*),nstart,nend,nstart_cm,nend_cm
      LOGICAL cpress,slv_exist,slt_exist
      REAL*8  tm4,vcax(*),vcay(*),vcaz(*),vpx(*),vpy(*),vpz(*),vco(3,3)
     &     ,etap(*)

*------------------------- LOCAL VARIABLES ----------------------------*


*----------------------- EXECUTABLE STATEMENTS ------------------------*

c----------  Correct cofm and barostat velocities using thermostat 1 ---
c --
      CALL correc_scale(nstart_cm,nend_cm,vcax,vcay,vcaz,etap(1),tm4)
               
      IF(cpress) THEN
         CALL correc_scale_co(9,vco,etap(1),tm4)
      END IF
      
c----------  Correct solute and solvent velocities using  --------------
c ---------  thermostat 2 and 3 

      CALL correc_scale_2(nstart,nend,ma,vpx,vpy,vpz,etap(2),tm4)
      
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
