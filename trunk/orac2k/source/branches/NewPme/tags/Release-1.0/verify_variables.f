      SUBROUTINE verify_variables(iret,errmsg)

************************************************************************
*   Time-stamp: <2005-03-13 22:17:32 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Jun 20 1996 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER iret
      CHARACTER*80 errmsg

*----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'parst.h'
      INCLUDE 'cpropar.h'
      INCLUDE 'unit.h'

*------------------------- LOCAL VARIABLES ----------------------------*


*----------------------- EXECUTABLE STATEMENTS ------------------------*
      

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

100   format(//'in MTSMD: ******** DISASTER ERROR ************'
     &    /'number of group too large. Cannot address integer',I7,'.',
     &    /'action: Change neighbor list in mtsmd to integer*4 or run',
     &    /'        with no RESPA') 
200   format(//'in MTSMD: ******** DISASTER ERROR ************'
     &    /'number of solvent molecules too large. Cannot address',
     &     ' integer',I7,'.',
     &    /'action: Change neighbor list in mtsmd to integer*4 or run',
     &    /'        with no RESPA') 
      RETURN
      END
