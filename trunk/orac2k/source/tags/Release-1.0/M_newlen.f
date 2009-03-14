      SUBROUTINE M_newlen(ip_point,length)

************************************************************************
*   Time-stamp: <99/02/12 17:23:18 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed Feb 10 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

#if defined _CRAY_ | defined T3E
      INTEGER*8 ip_point
      INTEGER length
#else
      INTEGER ip_point
      INTEGER length
#endif

*------------------------- LOCAL VARIABLES ----------------------------*

#if defined _CRAY_ | defined T3E
      INTEGER*8 len,ierr1
#else
      INTEGER len,ierr1
#endif

*----------------------- EXECUTABLE STATEMENTS ------------------------*

#if defined _CRAY_ | defined T3E
      len=length
      CALL hpclmove(ip_point,len,ierr1,1)
#endif

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
