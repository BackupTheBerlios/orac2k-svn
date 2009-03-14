      SUBROUTINE memory(ip_point,length,nbyte)

************************************************************************
*   Time-stamp: <2009-03-09 12:44:54 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Oct 21 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

#if defined _CRAY_ | defined T3E
      INTEGER*8 ip_point
#elif defined AXP
      INTEGER*8 ip_point
#else
      INTEGER ip_point
#endif
      INTEGER length,nbyte

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER len,M_get_length

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      len=M_get_length(length,nbyte)
#if defined T3E & defined PARALLEL
      CALL M_memory_s(ip_point,len)
#else
      CALL M_memory(ip_point,len)
#endif

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
