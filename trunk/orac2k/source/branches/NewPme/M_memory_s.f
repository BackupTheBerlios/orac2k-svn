      SUBROUTINE M_memory_s(ip_point,length)

************************************************************************
*   Time-stamp: <99/02/13 14:16:19 marchi>                             *
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
#else
      INTEGER ip_point
#endif
      INTEGER length

*------------------------- LOCAL VARIABLES ----------------------------*

#if defined _CRAY_ | defined T3E
      INTEGER*8 len,ierr1,one
#else
      INTEGER len,ierr1,one
#endif
      DATA one/1/

*----------------------- EXECUTABLE STATEMENTS ------------------------*

#if defined _CRAY_ | defined T3E
      len=length
      CALL shpalloc(ip_point,len,ierr1,one)
      IF(ierr1 .NE. 0) THEN
         WRITE(6,'('' Symmetric heap memory has been corrupted. '')')
      END IF
      CALL hpcheck(ierr1)
#endif

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
