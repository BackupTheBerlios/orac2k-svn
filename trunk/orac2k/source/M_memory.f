      SUBROUTINE M_memory(ip_point,length)

************************************************************************
*   Time-stamp: <1999-10-12 18:33:07 marchi>                             *
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

#elif defined AXP
      REAL*8  a(*)
      INTEGER*8 ip_point,malloc
      POINTER (ip_a,a)
#else
      REAL*8  a(*)
      INTEGER ip_point,malloc
      POINTER (ip_a,a)
#endif
      INTEGER length

*------------------------- LOCAL VARIABLES ----------------------------*

#if defined _CRAY_ | defined T3E
      INTEGER*8 len,ierr1,one
#else
      INTEGER ierr1,one,len
      INTEGER elsize,offset
#endif
      DATA one/1/

*----------------------- EXECUTABLE STATEMENTS ------------------------*

#if defined DYNAMIC_MEM
#if defined _CRAY_ | defined T3E
      len=length
      CALL hpalloc(ip_point,len,ierr1,one)
      CALL hpcheck(ierr1)
      IF(ierr1 .NE. 0) THEN
         WRITE(6,'('' Heap memory has been corrupted. '')')
      END IF
#else
      len=length+1
      len=MAX(len,400)
#ifdef AIX
      ip_a=malloc(%val(len))
#elif defined AXP | defined OSF1
      ip_a=malloc(len)
#endif
      ip_point=loc(a)
      IF(ip_point .EQ. 0) THEN
         WRITE(6,'('' Heap memory has been corrupted. '')')
         STOP
      END IF
#endif
#endif

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
