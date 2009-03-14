      SUBROUTINE dupl_real2(array,na,n1,ndim,nrep,iret)

************************************************************************
*   Time-stamp: <95/03/14 15:23:51 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Mar 14 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER na,n1,ndim,nrep,iret
      REAL*8  array(n1,*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER n,m,l,k

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO n=2,nrep
         DO m=1,na
            l=na*(n-1)+m
            IF(l .GT. ndim) THEN
               iret=1
               RETURN
            END IF
            DO k=1,n1
               array(k,l)=array(k,m)
            END DO
         END DO
      END DO
      
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
