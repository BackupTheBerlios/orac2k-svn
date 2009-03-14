      SUBROUTINE gen_abmd_kvect(vect,nvect,pvect,oc)

************************************************************************
*   Time-stamp: <98/06/14 16:05:46 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Aug 21 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER nvect
      REAL*8 vect,pvect(3,*),oc(3,3)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i
      REAL*8  u1,u2,u3,w1,w2,w3,w,ranf,sumx,sumy,sumz

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO i=1,nvect
100      CONTINUE
         u1=ranf()
         u2=ranf()
         u3=ranf()
         w1=1.0D0-2.0D0*u1
         w2=1.0D0-2.0D0*u2
         w3=1.0D0-2.0D0*u3
         w=w1**2+w2**2+w3**2
         IF(w .LT. 1.0D0) THEN
            sumx=vect*w1/DSQRT(w)
            sumy=vect*w2/DSQRT(w)
            sumz=vect*w3/DSQRT(w)
            pvect(1,i)=sumx
            pvect(2,i)=sumy
            pvect(3,i)=sumz
         ELSE
            GOTO 100
         END IF
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
