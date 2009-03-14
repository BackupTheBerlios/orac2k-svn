      SUBROUTINE get_atres(atres,res,nato,nbun)

************************************************************************
*   Time-stamp: <01/08/14 17:41:51 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Mar 13 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER atres(2,*),nato,res(*),nbun
      
*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,mao,res_o,res_n,map

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      res_o=0
      map=0
      DO i=1,nato
         res_n=res(i)
         IF(res_o .NE. res_n) THEN
            IF(i .NE. nato) THEN
               map=map+1
               atres(1,map)=i
            END IF
            IF(i .NE. 1) atres(2,map-1)=i-1
            res_o=res_n
         END IF
      END DO
      atres(2,map)=nato
      nbun=map

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*
      
      RETURN
      END
