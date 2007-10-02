      SUBROUTINE set_hhisto(hhisto_list,hhisto_bin,hhisto_dim
     &     ,hhisto_maxbin,iret,errmsg)

************************************************************************
*   Time-stamp: <98/04/10 15:37:26 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Fri Apr 10 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER hhisto_maxbin,hhisto_dim,hhisto_list(3,*)
      REAL*8  hhisto_bin
      INTEGER iret
      CHARACTER*80 errmsg
      
*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8  length,x1,x
      INTEGER i,GET_IND
      DATA length/5.0D0/
      GET_IND(x)=INT(x)+(SIGN(1.0D0,x-INT(x))-1.0D0)/2+1

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      iret=0
      errmsg=' '
      x1=length/hhisto_bin
      hhisto_dim=GET_IND(x1)
      IF(hhisto_dim .GT. hhisto_maxbin) THEN
         iret=1
         errmsg=
     &        ' Redimension _MAX_HHISTO_BIN_ in config.h. Dimensions '
     &        / /'too small.'
         RETURN
      END IF
      DO i=1,hhisto_dim
         hhisto_list(1,i)=0
         hhisto_list(2,i)=0
         hhisto_list(3,i)=0
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
