      SUBROUTINE write_contacts(kprint,fstep,listp,list,listp_o,list_o)

************************************************************************
*   Time-stamp: <99/03/14 17:03:56 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Mar 14 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER kprint,listp,list(2,*)
      REAL*8  fstep,listp_o,list_o(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,ii,jj

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      WRITE(kprint,'('' T ='',f11.3,f15.8)') fstep,listp_o

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
