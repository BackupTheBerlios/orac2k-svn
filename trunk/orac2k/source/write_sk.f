      SUBROUTINE write_sk(ksofk,fstep,sofk_cut,sofk_delta,sofk,nsofk
     &     ,maxsk)

************************************************************************
*   Time-stamp: <1998-11-16 15:18:34 marchi>                           *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Mon Nov 16 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER ksofk,maxsk
      REAL*8  sofk(maxsk,2),nsofk(maxsk,2),sofk_delta,sofk_cut,fstep

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER kbox,i
      REAL*8  aux

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      REWIND ksofk
      WRITE(ksofk,'(/''# Time Step ='',f15.5)') fstep
      kbox=MIN0(IDINT(sofk_cut/sofk_delta+0.5D0),maxsk)
      DO i=1,kbox
         IF(nsofk(i,2) .NE. 0) THEN
            aux=sofk(i,2)/nsofk(i,2)
            WRITE(ksofk,'(f10.4,e16.7)') DFLOAT(i)*sofk_delta,aux
         END IF
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
