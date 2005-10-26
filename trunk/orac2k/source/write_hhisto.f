      SUBROUTINE write_hhisto(hhisto_list,hhisto_bin,hhisto_dim
     &     ,hhisto_count)

************************************************************************
*   Time-stamp: <98/06/25 10:36:48 marchi>                             *
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

      INCLUDE 'unit.h'

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER hhisto_dim,hhisto_count,hhisto_list(3,*)
      REAL*8  hhisto_bin

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i
      REAL*8  x,y

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      REWIND(khhisto)
      WRITE(khhisto,11000)
      DO i=1,hhisto_dim
         x=DBLE(i)*hhisto_bin
         y=DBLE(hhisto_list(1,i))/DBLE(hhisto_count)/hhisto_bin
         WRITE(khhisto,10000) x,y
      END DO
      WRITE(khhisto,12000)
      DO i=1,hhisto_dim
         x=DBLE(i)*hhisto_bin
         y=DBLE(hhisto_list(2,i))/DBLE(hhisto_count)/hhisto_bin
         WRITE(khhisto,10000) x,y
      END DO
      WRITE(khhisto,13000)
      DO i=1,hhisto_dim
         x=DBLE(i)*hhisto_bin
         y=DBLE(hhisto_list(3,i))/DBLE(hhisto_count)/hhisto_bin
         WRITE(khhisto,10000) x,y
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

10000 FORMAT(f12.3,e16.6)
11000 FORMAT('#  Dist.    No. Contacts    SLT-SLT')
12000 FORMAT('#  Dist.    No. Contacts    SLT-SLV')
13000 FORMAT('#  Dist.    No. Contacts    SLV-SLV')
      RETURN
      END
