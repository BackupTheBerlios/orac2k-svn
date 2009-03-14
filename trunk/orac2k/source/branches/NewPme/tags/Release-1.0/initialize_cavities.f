      SUBROUTINE initialize_cavities(ss_index,ntap,res,bin_size_cav
     &     ,rmax_size_cav,cavity_h,maxcav_bin,maxcav_nres,jmax_cav,iret
     &     ,errmsg)

************************************************************************
*   Time-stamp: <01/07/25 15:09:43 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Mon Apr 13 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER maxcav_bin,maxcav_nres,iret,jmax_cav,ntap
      INTEGER res(*),ss_index(*)
      REAL*8  cavity_h(maxcav_bin,*),bin_size_cav,rmax_size_cav
      CHARACTER*80 errmsg

*----------------------- VARIABLES IN COMMON --------------------------*



*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER nmax,i,j,GET_IND
      REAL*8  dx,y
      GET_IND(y)=INT(y)+(SIGN(1.0D0,y-INT(y))-1.0D0)/2+1

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      iret=0
      dx=bin_size_cav
      nmax=GET_IND(rmax_size_cav/dx)
      IF(nmax+1 .GT. maxcav_bin) THEN
         iret=1
         errmsg=
     &        ' Dimension of _MAX_CAVITIES_BIN_ too small.'
     &     / /' Redimension config.h.'
         RETURN
      END IF
      jmax_cav=0
      DO i=1,ntap
         IF(ss_index(i) .EQ. 1) THEN
            j=res(i)
            IF(jmax_cav .LT. j) jmax_cav=j
         END IF
      END DO
      jmax_cav=jmax_cav+1
      IF(jmax_cav .GT. maxcav_nres) THEN
         iret=1
         WRITE(errmsg,'('' Dimension of _MAX_CAVITIES_NRES_ '',
     &        ''too small. need dimension .GE. to '',i8)') jmax_cav
         RETURN
      END IF
      DO i=1,maxcav_bin
         DO j=1,maxcav_nres
            cavity_h(i,j)=0.0D0
         END DO
      END DO
      
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
