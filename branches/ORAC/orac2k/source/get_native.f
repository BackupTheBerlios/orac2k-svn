      SUBROUTINE get_native(knative,listp,list,ma,atres,beta,map1,map2
     &     ,iret,errmsg)

************************************************************************
*   Time-stamp: <01/08/16 10:30:56 marchi>                             *
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

      INTEGER knative,ma,listp,list(2,*),iret,map1(*),map2(*),atres(2,*)
      CHARACTER*7 beta(*)
      CHARACTER*80 errmsg

*----------------------- VARIABLES IN COMMON --------------------------*

*------------------------- LOCAL VARIABLES ----------------------------*

      CHARACTER*80 line
      INTEGER map,i,j,ii,jj,count1,count2,m1,m2,ia

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      iret=0
100   READ(knative,'(a78)',END=600) line(1:78)
      IF(line(1:1) .EQ. '#') GOTO 100
      IF(line .EQ. ' ') GOTO 100
      READ(line,'(3i9)') map,list(1,map),list(2,map)
      IF(map .GT. ma) THEN
         iret=1
         errmsg=
     &        ' Physical dimensions of the native list exceeded.'
         RETURN
      END IF
      GOTO 100

600   CONTINUE
      listp=map

      count1=0
      count2=0
      DO ia=1,listp
         ii=list(1,ia)
         m1=0
         DO i=atres(1,ii),atres(2,ii)
            IF(beta(i)(1:1) .NE. 'h') THEN
               m1=m1+1
               map1(count1+1+m1)=i
            END IF
         END DO
         map1(count1+1)=m1
         jj=list(2,ia)
         m2=0
         DO j=atres(1,jj),atres(2,jj)
            IF(beta(j)(1:1) .NE. 'h') THEN
               m2=m2+1
               map2(count2+1+m2)=j
            END IF
         END DO
         map2(count2+1)=m2
         count1=count1+1+m1
         count2=count2+1+m2
      END DO


*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
