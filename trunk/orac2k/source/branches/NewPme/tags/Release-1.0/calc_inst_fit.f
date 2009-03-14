      SUBROUTINE calc_inst_fit(anprot,annpro,anpoint,protl,wca,xp_ini,
     &     yp_ini,zp_ini,xpo,ypo,zpo,qt,xp0,yp0,zp0,nato)

************************************************************************
*   Time-stamp: <98/03/20 15:32:27 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Jun 29 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none
      INCLUDE  'parst.h'

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8  wca(*),xp0(*),yp0(*),zp0(*),xp_ini(*),yp_ini(*)
     &     ,zp_ini(*),xpo(*),ypo(*),qt(4)
     &     ,zpo(*)
      INTEGER annpro,anpoint(2,*),protl(*),nato
      LOGICAL anprot

*----------------------- VARIABLES IN COMMON --------------------------*

      REAL*8  xyz(3,m1),xyz0(3,m1),xyzfit(3,m1),wca2(m1),work(m1)
      COMMON /rag1/ xyz,xyz0,xyzfit,work,wca2

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8 err,sum
      INTEGER n,m,l,i,ncount,nc2

*----------------------- EXECUTABLE STATEMENTS ------------------------*

*=======================================================================
*----- Calculate XRMS from template ------------------------------------
*=======================================================================

      IF(anprot) THEN
         DO n=1,annpro
            m=anpoint(2,n)-anpoint(1,n)+1
            DO l=1,m
               i=anpoint(1,n)-1+l
               xyz(1,l)=xp0(i)
               xyz(2,l)=yp0(i)
               xyz(3,l)=zp0(i)
               xyz0(1,l)=xp_ini(i)
               xyz0(2,l)=yp_ini(i)
               xyz0(3,l)=zp_ini(i)
               wca2(l)=wca(i)
            END DO
            sum=0.0D0
            DO i=1,m
               sum = sum + wca2(i)
            END DO
            CALL normal(wca2,m)
            IF(DABS(sum) .GT. 1.0D-3) THEN
               CALL xfit(xyz0,xyz,xyzfit,qt,wca2,work,m,err)
               DO l=1,m
                  i=anpoint(1,n)-1+l
                  xpo(i)=xyzfit(1,l)
                  ypo(i)=xyzfit(2,l)
                  zpo(i)=xyzfit(3,l)
               END DO
            END IF
         END DO

      ELSE

         n=0
         nc2=0
         ncount=1
100      CONTINUE
         nc2=nc2+1
         m=protl(ncount)
         DO l=1,m
            i=protl(ncount+l)
            xyz(1,l)=xp0(i)
            xyz(2,l)=yp0(i)
            xyz(3,l)=zp0(i)
            xyz0(1,l)=xp_ini(i)
            xyz0(2,l)=yp_ini(i)
            xyz0(3,l)=zp_ini(i)
            wca2(l)=wca(i)
         END DO
         sum=0.0D0
         DO l=1,m
            i=protl(ncount+l)
            sum=sum+wca(i)
         END DO
         CALL normal(wca2,m)
         IF(DABS(sum) .GT. 1.0D-3) THEN
            CALL xfit(xyz0,xyz,xyzfit,qt,wca2,work,m,err)
            DO l=1,m
               i=protl(ncount+l)
                  xpo(i)=xyzfit(1,l)
                  ypo(i)=xyzfit(2,l)
                  zpo(i)=xyzfit(3,l)
            END DO
         END IF
         ncount=ncount+m+1
         n=n+m
         IF(n .LT. nato) GOTO 100
      END IF
         
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*
         
      RETURN
      END

