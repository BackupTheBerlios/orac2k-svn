      SUBROUTINE calccofm(kout,latms,natms,patms,wca,xp_ini,
     &     yp_ini,zp_ini,qt,xp0,yp0,zp0,mass,nato,fstep)

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

      REAL*8  wca(*),xp0(*),yp0(*),zp0(*),xp_ini(*),yp_ini(*),zp_ini(*)
      real*8  mass(*),qt(4),fstep
      INTEGER natms,patms(2,*),nato,kout
      LOGICAL latms

*----------------------- VARIABLES IN COMMON --------------------------*

      REAL*8  xyz(3,m1),xyz0(3,m1),xyzfit(3,m1),wca2(m1),work(m1)
      COMMON /rag1/ xyz,xyz0,xyzfit,work,wca2

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8 err,sum,cofm(3),mtot,dist
      INTEGER n,i,j

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      sum=0.0D0
      DO i=1,nato
         xyz(1,i)=xp0(i)
         xyz(2,i)=yp0(i)
         xyz(3,i)=zp0(i)
         xyz0(1,i)=xp_ini(i)
         xyz0(2,i)=yp_ini(i)
         xyz0(3,i)=zp_ini(i)
         wca2(i)=wca(i)
         sum = sum + wca2(i)
      END DO

      CALL normal(wca2,nato)

      IF(DABS(sum) .GT. 1.0D-3) THEN
         CALL xfit(xyz0,xyz,xyzfit,qt,wca2,work,nato,err)
         
         write(kout,'(a7,2x,f12.3)') 'Tstep =',fstep 
         DO n=1,natms
            mtot = 0.0d0
            DO j=1,3
               cofm(j) = 0.0d0
            END DO

            DO i=patms(1,n),patms(2,n)
               DO j=1,3
                  cofm(j) = cofm(j) + mass(i)*xyzfit(j,i)
               END DO
               mtot = mtot + mass(i)
            END DO

            dist = 0.0d0
            DO j=1,3
               cofm(j) = cofm(j)/mtot
               dist = dist + cofm(j)**2
            END DO
            dist = dsqrt(dist)
          write(kout,'(i6,4(2x,f8.3))') patms(1,n),dist,(cofm(j),j=1,3)
         END DO

      ELSE
         write(kout,*) 'No atoms to fit'
      END IF

         
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END

