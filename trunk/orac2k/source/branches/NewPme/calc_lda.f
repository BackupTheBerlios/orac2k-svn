      SUBROUTINE calc_lda(kout,nlda_mol,nlda_atm,nlda_zero,wca,xp_ini,
     &     yp_ini,zp_ini,qt,xp0,yp0,zp0,mass,nato,fstep)

************************************************************************
*                                                                      *
*              Author:  Matteo Ceccarelli                              *
*              CECAM/ENS Lyon, FRANCE                                  *
*                                                                      *
*              - May 1999 -                                            *

************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none
      INCLUDE  'parst.h'

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8  wca(*),xp0(*),yp0(*),zp0(*),xp_ini(*),yp_ini(*),zp_ini(*)
      real*8  mass(*),qt(4),fstep
      INTEGER nlda_mol,nlda_atm,nlda_zero,nato,kout

*----------------------- VARIABLES IN COMMON --------------------------*

      REAL*8  xyz(3,m1),xyz0(3,m1),xyzfit(3,m1),wca2(m1),work(m1)
      COMMON /rag1/ xyz,xyz0,xyzfit,work,wca2

*------------------------- LOCAL VARIABLES ----------------------------*
      integer LWORK
      parameter(LWORK=9)

      REAL*8 err,sum,cofm(3),mtot,dist,cofm0(3)
      REAL*8 mofi(3),INZ(3,3),VAUX(LWORK),norm2,mvect(3,3)
      REAL*8 axex,axey,axez
      INTEGER n,i,j,k,idx,INFO

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

         DO j=1,3
            cofm0(j) = 0.0d0
            cofm(j) = 0.0d0
         END DO

c------  calculate RCs cofm
         mtot = 0.0d0

         do i=1,nlda_zero
              DO j=1,3
                 cofm0(j) = cofm0(j) + mass(i)*xyzfit(j,i)
              END DO
              mtot = mtot + mass(i)
         end do

         DO j=1,3
            cofm0(j) = cofm0(j)/mtot
         END DO

c------ new coordinates
         do i=1,nato
            do j=1,3
               xyzfit(j,i) = xyzfit(j,i) - cofm0(j)
            end do
         end do

c------ calculate LDAs cofm
         mtot = 0.0d0

         DO k=1,nlda_mol
            do i=1,nlda_atm
               idx = nlda_zero + (k-1)*nlda_atm + i
               DO j=1,3
                  cofm(j) = cofm(j) + mass(idx)*xyzfit(j,idx)
               END DO
               mtot = mtot + mass(idx)
            end do
         END DO

         dist = 0.0d0
         DO j=1,3
            cofm(j) = cofm(j)/mtot
            dist = dist + cofm(j)**2
         END DO
         dist = dsqrt(dist)

         write(kout,'(a10,4(2x,f8.3))') 
     &              'COFM-RC   ',dist,(cofm(j),j=1,3)

c------ LDA
     
      DO i=nlda_zero+1,nlda_zero+nlda_mol*nlda_atm
         do j=1,3
            xyzfit(j,i) = xyzfit(j,i) - cofm(j)
         end do
      END DO

      do j=1,3
      do k=1,3
         INZ(j,k) = 0.0d0
      end do
      end do

      do i=nlda_zero+1,nlda_zero+nlda_mol*nlda_atm
         norm2 = xyzfit(1,i)**2+xyzfit(2,i)**2+xyzfit(3,i)**2
         do j=1,3
            do k=1,3
               INZ(j,k) = INZ(j,k) - mass(i)*xyzfit(j,i)*xyzfit(k,i)
            end do
            INZ(j,j) = INZ(j,j) + mass(i)*norm2
         end do
      end do

      INZ(1,2) = 0.50d0*( INZ(1,2) + INZ(2,1) )
      INZ(1,3) = 0.50d0*( INZ(1,3) + INZ(3,1) )
      INZ(2,3) = 0.50d0*( INZ(2,3) + INZ(3,2) )
      INZ(2,1) =  INZ(1,2)
      INZ(3,1) =  INZ(1,3)
      INZ(3,2) =  INZ(2,3)


      call eigrs(INZ,3,11,mofi,mvect,3,VAUX,INFO)

      axex = dsqrt(5.0d0/mtot)*dsqrt(mofi(1)+mofi(3)-mofi(2))
      axey = dsqrt(5.0d0/mtot)*dsqrt(mofi(2)+mofi(3)-mofi(1))
      axez = dsqrt(5.0d0/mtot)*dsqrt(mofi(1)+mofi(2)-mofi(3))

         write(kout,'(a7,3(2x,f8.3))') 'INZ    ',(mofi(j)/10.d+6,j=1,3)
c         write(kout,'(a7,2(2x,f8.3))') 'INZ    ',mofi(3)/mofi(2),
c     &                                           mofi(2)/mofi(1)
c         write(kout,'(9f8.4)') ((mvect(i,j),j=1,3),i=1,3)
c         write(kout,'(a7,3(2x,f8.3))') 'AXES   ', axex,axey,axez

      ELSE
         write(kout,*) 'No atoms to fit'
      END IF

         
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END

