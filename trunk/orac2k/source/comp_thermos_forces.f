      SUBROUTINE comp_thermos_forces(nprocs,nstart,nend,nstart_cm
     &     ,nend_cm,ndf_thermos,ss_index,co,vpx,vpy,vpz,vcmx,vcmy
     &     ,vcmz,mass,masscm,t,fth)

************************************************************************
*   Time-stamp: <99/02/18 18:36:46 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Apr  5 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nprocs,nstart,nend,nstart_cm,nend_cm,ndf_thermos(*)
     &     ,ss_index(*)
      REAL*8  t,vpx(*),vpy(*),vpz(*),vcmx(*),vcmy(*),vcmz(*)
     &     ,mass(*),masscm(*),fth(*),co(3,3)

*----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'unit.h'

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,k,j
      REAL*8  tto(3),xd,yd,zd

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      tto(1)=0.0D0
      tto(2)=0.0D0
      tto(3)=0.0D0
      DO i=nstart_cm,nend_cm
         xd=co(1,1)*vcmx(i)+co(1,2)*vcmy(i)+co(1,3)*vcmz(i)
         yd=                co(2,2)*vcmy(i)+co(2,3)*vcmz(i)
         zd=                                co(3,3)*vcmz(i)
         tto(1)=tto(1)+masscm(i)*(xd**2+yd**2+zd**2)
      END DO

      DO i=nstart,nend
         k=ss_index(i)+1
         tto(k)=tto(k)+mass(i)*(vpx(i)**2+vpy(i)**2+vpz(i)**2)
      END DO
      
      DO i=1,3
         fth(i)=tto(i)-DFLOAT(ndf_thermos(i))*boltz*t/unite
     &        /DFLOAT(nprocs)
       END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
