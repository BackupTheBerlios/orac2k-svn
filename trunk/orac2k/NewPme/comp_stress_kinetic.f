      SUBROUTINE comp_stress_kinetic(vpcmx,vpcmy,vpcmz,massp,co,nstart
     &     ,nend,volume,unitp,stressk,pressk)

************************************************************************
*   Time-stamp: <99/02/19 13:05:33 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Dec  8 1996 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8  vpcmx(*),vpcmy(*),vpcmz(*),massp(*),stressk(3,3),co(3,3)
     &     ,unitp,volume,pressk
      INTEGER nstart,nend

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,k
      REAL*8  vx,vy,vz,st(3,3)

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO i=1,3
         DO j=1,3
            stressk(i,j)=0.0D0
            st(i,j)=0.0D0
         END DO
      END DO

      DO i=nstart,nend
         vx=co(1,1)*vpcmx(i)+co(1,2)*vpcmy(i)+co(1,3)*vpcmz(i)
         vy=co(2,1)*vpcmx(i)+co(2,2)*vpcmy(i)+co(2,3)*vpcmz(i)
         vz=co(3,1)*vpcmx(i)+co(3,2)*vpcmy(i)+co(3,3)*vpcmz(i)
         st(1,1) = st(1,1)+massp(i)*vx*vpcmx(i)
         st(1,2) = st(1,2)+massp(i)*vx*vpcmy(i)
         st(1,3) = st(1,3)+massp(i)*vx*vpcmz(i)
         st(2,1) = st(2,1)+massp(i)*vy*vpcmx(i)
         st(2,2) = st(2,2)+massp(i)*vy*vpcmy(i)
         st(2,3) = st(2,3)+massp(i)*vy*vpcmz(i)
         st(3,1) = st(3,1)+massp(i)*vz*vpcmx(i)
         st(3,2) = st(3,2)+massp(i)*vz*vpcmy(i)
         st(3,3) = st(3,3)+massp(i)*vz*vpcmz(i)
      END DO

      DO i=1,3
         DO j=1,3
            DO k=1,3
               stressk(i,j) = stressk(i,j)+st(i,k)*co(j,k)
            END DO
         END DO
      END DO

      pressk = 0.0D0 
      DO i=1,3
         pressk=pressk+stressk(i,i)
      END DO
      pressk = pressk*unitp/(3.0D0*volume)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
