      SUBROUTINE comp_ang_vel(start,sp,mass,vx,vy,vz,x0,y0,z0,omega
     &     ,iner,i_iner,erot)

************************************************************************
*   Time-stamp: <99/10/06 23:08:09 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Fri Oct  1 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER sp(*)
      LOGICAL start
      REAL*8  mass(*),vx(*),vy(*),vz(*),x0(*),y0(*),z0(*)
     &     ,omega(3),iner(3,3),i_iner(3,3),erot

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER ii,i,j,lh(3),mh(3),m
      REAL*8  det,xc,yc,zc,lm(3),am(3),cm(3),massa,tmass

*----------------------- EXECUTABLE STATEMENTS ------------------------*

c$$$  
c$$$  Compute total linear moment and center of mass
c$$$      
      
      DO i=1,3
         lm(i)=0.0D0
      END DO
      m=sp(1)
      DO ii=1,m
         i=ii
         IF(start) i=sp(ii+1)
         massa=mass(i)
         lm(1)=lm(1)+massa*vx(i)
         lm(2)=lm(2)+massa*vy(i)
         lm(3)=lm(3)+massa*vz(i)
      END DO
      DO i=1,3
         cm(i)=0.0D0
      END DO
      tmass=0.0D0
      DO ii=1,m
         i=ii
         IF(start) i=sp(ii+1)
         massa=mass(i)
         xc=x0(i)
         yc=y0(i)
         zc=z0(i)
         cm(1)=cm(1)+massa*xc
         cm(2)=cm(2)+massa*yc
         cm(3)=cm(3)+massa*zc
         tmass=tmass+massa
      END DO
      DO i=1,3
         cm(i)=cm(i)/tmass
      END DO
c$$$  
c$$$  Compute total angular momentum  and tensor of inertia
c$$$      
      DO i=1,3
         iner(i,1)=0.0D0
         iner(i,2)=0.0D0
         iner(i,3)=0.0D0
      END DO
      DO i=1,3
         am(i)=0.0D0
      END DO
      DO ii=1,m
         i=ii
         IF(start) i=sp(ii+1)
         massa=mass(i)
         xc=x0(i)
         yc=y0(i)
         zc=z0(i)
         am(1)=am(1)+massa*(yc*vz(i)-zc*vy(i))
         am(2)=am(2)+massa*(zc*vx(i)-xc*vz(i))
         am(3)=am(3)+massa*(xc*vy(i)-yc*vx(i))
         iner(1,1)=iner(1,1)+massa*(yc**2+zc**2)
         iner(2,2)=iner(2,2)+massa*(zc**2+xc**2)
         iner(3,3)=iner(3,3)+massa*(xc**2+yc**2)
         iner(1,2)=iner(1,2)-massa*xc*yc
         iner(1,3)=iner(1,3)-massa*xc*zc
         iner(2,3)=iner(2,3)-massa*yc*zc
      END DO
      iner(2,1)=iner(1,2)
      iner(3,1)=iner(1,3)
      iner(3,2)=iner(2,3)

      DO i=1,3
         i_iner(i,1)=iner(i,1)
         i_iner(i,2)=iner(i,2)
         i_iner(i,3)=iner(i,3)
      END DO

c$$$  
c$$$  Subtract linear momentum contribution
c$$$      

      am(1)=am(1)-(cm(2)*lm(3)-cm(3)*lm(2))
      am(2)=am(2)-(cm(3)*lm(1)-cm(1)*lm(3))
      am(3)=am(3)-(cm(1)*lm(2)-cm(2)*lm(1))

c$$$  
c$$$  Invert the tensor of inertia and compute angular velocity
c$$$      

      CALL minv(3,i_iner,det,lh,mh)
      omega(1)=i_iner(1,1)*am(1)+i_iner(1,2)*am(2)+i_iner(1,3)*am(3)
      omega(2)=i_iner(2,1)*am(1)+i_iner(2,2)*am(2)+i_iner(2,3)*am(3)
      omega(3)=i_iner(3,1)*am(1)+i_iner(3,2)*am(2)+i_iner(3,3)*am(3)

c$$$  
c$$$  Compute rotational kinetic energy
c$$$      

      erot=omega(1)*am(1)+omega(2)*am(2)+omega(3)*am(3)
      erot=0.5D0*erot

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
