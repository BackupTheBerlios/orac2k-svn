      SUBROUTINE grad_sum2(node,nodey,nodez,nay,naz,numatoms,charge
     &     ,recip,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3,fx,fy,fz
     &     ,ene,fr1,fr2,fr3,order,nfft1,nfft2,nfft3,nfftdim1,nfftdim2
     &     ,nfftdim3,Q,indk1,indk2,indj1,indj2,mk,mj)

************************************************************************
*   Time-stamp: <04/12/16 14:43:01 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Tom Darden    NIHS                             *
*              Modified by: Massimo Marchi                             *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Feb  9 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER numatoms,order,nfft1,nfft2,nfft3,node,nodey,nodez,nay,naz
      INTEGER nfftdim1,nfftdim2,nfftdim3
      REAL*8  recip(3,3)
      REAL*8  fr1(numatoms),fr2(numatoms),fr3(numatoms)
      REAL*8  ene(*)
      REAL*8  fx(numatoms),fy(numatoms),fz(numatoms)
      REAL*8  theta1(order,numatoms),theta2(order,numatoms),theta3(order
     &     ,numatoms),charge(numatoms)
      REAL*8  dtheta1(order,numatoms),dtheta2(order,numatoms)
     &     ,dtheta3(order,numatoms)
      REAL*8  Q(2,nfftdim1,nfftdim2,nfftdim3)
      INTEGER indk1(order,*),indk2(order,*),indj1(order,*),indj2(order
     &     ,*),mk(*),mj(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER n,ntot,ith1,ith2,ith3,i0,j0,k0,i,j,k
      INTEGER ka,ja,mkk,mjj,it2,it3
      REAL*8  prod,chg,t2,t3,dt2,dt3,f1,f2,f3
      REAL*8  gcpu_ll,vfcp_ll,tfcp_ll,tdelta_ll,elapse,term

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      DO n = 1,numatoms
         f1 = 0.d0
         f2 = 0.d0
         f3 = 0.d0
         mkk=mk(n)
         mjj=mj(n)
         chg=charge(n)
         DO it3 = 1,mkk
            ith3=indk1(it3,n)
            k=indk2(it3,n)
            ka=k-nodez*naz
            t3=theta3(ith3,n)
            dt3=dtheta3(ith3,n)
            DO it2 = 1,mjj
               ith2=indj1(it2,n)
               j=indj2(it2,n)
               ja=j-nodey*nay
               t2=theta2(ith2,n)
               dt2=dtheta2(ith2,n)
               i0 = int(fr1(n)) - order
               DO ith1 = 1,order
                  i0 = i0 + 1
                  i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
                  term = chg*Q(1,i,ja,ka)
                  ene(n)=ene(n)+term*theta1(ith1,n)*t2*t3
                  
c force is negative of grad
                  f1 = f1 - nfft1*term*dtheta1(ith1,n)*t2*t3
                  f2 = f2 - nfft2*term*theta1(ith1,n)*dt2*t3
                  f3 = f3 - nfft3*term*theta1(ith1,n)*t2*dt3
               END DO
            END DO
         END DO
         fx(n) = fx(n) + recip(1,1)*f1+recip(1,2)*f2+recip(1,3)*f3
         fy(n) = fy(n) + recip(2,1)*f1+recip(2,2)*f2+recip(2,3)*f3
         fz(n) = fz(n) + recip(3,1)*f1+recip(3,2)*f2+recip(3,3)*f3
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
