      SUBROUTINE grad_sum_dipole(node,nodey,nodez,nay,naz,numatoms
     &     ,charge,dipole,recip,theta1,theta2,theta3,dtheta1,dtheta2
     &     ,dtheta3,d2theta1,d2theta2,d2theta3,fpx,fpy,fpz,ene,pphi,fx
     &     ,fy,fz,fr1,fr2,fr3,order,nfft1,nfft2,nfft3,nfftdim1,nfftdim2
     &     ,nfftdim3,Q,indk1,indk2,indj1,indj2,mk,mj)

************************************************************************
*   Time-stamp: <2005-02-03 17:22:12 marchi>                           *
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
      REAL*8  fr1(numatoms),fr2(numatoms),fr3(numatoms),fpx(*),fpy(*)
     &     ,fpz(*)
      REAL*8  fx(numatoms),fy(numatoms),fz(numatoms),pphi(*),ene
      REAL*8  theta1(order,numatoms),theta2(order,numatoms),theta3(order
     &     ,numatoms),charge(numatoms),dipole(3,numatoms)
      REAL*8  dtheta1(order,numatoms),dtheta2(order,numatoms)
     &     ,dtheta3(order,numatoms),d2theta1(order,numatoms)
     &     ,d2theta2(order,numatoms),d2theta3(order,numatoms)
      REAL*8  Q(nfftdim1,nfftdim2,nfftdim3)
      INTEGER indk1(order,*),indk2(order,*),indj1(order,*),indj2(order
     &     ,*),mk(*),mj(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER n,ntot,ith1,ith2,ith3,i0,j0,k0,i,j,k
      INTEGER ka,ja,mkk,mjj,it2,it3
      REAL*8  prod,chg,t2,t3,dt2,dt3,f1,f2,f3
      REAL*8  du1_dx,du1_dy,du1_dz,du2_dx,du2_dy,du2_dz,du3_dx,du3_dy
     &     ,du3_dz
      REAL*8  mu1,mu2,mu3
      REAL*8  phi,dphi_du1,dphi_du2,dphi_du3,d2phi_du1du1,d2phi_du1du2
     &     ,d2phi_du1du3,d2phi_du2du2,d2phi_du2du3,d2phi_du3du3
      REAL*8  de_du1,de_du2,de_du3,enee,
     $       ft1,ft2,ft3,s1,s2,ft22,ft23,ft33,s3,dfx,dfy,dfz
      REAL*8  gcpu_ll,vfcp_ll,tfcp_ll,tdelta_ll,elapse,term

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      du1_dx = nfft1*recip(1,1)
      du1_dy = nfft1*recip(2,1)
      du1_dz = nfft1*recip(3,1)
      du2_dx = nfft2*recip(1,2)
      du2_dy = nfft2*recip(2,2)
      du2_dz = nfft2*recip(3,2)
      du3_dx = nfft3*recip(1,3)
      du3_dy = nfft3*recip(2,3)
      du3_dz = nfft3*recip(3,3)
      ene=0.0D0
      DO n = 1,numatoms
         mkk=mk(n)
         mjj=mj(n)
         phi=0.0D0
         dphi_du1 = 0.d0
         dphi_du2 = 0.d0
         dphi_du3 = 0.d0
         d2phi_du1du1 = 0.0D0
         d2phi_du1du2 = 0.0D0
         d2phi_du1du3 = 0.0D0
         d2phi_du2du2 = 0.0D0
         d2phi_du2du3 = 0.0D0
         d2phi_du3du3 = 0.0D0
         mu1 = nfft1*(recip(1,1)*dipole(1,n)+recip(2,1)*dipole(2,n)
     &        +recip(3,1)*dipole(3,n))
         mu2 = nfft2*(recip(1,2)*dipole(1,n)+recip(2,2)*dipole(2,n)
     &        +recip(3,2)*dipole(3,n))
         mu3 = nfft3*(recip(1,3)*dipole(1,n)+recip(2,3)*dipole(2,n)
     &        +recip(3,3)*dipole(3,n))
         chg=charge(n)
         DO it3 = 1,mkk
            ith3=indk1(it3,n)
            k=indk2(it3,n)
            ka=k-nodez*naz
            DO it2 = 1,mjj
               ith2=indj1(it2,n)
               j=indj2(it2,n)
               ja=j-nodey*nay
               i0 = int(fr1(n)) - order
               ft1 = theta2(ith2,n) * theta3(ith3,n)
               ft2 = dtheta2(ith2,n) * theta3(ith3,n)
               ft3 = theta2(ith2,n) * dtheta3(ith3,n)
               ft22 = d2theta2(ith2,n) * theta3(ith3,n)
               ft23 = dtheta2(ith2,n) * dtheta3(ith3,n)
               ft33 = theta2(ith2,n) * d2theta3(ith3,n)
               s1 = 0.d0
               s2 = 0.d0
               s3 = 0.d0
               DO ith1 = 1,order
                  i0 = i0 + 1
                  i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
                  s1 = s1 + theta1(ith1,n) * Q(i,ja,ka)
                  s2 = s2 + dtheta1(ith1,n) * Q(i,ja,ka)
                  s3 = s3 + d2theta1(ith1,n) * Q(i,j,k)
               END DO
c force is negative of grad
               phi = phi + s1 * ft1
               dphi_du1 = dphi_du1 + s2 * ft1
               dphi_du2 = dphi_du2 + s1 * ft2
               dphi_du3 = dphi_du3 + s1 * ft3
               d2phi_du1du1 = d2phi_du1du1 + s3 * ft1
               d2phi_du1du2 = d2phi_du1du2 + s2 * ft2
               d2phi_du1du3 = d2phi_du1du3 + s2 * ft3
               d2phi_du2du2 = d2phi_du2du2 + s1 * ft22
               d2phi_du2du3 = d2phi_du2du3 + s1 * ft23
               d2phi_du3du3 = d2phi_du3du3 + s1 * ft33
            END DO
         END DO
         pphi(n)=phi
         ene = ene + chg*phi + mu1*dphi_du1+mu2*dphi_du2+mu3
     &        *dphi_du3

         de_du1 = chg*dphi_du1+mu1*d2phi_du1du1+mu2*d2phi_du1du2 +
     &        mu3*d2phi_du1du3
         de_du2 = chg*dphi_du2+mu1*d2phi_du1du2+mu2*d2phi_du2du2 +
     &        mu3*d2phi_du2du3
         de_du3 = chg*dphi_du3+mu1*d2phi_du1du3+mu2*d2phi_du2du3 +
     &        mu3*d2phi_du3du3

c field is negative of grad of phi

         fx(n) = fx(n) - (dphi_du1*du1_dx+dphi_du2*du2_dx+dphi_du3
     &        *du3_dx)
         fy(n) = fy(n) - (dphi_du1*du1_dy+dphi_du2*du2_dy+dphi_du3
     &        *du3_dy)
         fz(n) = fz(n) - (dphi_du1*du1_dz+dphi_du2*du2_dz+dphi_du3
     &        *du3_dz)
         dfx = de_du1*du1_dx+de_du2*du2_dx+de_du3*du3_dx
         dfy = de_du1*du1_dy+de_du2*du2_dy+de_du3*du3_dy
         dfz = de_du1*du1_dz+de_du2*du2_dz+de_du3*du3_dz
         fpx(n)=fpx(n)-dfx
         fpy(n)=fpy(n)-dfy
         fpz(n)=fpz(n)-dfz
      END DO
      ene=ene*0.5D0
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
