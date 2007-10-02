      SUBROUTINE fill_dipole_grid(node,nodey,nodez,nay,naz,numatoms
     &     ,charge,dipole,recip,dtheta1,dtheta2,dtheta3,theta1,theta2
     &     ,theta3,fr1,fr2,fr3,order,nfft1,nfft2,nfft3,nfftdim1,nfftdim2
     &     ,nfftdim3,Q,indk1,indk2,indj1,indj2,mk,mj)

************************************************************************
*   Time-stamp: <2005-02-03 17:19:29 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Tom Darden   NIHS                              *
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

      INTEGER node,nodey,nodez,nay,naz,numatoms,order,nfft1,nfft2,nfft3
      INTEGER nfftdim1,nfftdim2,nfftdim3
      REAL*8  fr1(numatoms),fr2(numatoms),fr3(numatoms),recip(3,3)
      REAL*8  theta1(order,numatoms),theta2(order,numatoms),theta3(order
     &     ,numatoms),dtheta1(order,numatoms),dtheta2(order,numatoms)
     &     ,dtheta3(order,numatoms),charge(numatoms),dipole(3,numatoms)
      REAL*8  Q(nfftdim1,nfftdim2,nfftdim3)
      INTEGER indk1(order,*),indk2(order,*),indj1(order,*),indj2(order
     &     ,*),mk(*),mj(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER n,ntot,ith1,ith2,ith3,i0,j0,k0,i,j,k
      INTEGER ka,ja,mkk,mjj,it2,it3
      REAL*8  prod1,prod2,chg,aux,mu1,mu2,mu3
      REAL*8  gcpu_ll,vfcp_ll,tfcp_ll,tdelta_ll,elapse

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      ntot = nfftdim1*nfftdim2*nfftdim3
      call clearQ(Q,ntot)

      DO n = 1,numatoms
         mkk=mk(n)
         mjj=mj(n)
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
               prod1 = theta2(ith2,n)*theta3(ith3,n)*chg
     &              +dtheta2(ith2,n)*theta3(ith3,n)*mu2 +theta2(ith2,n)
     &              *dtheta3(ith3,n)*mu3
               prod2 = theta2(ith2,n)*theta3(ith3,n)*mu1
               i0 = int(fr1(n)) - order
               DO ith1 = 1,order
                  i0 = i0 + 1
                  i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
                  Q(i,ja,ka) = Q(i,ja,ka) + theta1(ith1,n)*prod1 +
     &                 dtheta1(ith1,n)*prod2
               END DO
            END DO
         END DO
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
