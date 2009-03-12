      SUBROUTINE fill_charge_grid(node,nodey,nodez,jstart,kstart
     &     ,numatoms,charge,theta1,theta2,theta3,fr1,fr2,fr3,order,nfft1
     &     ,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,Q,indk1,indk2,indj1
     &     ,indj2,mk,mj)

************************************************************************
*   Time-stamp: <2009-03-12 12:27:59 marchi>                             *
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

      INTEGER node,nodey,nodez,jstart,kstart,numatoms,order,nfft1,nfft2
     &     ,nfft3
      INTEGER nfftdim1,nfftdim2,nfftdim3
      REAL*8  fr1(numatoms),fr2(numatoms),fr3(numatoms)
      REAL*8  theta1(order,numatoms),theta2(order,numatoms),theta3(order
     &     ,numatoms),charge(numatoms)
      REAL*8  Q(nfftdim1,nfftdim2,nfftdim3)
      INTEGER indk1(order,*),indk2(order,*),indj1(order,*),indj2(order
     &     ,*),mk(*),mj(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER n,ntot,ith1,ith2,ith3,i0,j0,k0,i,j,k
      INTEGER ka,ja,mkk,mjj,it2,it3
      REAL*8  prod,chg,aux
      REAL*8  gcpu_ll,vfcp_ll,tfcp_ll,tdelta_ll,elapse

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      ntot = nfftdim1*nfftdim2*nfftdim3
      call clearQ(Q,ntot)

      DO n = 1,numatoms
         mkk=mk(n)
         mjj=mj(n)
         chg=charge(n)
         DO it3 = 1,mkk
            ith3=indk1(it3,n)
            k=indk2(it3,n)
            ka=k-kstart+1
            aux=chg*theta3(ith3,n)
            DO it2 = 1,mjj
               ith2=indj1(it2,n)
               j=indj2(it2,n)
               ja=j-jstart+1
               prod = theta2(ith2,n)*aux
               i0 = int(fr1(n)) - order
               DO ith1 = 1,order
                  i0 = i0 + 1
                  i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
                  Q(i,ja,ka) = Q(i,ja,ka) + theta1(ith1,n)
     &                 * prod
               END DO
            END DO
         END DO
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
