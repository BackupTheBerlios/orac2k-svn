      SUBROUTINE scalar_sum_transp(node,kstart,naz,Q,ewaldcof,volume
     &     ,recip,bsp_mod1,bsp_mod2,bsp_mod3,nfft1,nfft2,nfft3,nfftdim1
     &     ,nfftdim2,nfftdim3,eer,vir,rkcut)

************************************************************************
*   Time-stamp: <2009-03-13 16:33:11 marchi>                             *
*                                                                      *
*   New scalar_sum. Works with FFTW_TRANSPOSED_ORDER. Dimensions       *
*   y and z are transposed                                             *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Tom Darden                                     *
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

      INTEGER nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,node,naz
      REAL*8  Q(nfftdim1,nfftdim2,nfftdim3)
      REAL*8  bsp_mod1(*),bsp_mod2(*),bsp_mod3(*)
     &     ,ewaldcof,volume,rkcut
      REAL*8  eer,vir(3,3)
      REAL*8  recip(3,3)

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8  pi,fac,denom,eterm,vterm,energy,fact
      INTEGER  k,k1,k2,k3,m1,m2,m3,nff,ind,jnd,indtop
      INTEGER  nf1,nf2,nf3,ka3,ka2,i1,i2,i,ja,ka,kstart
      REAL*8  mhat1,mhat2,mhat3,msq,struc2,rkcut2
      LOGICAL k_c1,k_c2,j_c1,j_c2,i_c1,ok_1,ok_2,ok_3
      
*----------------------- EXECUTABLE STATEMENTS ------------------------*

      pi = 3.14159265358979323846
      rkcut2=0.25*rkcut*rkcut/pi**2
      fac = pi**2/ewaldcof**2
      nff = nfft1*nfft2
      nf1 = nfft1/2
      if ( 2*nf1 .lt. nfft1 )nf1 = nf1+1
      nf2 = nfft2/2
      if ( 2*nf2 .lt. nfft2 )nf2 = nf2+1
      nf3 = nfft3/2
      if ( 2*nf3 .lt. nfft3 )nf3 = nf3+1
      energy = 0.d0
      fact=2.0D0
      DO k1 = 1,3
         DO k2 = 1,3
            vir(k1,k2) = 0.0D0
         END DO
      END DO

      DO ka3=1,naz
         k3=kstart-1+ka3
         m3 = k3 - 1
         if ( k3 .gt. nf3 )m3 = k3 - 1 - nfft3
         k_c1=k3 .GE. 2 
         k_c2=k3 .LE. nfft3/2+1
         DO ka2=1,nfft2
            k2=ka2
            m2 = k2 - 1
            j_c1=k2 .GE. 2 
            j_c2=k2 .LE. nfft2/2+1

            if ( k2 .gt. nf2 )m2 = k2 - 1 - nfft2
            DO k1=1,nfft1/2+1
               i1=(k1-1)*2+1
               i2=(k1-1)*2+2
               i_c1=k1 .GT. 1

               ok_1=i_c1
               ok_2=(.NOT. ok_1) .AND. (j_c1 .AND. j_c2)
               ok_3=(.NOT. ok_1) .AND. (.NOT. j_c1)
     &              .AND. (k_c1 .AND. k_c2)
               m1 = k1 - 1
               if ( k1 .gt. nf1 )m1 = k1 - 1 - nfft1
               mhat1 = recip(1,1)*m1+recip(1,2)*m3+recip(1,3)*m2
               mhat2 = recip(2,1)*m1+recip(2,2)*m3+recip(2,3)*m2
               mhat3 = recip(3,1)*m1+recip(3,2)*m3+recip(3,3)*m2
               msq = mhat1*mhat1+mhat2*mhat2+mhat3*mhat3
               IF(ok_1 .OR. ok_2 .OR. ok_3) THEN
                  IF(msq .GT. rkcut2 .OR. msq .EQ. 0.0D0) THEN 
                     eterm = 0.d0 
                  ELSE
                     denom = pi*volume*bsp_mod1(k1)*bsp_mod2(k3)
     &                    *bsp_mod3(k2)*msq
                     eterm = dexp(-fac*msq)/denom
                     vterm = 2.d0*(fac*msq + 1.d0)/msq
                     struc2 = Q(i1,ka2,ka3)**2 + Q(i2,ka2,ka3)**2
                     energy = energy + eterm * struc2
                     vir(1,1) = vir(1,1) + eterm * struc2 * (vterm
     &                    *mhat1*mhat1 -1.d0)
                     vir(1,2) = vir(1,2) + eterm * struc2 * (vterm
     &                    *mhat1*mhat2)
                     vir(1,3) = vir(1,3) + eterm * struc2 * (vterm
     &                    *mhat1*mhat3)
                     vir(2,2) = vir(2,2) + eterm * struc2 * (vterm
     &                    *mhat2*mhat2 -1.d0)
                     vir(2,3) = vir(2,3) + eterm * struc2 * (vterm
     &                    *mhat2*mhat3)
                     vir(3,3) = vir(3,3) + eterm * struc2 * (vterm
     &                    *mhat3*mhat3 -1.d0)
                  END IF
                  Q(i1,ka2,ka3) = eterm * Q(i1,ka2,ka3)
                  Q(i2,ka2,ka3) = eterm * Q(i2,ka2,ka3)
               ELSE
                  IF(msq .GT. rkcut2 .OR. msq .EQ. 0.0D0) THEN 
                     eterm = 0.d0 
                  ELSE
                     denom = pi*volume*bsp_mod1(k1)*bsp_mod2(k3)
     &                    *bsp_mod3(k2)*msq
                     eterm = dexp(-fac*msq)/denom
                  END IF
                  Q(i1,ka2,ka3) = eterm * Q(i1,ka2,ka3)
                  Q(i2,ka2,ka3) = eterm * Q(i2,ka2,ka3)
               END IF
            END DO
         END DO
      END DO

      eer =  energy
      vir(2,1)=vir(1,2)
      vir(3,1)=vir(1,3)
      vir(3,2)=vir(2,3)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*
      RETURN
      END
