SUBROUTINE ScalarSum(Q,Rho_k,eer)

!!$***********************************************************************
!!$   Time-stamp: <2004-11-13 16:55:32 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Apr  6 2001 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*

  USE ElecPotential_Mod
  IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*
  
  REAL(8), DIMENSION (:,:,:,:) :: Q,Rho_k
  REAL(8) :: eer
 
!!$------------------------- LOCAL VARIABLES ----------------------------*

  REAL(8)  ::  pi=3.14159265358979323846,fac,denom,eterm,eterm2,vterm,energy
  INTEGER  ::  k,k1,k2,k3,m1,m2,m3,nff,ind,jnd
  INTEGER  ::  ka3,ka2
  REAL(8)  ::  mhat1,mhat2,mhat3,msq,struc2,rkcut2,DensFact
  REAL(8)  ::  TempQ1,TempQ2
  DOUBLE COMPLEX :: TempA,TempB,TempC,TempD
  
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*
  
  rkcut2=0.25*rkcut*rkcut/pi**2
  fac = pi**2/ewald_coeff**2
  energy = 0.d0
  DO ka3=1,naz
     k3=nodez*naz+ka3
     m3 = k3 - 1
     if ( k3 .gt. nf3 )m3 = k3 - 1 - nfft3
     
     TempA=CMPLX(bsp_mod3(1,k3),bsp_mod3(2,k3))
     DO ka2=1,nay
        k2=nodey*nay+ka2
        m2 = k2 - 1
        if ( k2 .gt. nf2 )m2 = k2 - 1 - nfft2
        
        TempB=CMPLX(bsp_mod2(1,k2),bsp_mod2(2,k2))
        DO k1=1,nfft1
           
!!$c get k1,k2,k3 from the relationship
!!$c           ind = (k1-1) + (k2-1)*nfft1 + (k3-1)*nfft2*nfft1
           m1 = k1 - 1
           if ( k1 .gt. nf1 )m1 = k1 - 1 - nfft1
           
           TempC=CMPLX(bsp_mod1(1,k1),bsp_mod1(2,k1))
           TempD=CMPLX(Q(1,k1,ka2,ka3),Q(2,k1,ka2,ka3))
           
           mhat1 = recip(1,1)*m1+recip(1,2)*m2+recip(1,3)*m3
           mhat2 = recip(2,1)*m1+recip(2,2)*m2+recip(2,3)*m3
           mhat3 = recip(3,1)*m1+recip(3,2)*m2+recip(3,3)*m3
           msq = mhat1*mhat1+mhat2*mhat2+mhat3*mhat3
           
           TempQ1=DBLE(TempA*TempB*TempC*TempD)
           TempQ2=AIMAG(TempA*TempB*TempC*TempD)
           
           IF(msq .GT. rkcut2) THEN 
              eterm = 0.d0
              eterm2 = 0.d0
              DensFact=0.0D0
           ELSE IF(msq .EQ. 0.0D0) THEN
              eterm = 0.d0
              eterm2 = 0.d0
              denom=volume
              DensFact=1.0D0/denom
           ELSE
              denom = pi*msq
              eterm = dexp(-fac*msq)/denom/volume
              eterm2 = dexp(-0.5D0*fac*msq)/denom/volume
              DensFact=dexp(-0.5D0*fac*msq)/volume
              
              struc2 = TempQ1**2+TempQ2**2
              energy = energy + eterm * struc2
           END IF
           Rho_k(1,k1,ka2,ka3) = DensFact*TempQ1
           Rho_k(2,k1,ka2,ka3) = DensFact*TempQ2
           Q(1,k1,ka2,ka3) = eterm2 * TempQ1
           Q(2,k1,ka2,ka3) = eterm2 * TempQ2
         
        END DO
     END DO
  END DO
  eer = 0.5d0 * energy

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE ScalarSum
