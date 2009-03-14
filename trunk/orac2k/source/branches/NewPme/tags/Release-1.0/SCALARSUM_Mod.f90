MODULE SCALARSUM_Mod

!!$***********************************************************************
!!$   Time-stamp: <2005-10-25 23:13:33 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Tue Oct 25 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*

CONTAINS
  SUBROUTINE Normal(phi,Ex,Ey,Ez,node,kstart,naz,ewaldcof&
       &,volume,recip,bsp_mod1,bsp_mod2,bsp_mod3,nfft1,nfft2,nfft3&
       &,nfftdim1,nfftdim2,nfftdim3,eer,rkcut)


!!$---- This subroutine is part of the program ORAC ----*

!!$======================== DECLARATIONS ================================*

      IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

      INTEGER :: nfft1,nfft2,nfft3,node,naz
      REAL(8), DIMENSION(nfftdim1,nfftdim2,nfftdim3) ::  phi,Ex,Ey,Ez
      REAL(8) ::  bsp_mod1(:),bsp_mod2(:),bsp_mod3(:)&
           &,ewaldcof,volume,rkcut 
      REAL(8) :: eer,vir(3,3)
      REAL(8) ::  recip(3,3)

!!$------------------------- LOCAL VARIABLES ----------------------------*

      REAL(8) :: fac,denom,eterm,vterm,energy,fact
      INTEGER :: k,k1,k2,k3,m1,m2,m3,nff,ind,jnd,indtop
      INTEGER :: nf1,nf2,nf3,ka3,ka2,i1,i2,i,ja,ka,kstart
      REAL(8) :: mhat1,mhat2,mhat3,msq,struc2,rkcut2
      LOGICAL :: k_c1,k_c2,j_c1,j_c2,i_c1,ok_1,ok_2,ok_3
      REAL(8), SAVE :: pi=3.14159265358979323846D0

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*
      
      rkcut2=0.25*rkcut*rkcut/pi**2
      fac = pi**2/alphal**2

      nf1 = nfft1/2
      if ( 2*nf1 .lt. nfft1 )nf1 = nf1+1
      nf2 = nfft2/2
      if ( 2*nf2 .lt. nfft2 )nf2 = nf2+1
      nf3 = nfft3/2
      if ( 2*nf3 .lt. nfft3 )nf3 = nf3+1

      energy = 0.d0
      fact=2.0D0

      DO ka3=1,naz
         k3=kstart-1+ka3
         m3 = k3 - 1
         IF ( k3 .GT. nf3 )m3 = k3 - 1 - nfft3
         k_c1=k3 .GE. 2 
         k_c2=k3 .LE. nfft3/2+1
         DO ka2=1,nfft2
            k2=ka2
            m2 = k2 - 1
            j_c1=k2 .GE. 2 
            j_c2=k2 .LE. nfft2/2+1
            
            IF ( k2 .GT. nf2 )m2 = k2 - 1 - nfft2
            DO k1=1,nfft1/2+1
               i1=(k1-1)*2+1
               i2=(k1-1)*2+2
               i_c1=k1 .GT. 1
               
               ok_1=i_c1
               ok_2=(.NOT. ok_1) .AND. (j_c1 .AND. j_c2)
               ok_3=(.NOT. ok_1) .AND. (.NOT. j_c1) .AND. (k_c1 .AND.&
                    & k_c2) 
               m1 = k1 - 1
               IF( k1 .gt. nf1 )m1 = k1 - 1 - nfft1
               mhat1 = oca(1,1)*m1+oca(1,2)*m2+oca(1,3)*m3
               mhat2 = oca(2,1)*m1+oca(2,2)*m2+oca(2,3)*m3
               mhat3 = oca(3,1)*m1+oca(3,2)*m2+oca(3,3)*m3
               msq = mhat1*mhat1+mhat2*mhat2+mhat3*mhat3
               IF(ok_1 .OR. ok_2 .OR. ok_3) THEN
                  IF(msq .GT. rkcut2 .OR. msq .EQ. 0.0D0) THEN 
                     eterm = 0.d0 
                  ELSE
                     denom = pi*volume*bsp_mod1(k1)*bsp_mod2(k2)&
                          &*bsp_mod3(k3)*msq 
                     eterm = dexp(-fac*msq)/denom
                     struc2 = phi(i1,ka2,ka3)**2 + phi(i2,ka2,ka3)**2
                     energy = energy + eterm * struc2
                  END IF
               ELSE
                  IF(msq .GT. rkcut2 .OR. msq .EQ. 0.0D0) THEN 
                     eterm = 0.d0 
                  ELSE
                     denom = pi*volume*bsp_mod1(k1)*bsp_mod2(k2)&
                          &*bsp_mod3(k3)*msq 
                     eterm = dexp(-fac*msq)/denom
                  END IF
               END IF
               phi(i1,ka2,ka3) = eterm * phi(i1,ka2,ka3) 
               phi(i2,ka2,ka3) = eterm * phi(i2,ka2,ka3) 
               Ex (i2,ka2,ka3) = mhat1 * phi(i1,ka2,ka3)
               Ex (i1,ka2,ka3) = - mhat1 * phi(i2,ka2,ka3)
               Ey (i2,ka2,ka3) = mhat2 * phi(i1,ka2,ka3)
               Ey (i1,ka2,ka3) = - mhat2 * phi(i2,ka2,ka3)
               Ez (i2,ka2,ka3) = mhat3 * phi(i1,ka2,ka3)
               Ez (i1,ka2,ka3) = - mhat3 * phi(i2,ka2,ka3)
            END DO
         END DO
      END DO
      eer =  energy

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

    END SUBROUTINE Normal
    SUBROUTINE Transpose(phi,Ex,Ey,Ez,node,kstart,naz,ewaldcof,volume&
         &,recip,bsp_mod1,bsp_mod2,bsp_mod3,nfft1,nfft2,nfft3&
         &,nfftdim1,nfftdim2,nfftdim3,eer,rkcut)

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*

      IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

      INTEGER :: nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,node,naz
      REAL(8), DIMENSION(nfftdim1,nfftdim2,nfftdim3) :: phi, Ex, Ey, Ez
      REAL(8)  bsp_mod1(nfft1),bsp_mod2(nfft2),bsp_mod3(nfft3)&
           &,ewaldcof,volume,rkcut
      REAL(8) :: eer
      REAL(8) :: recip(3,3)

!!$------------------------- LOCAL VARIABLES ----------------------------*

      REAL(8)  :: fac,denom,eterm,vterm,energy,fact
      INTEGER  :: k,k1,k2,k3,m1,m2,m3,nff,ind,jnd,indtop
      INTEGER  :: nf1,nf2,nf3,ka3,ka2,i1,i2,i,ja,ka,kstart
      REAL(8)  :: mhat1,mhat2,mhat3,msq,struc2,rkcut2
      LOGICAL  :: k_c1,k_c2,j_c1,j_c2,i_c1,ok_1,ok_2,ok_3
      REAL(8), SAVE :: pi=3.14159265358979323846D0
      
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

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
               ok_3=(.NOT. ok_1) .AND. (.NOT. j_c1) .AND. (k_c1 .AND. k_c2)
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
                     denom = pi*volume*bsp_mod1(k1)*bsp_mod2(k3)&
                          &*bsp_mod3(k2)*msq 
                     eterm = dexp(-fac*msq)/denom
                     struc2 = phi(i1,ka2,ka3)**2 + phi(i2,ka2,ka3)**2
                     energy = energy + eterm * struc2
                  END IF
               ELSE
                  IF(msq .GT. rkcut2 .OR. msq .EQ. 0.0D0) THEN 
                     eterm = 0.d0 
                  ELSE
                     denom = pi*volume*bsp_mod1(k1)*bsp_mod2(k3)&
                          &*bsp_mod3(k2)*msq
                     eterm = dexp(-fac*msq)/denom
                  END IF
               END IF
               phi(i1,ka2,ka3) = eterm * phi(i1,ka2,ka3)
               phi(i2,ka2,ka3) = eterm * phi(i2,ka2,ka3)
               Ex (i2,ka2,ka3) = mhat1 * phi(i1,ka2,ka3)
               Ex (i1,ka2,ka3) = - mhat1 * phi(i2,ka2,ka3)
               Ey (i2,ka2,ka3) = mhat2 * phi(i1,ka2,ka3)
               Ey (i1,ka2,ka3) = - mhat2 * phi(i2,ka2,ka3)
               Ez (i2,ka2,ka3) = mhat3 * phi(i1,ka2,ka3)
               Ez (i1,ka2,ka3) = - mhat3 * phi(i2,ka2,ka3)
            END DO
         END DO
      END DO

      eer =  energy

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

    END SUBROUTINE Transpose

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE SCALARSUM_Mod
