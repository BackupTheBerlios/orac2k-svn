SUBROUTINE LoadBspModuli

!!$***********************************************************************
!!$   Time-stamp: <01/04/12 14:31:59 marchi>                           *
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

!!$------------------------- LOCAL VARIABLES ----------------------------*

  REAL(8), DIMENSION (:), ALLOCATABLE :: array,darray,bsp_arr
  INTEGER :: MaxFFT,i
  REAL(8) :: w

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

  MaxFFT=MAX(nfft1,nfft2,nfft3)
  ALLOCATE(array(order),darray(order),bsp_arr(MaxFFT))
  w = 0.d0
  call fill_bspline(w,order,array,darray)
  DO i = 1,MaxFFT
     bsp_arr(i) = 0.0D0
  END DO
  DO i = 2,order+1
     bsp_arr(i) = array(i-1)
  END DO
  CALL DftModulus(bsp_mod1,bsp_arr,nfft1)
  CALL DftModulus(bsp_mod2,bsp_arr,nfft2)
  CALL DftModulus(bsp_mod3,bsp_arr,nfft3)

  DO i=1,nfft1
     w=1.0D0/(bsp_mod1(1,i)**2+bsp_mod1(2,i)**2)
     bsp_mod1(1,i)=w*bsp_mod1(1,i)
     bsp_mod1(2,i)=-w*bsp_mod1(2,i)
  END DO
  DO i=1,nfft2
     w=1.0D0/(bsp_mod2(1,i)**2+bsp_mod2(2,i)**2)
     bsp_mod2(1,i)=w*bsp_mod2(1,i)
     bsp_mod2(2,i)=-w*bsp_mod2(2,i)
  END DO
  DO i=1,nfft3
     w=1.0D0/(bsp_mod3(1,i)**2+bsp_mod3(2,i)**2)
     bsp_mod3(1,i)=w*bsp_mod3(1,i)
     bsp_mod3(2,i)=-w*bsp_mod3(2,i)
  END DO
  DEALLOCATE(array,darray,bsp_arr)
CONTAINS
  SUBROUTINE DftModulus(bsp_mod,bsp_arr,nfft)
    implicit none
    integer nfft
    REAL(8), DIMENSION (:,:) :: bsp_mod
    REAL(8), DIMENSION (:)   :: bsp_arr
    INTEGER :: j,k
    REAL(8) :: sum1,sum2,twopi=2.0D0*3.14159265358979323846D0,arg,tiny=1.0D-7

!!$c Computes the modulus of the discrete fourier transform of bsp_arr,
!!$c  storing it into bsp_mod

    DO k = 1,nfft
       sum1 = 0.d0
       sum2 = 0.d0
       DO j = 1,nfft
          arg = twopi*(k-1)*(j-1)/nfft
          sum1 = sum1 + bsp_arr(j)*dcos(arg)
          sum2 = sum2 + bsp_arr(j)*dsin(arg)
       END DO
       bsp_mod(1,k) = sum1
       bsp_mod(2,k) = -sum2
    END DO
    DO k = 1,nfft
       IF( bsp_mod(1,k)**2+bsp_mod(2,k)**2 .lt. tiny ) THEN
          bsp_mod(1,k) = 0.5d0*(bsp_mod(1,k-1) + bsp_mod(1,k+1))
          bsp_mod(2,k) = 0.5d0*(bsp_mod(2,k-1) + bsp_mod(2,k+1))
       END if
    END DO
  END SUBROUTINE DftModulus

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE LoadBspModuli
