!!$/---------------------------------------------------------------------\
!!$   Copyright  © 2006-2007 Massimo Marchi <Massimo.Marchi at cea.fr>   |
!!$                                                                      |
!!$    This software is a computer program named oracDD whose            |
!!$    purpose is to simulate and model complex molecular systems.       |
!!$    The code is written in fortran 95 compliant with Technical        |
!!$    Report TR 15581, and uses MPI-1 routines for parallel             |
!!$    coding.                                                           |
!!$                                                                      |
!!$    This software is governed by the CeCILL license under             |
!!$    French law and abiding by the rules of distribution of            |
!!$    free software.  You can  use, modify and/ or redistribute         |
!!$    the software under the terms of the CeCILL icense as              |
!!$    circulated by CEA, CNRS and INRIA at the following URL            |
!!$    "http://www.cecill.info".                                         |
!!$                                                                      |
!!$    As a counterpart to the access to the source code and rights      |
!!$    to copy, modify and redistribute granted by the license,          |
!!$    users are provided only with a limited warranty and the           |
!!$    software's author, the holder of the economic rights, and         |
!!$    the successive licensors have only limited liability.             |
!!$                                                                      |
!!$    The fact that you are presently reading this means that you       |
!!$    have had knowledge of the CeCILL license and that you accept      |
!!$    its terms.                                                        |
!!$                                                                      |
!!$    You should have received a copy of the CeCILL license along       |
!!$    with this program; if not, you can collect copies on the URL's    |
!!$    "http://www.cecill.info/licences/Licence_CeCILL_V2-en.html"       |
!!$    "http://www.cecill.info/licences/Licence_CeCILL_V2-fr.html"       |
!!$                                                                      |
!!$----------------------------------------------------------------------/

SUBROUTINE Get_Bsplines
  INTEGER :: n,i,j,k,ivect(7),count_out,count_A
  INTEGER :: ith1,ith2,ith3,k0,j0,i0,count0
  INTEGER, POINTER :: p_vect(:)=>NULL()
  REAL(8) :: w
  
  count0=0
  DO n=1,natom
     w = fr1(n)-int(fr1(n))
     CALL fill_bspline(w,order,theta1(1,n),dtheta1(1,n))
     w = fr2(n)-int(fr2(n))
     CALL fill_bspline(w,order,theta2(1,n),dtheta2(1,n))
     w = fr3(n)-int(fr3(n))
     CALL fill_bspline(w,order,theta3(1,n),dtheta3(1,n))
  END DO
  count_out=count0
END SUBROUTINE Get_Bsplines
SUBROUTINE Charges_onGrid(Cq)
  INTEGER :: nmaps,m,n,ith1,ith2,ith3,i,j,k,ka,ia,ja,k0,j0,i0
  REAL(8) :: aux,prod,chg0
  REAL(8)  :: Cq(:,:,:)
  Cq=0.0D0
  DO n=1,natom
     chg0=chg(n)
     IF(ABS(chg0) < Tol_q) CYCLE
     k0 = int(fr3(n)) - order
     DO ith3 = 1,order
        k0 = k0 + 1
        k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
        IF(k < kstart .OR. k > kend) CYCLE
        ka=k-kstart+1
        aux=chg0*theta3(ith3,n)

        j0 = int(fr2(n)) - order
        DO ith2 = 1,order
           j0 = j0 + 1
           j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2
           IF(j < jstart .OR. j > jend) CYCLE
           ja=j-jstart+1
           prod=theta2(ith2,n)*aux

           i0 = int(fr1(n)) - order
           DO ith1 = 1,order
              i0 = i0 + 1
              i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
              IF(i < istart .OR. i > iend) CYCLE
              ia=i-istart+1
              Cq(ia,ja,ka)=Cq(ia,ja,ka)+theta1(ith1,n)*prod
           END DO
        END DO
     END DO
  END DO
!!$  DO ia=1,64
!!$     DO ja=1,64
!!$        DO ka=1,64
!!$           WRITE(60,*) ia,ja,ka,Cq(ia,ja,ka)
!!$        END DO
!!$     END DO
!!$  END DO
!!$  STOP
END SUBROUTINE Charges_onGrid
SUBROUTINE ScalarSum_Transposed(Cq,nstart,nlocal,ndim1,ndim2,ndim3,nb1,nb2,nb3)

!!$************************************************************************
!!$*   Time-stamp: <04/12/15 13:02:54 marchi>                             *
!!$*                                                                      *
!!$*   New scalar_sum. Works with FFTW_TRANSPOSED_ORDER. Dimensions       *
!!$*   y and z are transposed                                             *
!!$*                                                                      *
!!$*======================================================================*
!!$*                                                                      *
!!$*              Author:  Tom Darden                                     *
!!$*              Modified by: Massimo Marchi                             *
!!$*              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$*                                                                      *
!!$*              - Tue Feb  9 1999 -                                     *
!!$*                                                                      *
!!$************************************************************************

  
  INTEGER :: nb1,nb2,nb3,ndim1,ndim2,ndim3,nstart,nlocal
  REAL(8) :: Cq(ndim1,ndim2,ndim3)

  INTEGER :: nfft1,nfft2,nfft3,naz,kwstart
  
  REAL(8) ::   fac,denom,eterm,vterm,energy,fact
  INTEGER  :: k,k1,k2,k3,m1,m2,m3,nff,ind,jnd,indtop
  INTEGER  :: nf1,nf2,nf3,ka3,ka2,i1,i2,i,ja,ka
  REAL(8) ::   mhat1,mhat2,mhat3,msq,struc2,rkcut2,ewaldcof
  LOGICAL  :: k_c1,k_c2,j_c1,j_c2,i_c1,ok_1,ok_2,ok_3
  
      
  ewaldcof=alphal
  naz=nlocal
  kwstart=nstart

  nfft1=nb1
  nfft2=nb2
  nfft3=nb3

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
     k3=kwstart-1+ka3
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
           ok_3=(.NOT. ok_1) .AND. (.NOT. j_c1)&
                &.AND. (k_c1 .AND. k_c2)
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
                 vterm = 2.d0*(fac*msq + 1.d0)/msq
                 struc2 = Cq(i1,ka2,ka3)**2 + Cq(i2,ka2,ka3)**2
                 energy = energy + eterm * struc2
                 vir(1,1) = vir(1,1) + eterm * struc2 * (vterm*mhat1&
                      &*mhat1 -1.d0)
                 vir(1,2) = vir(1,2) + eterm * struc2 * (vterm*mhat1&
                      &*mhat2)
                 vir(1,3) = vir(1,3) + eterm * struc2 * (vterm*mhat1&
                      &*mhat3)
                 vir(2,2) = vir(2,2) + eterm * struc2 * (vterm*mhat2&
                      &*mhat2 -1.d0)
                 vir(2,3) = vir(2,3) + eterm * struc2 * (vterm*mhat2&
                      &*mhat3)
                 vir(3,3) = vir(3,3) + eterm * struc2 * (vterm*mhat3&
                      &*mhat3 -1.d0)
              END IF
              Cq(i1,ka2,ka3) = eterm * Cq(i1,ka2,ka3)
              Cq(i2,ka2,ka3) = eterm * Cq(i2,ka2,ka3)
           ELSE
              IF(msq .GT. rkcut2 .OR. msq .EQ. 0.0D0) THEN 
                 eterm = 0.d0 
              ELSE
                 denom = pi*volume*bsp_mod1(k1)*bsp_mod2(k3)&
                      &*bsp_mod3(k2)*msq
                 eterm = dexp(-fac*msq)/denom
              END IF
              Cq(i1,ka2,ka3) = eterm * Cq(i1,ka2,ka3)
              Cq(i2,ka2,ka3) = eterm * Cq(i2,ka2,ka3)
           END IF
        END DO
     END DO
  END DO
  
  eer =  energy
  vir(2,1)=vir(1,2)
  vir(3,1)=vir(1,3)
  vir(3,2)=vir(2,3)
END SUBROUTINE ScalarSum_Transposed
SUBROUTINE Grad_Sum(Cq,Energy)
  REAL(8) :: Energy
  INTEGER :: nmaps,m,n,ith1,ith2,ith3,i,j,k,ka,ia,ja,k0,j0,i0
  REAL(8) :: f1,f2,f3,phid,chg0,t3,t2,dt3,dt2,term
  REAL(8) :: Cq(:,:,:)

  energy=0.0D0
  DO n=1,natom
     f1 = 0.d0
     f2 = 0.d0
     f3 = 0.d0
     chg0=chg(n)
     IF(ABS(chg0) < Tol_q) CYCLE
     phid=0.0D0
     k0 = int(fr3(n)) - order
     DO ith3 = 1,order
        k0 = k0 + 1
        k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
        IF(k < kstart .OR. k > kend) CYCLE
        t3=theta3(ith3,n)
        dt3=dtheta3(ith3,n)
        ka=k-kstart+1

        j0 = int(fr2(n)) - order
        DO ith2 = 1,order
           j0 = j0 + 1
           j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2
           IF(j < jstart .OR. j > jend) CYCLE
           t2=theta2(ith2,n)
           dt2=dtheta2(ith2,n)
           ja=j-jstart+1

           i0 = int(fr1(n)) - order
           DO ith1 = 1,order
              i0 = i0 + 1
              i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
              IF(i < istart .OR. i > iend) CYCLE
              ia=i-istart+1

              term=Cq(ia,ja,ka)
              phid=phid+term*theta1(ith1,n)*t2*t3
              term=term*chg0
!!$ force is negative of grad
              energy=energy+term*theta1(ith1,n)*t2*t3
              f1 = f1 - nfft1*term*dtheta1(ith1,n)*t2*t3
              f2 = f2 - nfft2*term*theta1(ith1,n)*dt2*t3
              f3 = f3 - nfft3*term*theta1(ith1,n)*t2*dt3
           END DO
        END DO
     END DO
     fx(n) = fx(n) + recip(1,1)*f1+recip(1,2)*f2+recip(1,3)*f3
     fy(n) = fy(n) + recip(2,1)*f1+recip(2,2)*f2+recip(2,3)*f3
     fz(n) = fz(n) + recip(3,1)*f1+recip(3,2)*f2+recip(3,3)*f3
     phi(n)=phid
  END DO
  Energy=Energy*0.5D0
END SUBROUTINE Grad_Sum
