      subroutine pmesh_kspace_get_sizes(
     $     nfft1,nfft2,nfft3,numatoms,order,
     $     sizfftab,sizffwrk,siztheta,siz_Q,sizheap,sizstack)
      implicit none
      integer nfft1,nfft2,nfft3,numatoms,order,
     $     sizfftab,sizffwrk,siztheta,siz_Q,sizheap,sizstack

c INPUT  
c      nfft1,nfft2,nfft3,numatoms,order
c      nfft1,nfft2,nfft3 are the dimensions of the charge grid array
c      numatoms is number of atoms
c      order is the order of B-spline interpolation

c OUTPUT
c      sizfftab,sizffwrk,siztheta,siz_Q
c      sizfftab is permanent 3d fft table storage
c      sizffwrk is temporary 3d fft work storage
c      siztheta is size of arrays theta1-3 dtheta1-3
c      sizheap is total size of permanent storage
c      sizstack is total size of temporary storage


c This routine computes the above output parameters needed for 
c heap or stack allocation.

      integer nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork

      call get_fftdims(nfft1,nfft2,nfft3,
     $       nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,
     $       sizfftab,sizffwrk)
      siztheta = numatoms*order
      siz_Q = 2*nfftdim1*nfftdim2*nfftdim3
      sizheap = nfft1+nfft2+nfft3+sizfftab
      sizstack = siz_Q+6*siztheta+sizffwrk+3*numatoms
C      write(6,*)'total HEAP storage needed = ',sizheap
C      write(6,*)'total STACK storage needed = ',sizstack
      return
      end
c----------------------------------------------------
      subroutine pmesh_kspace_setup(
     $    bsp_mod1,bsp_mod2,bsp_mod3,fftable,ffwork,
     $    nfft1,nfft2,nfft3,order,sizfftab,sizffwrk)
      implicit none

c  see DO_PMESH_KSPACE for explanation of arguments

      integer nfft1,nfft2,nfft3,order,sizfftab,sizffwrk
      double precision bsp_mod1(nfft1),bsp_mod2(nfft2),
     +   bsp_mod3(nfft3)
      double precision fftable(sizfftab),ffwork(sizffwrk)
   
      double precision dummy
      integer nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw

      call get_fftdims(nfft1,nfft2,nfft3,
     $       nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw)
      call load_bsp_moduli(bsp_mod1,bsp_mod2,bsp_mod3,
     $   nfft1,nfft2,nfft3,order)
      call fft_setup(dummy,fftable,ffwork,
     $      nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,
     $      nfftable,nffwork)
      return
      end
c----------------------------------------------------------------------
      subroutine get_scaled_fractionals(nato,x,y,z,recip,nfft1,nfft2
     &     ,nfft3,fr1,fr2,fr3)
      implicit none

c INPUT:
c      numatoms: number of atoms
c      x,y,z: arrays of cartesian coords
c      recip: the 3x3 array of reciprocal vectors stored as columns
c OUTPUT:
c     fr1,fr2,fr3 the scaled and shifted fractional coords

      integer nato,nfft1,nfft2,nfft3
      double precision x(*),y(*),z(*),recip(3,3)
      double precision fr1(*),fr2(*),fr3(*)

      integer n
      double precision w1,w2,w3

      do 100 n = 1,nato
        w1 = x(n)*recip(1,1)+y(n)*recip(2,1)+z(n)*recip(3,1)
        w2 = x(n)*recip(1,2)+y(n)*recip(2,2)+z(n)*recip(3,2)
        w3 = x(n)*recip(1,3)+y(n)*recip(2,3)+z(n)*recip(3,3)
        fr1(n) = nfft1*(w1 - anint(w1) + 0.5d0)
        fr2(n) = nfft2*(w2 - anint(w2) + 0.5d0)
        fr3(n) = nfft3*(w3 - anint(w3) + 0.5d0)
100   continue
      return
      end
c---------------------------------------------------------------
      subroutine load_bsp_moduli(bsp_mod1,bsp_mod2,bsp_mod3,
     $   nfft1,nfft2,nfft3,order)
      implicit none
      integer nfft1,nfft2,nfft3,order
      double precision bsp_mod1(nfft1),bsp_mod2(nfft2),
     +   bsp_mod3(nfft3)
      integer MAXORDER
      parameter (MAXORDER=25)
      integer MAXNFFT
      parameter (MAXNFFT=1000)
      double precision array(MAXORDER),darray(MAXORDER),w
      double precision bsp_arr(MAXNFFT)
      integer i,maxn,opt_infl
c this routine loads the moduli of the inverse DFT of the B splines
c bsp_mod1-3 hold these values, nfft1-3 are the grid dimensions,
c Order is the order of the B spline approx.
      opt_infl=1
      if ( order .gt. MAXORDER )then
       write(6,*)'order too large! check on MAXORDER'
       stop
      endif
      maxn = max(nfft1,nfft2,nfft3)
      if ( maxn .gt. MAXNFFT )then 
       write(6,*)'nfft1-3 too large! check on MAXNFFT'
       stop
      endif
      w = 0.d0
      call fill_bspline(w,order,array,darray)
      do 100 i = 1,maxn
        bsp_arr(i) = 0.d0
100   continue
      do 150 i = 2,order+1
       bsp_arr(i) = array(i-1)
150   continue
      call DFTMOD(bsp_mod1,bsp_arr,nfft1)
      call factor_lambda(bsp_mod1,nfft1,order,bsp_mod1,opt_infl)
      call DFTMOD(bsp_mod2,bsp_arr,nfft2)
      call factor_lambda(bsp_mod2,nfft2,order,bsp_mod2,opt_infl)
      call DFTMOD(bsp_mod3,bsp_arr,nfft3)
      call factor_lambda(bsp_mod3,nfft3,order,bsp_mod3,opt_infl)
      return
      end
      subroutine factor_lambda(bsp_mod,nfft,order,prefac,opt_infl)
      implicit none
      integer nfft,order,opt_infl
      REAL*8 bsp_mod(nfft),prefac(nfft)
      integer KCUT
      parameter (KCUT=50)
C factor in optimal lambda coefficient for bspline approximant
c of complex exponentials, thus modify influence function
      integer k,nf,m,order2
      REAL*8 lambda,ratio,x,pi,gsum,gsum2
      nf = nfft / 2
      pi= 3.14159265358979323846
      order2 = 2*order
c something wrong with get_tarray for large order.
c use clunky but reliable gamma_sum
      do k = 1,nfft
        if ( opt_infl .eq. 0 )then
          lambda = 1.d0
        else
         m = k-1
         if ( k .gt. nf )m = k - 1 - nfft
         x = pi*DBLE(m)/nfft
         order2 = 2*order
         if ( m .eq. 0)then
           lambda = 1.d0
         else
          call gamma_sum(gsum,m,nfft,order,KCUT)
          call gamma_sum(gsum2,m,nfft,order2,KCUT)
          lambda = gsum/gsum2
         endif
        endif
        prefac(k) = bsp_mod(k)/lambda**2
      enddo
      return
      end
      subroutine gamma_sum(gsum,m,nfft,order,KCUT)
      implicit none
      double precision gsum
      integer m,nfft,order,KCUT
      double precision frac,term,x,pi
      integer k
      pi= 3.14159265358979323846
      if ( m .eq. 0 )then
        gsum = 1.d0
        return
      endif
      frac = DBLE(m)/nfft 
      x = pi*frac
      gsum = 1.d0
      do k = 1,kcut
        gsum = gsum + (x/(x + pi*k))**order
      enddo 
      do k = 1,kcut
        gsum = gsum + (x/(x - pi*k))**order
      enddo 
      return
      end

c------------------------------------------------------------------------
      subroutine DFTMOD(bsp_mod,bsp_arr,nfft)
      implicit none
      integer nfft
      double precision bsp_mod(nfft),bsp_arr(nfft)
c Computes the modulus of the discrete fourier transform of bsp_arr,
c  storing it into bsp_mod

      integer j,k
      double precision sum1,sum2,twopi,arg,tiny
      twopi = 2.d0*3.14159265358979323846
      tiny = 1.d-7
      do 300 k = 1,nfft
       sum1 = 0.d0
       sum2 = 0.d0
       do 250 j = 1,nfft
         arg = twopi*(k-1)*(j-1)/nfft
         sum1 = sum1 + bsp_arr(j)*dcos(arg)
         sum2 = sum2 + bsp_arr(j)*dsin(arg)
250    continue
       bsp_mod(k) = sum1**2 + sum2**2
300   continue
      do 400 k = 1,nfft
       if ( bsp_mod(k) .lt. tiny )
     $     bsp_mod(k) = 0.5d0*(bsp_mod(k-1) + bsp_mod(k+1))
400   continue
      return
      end
c-----------------------------------------------------------
      subroutine clearQ(Q,ntot)
      integer ntot
      double precision Q(ntot)
      integer i
      do 10 i = 1,ntot
        Q(i) = 0.d0
10    continue
      return
      end
c-------------------------------------------------------------
      subroutine check_virial(self_ene,adj_ene,dir_ene,rec_ene,
     $       adj_vir,rec_vir,dir_vir)
      implicit none
      double precision self_ene,adj_ene,dir_ene,rec_ene
      double precision adj_vir(6),rec_vir(6),dir_vir(6)

      double precision etot,svir,relerr
      etot = self_ene+adj_ene+dir_ene+rec_ene
      svir = adj_vir(1)+rec_vir(1)+dir_vir(1)+
     $       adj_vir(4)+rec_vir(4)+dir_vir(4)+
     $       adj_vir(6)+rec_vir(6)+dir_vir(6)
      relerr = 2.d0*abs(etot+svir)/(abs(etot)+abs(svir))
      write(6,*)'tot ene =   ',etot
      write(6,*)'trace vir = ',svir
      write(6,*)'rel error = ',relerr
      return
      end
c-------------------------------------------------------------
      subroutine check_force(numatoms,
     $      fx1,fy1,fz1,fx2,fy2,fz2,fdx,fdy,fdz)
      integer numatoms
      double precision fx1(*),fy1(*),fz1(*),fx2(*),fy2(*),fz2(*),
     $      fdx(*),fdy(*),fdz(*)

      double precision rms_num,rms_den,rms
      integer i
      rms_num = 0.d0
      rms_den = 0.d0
      do 100 i = 1,numatoms
       rms_num = rms_num + (fx1(i)-fx2(i))**2 + (fy1(i)-fy2(i))**2 +
     $          (fz1(i)-fz2(i))**2
       rms_den = rms_den + fdx(i)**2 + fdy(i)**2 + fdz(i)**2
100   continue
      rms = dsqrt(rms_num/rms_den)
      write(6,*)'rms force err = ',rms
      return
      end
c-------------------------------------------------------------

      subroutine pubz3di(n1,n2,n3,table,ntable)
      implicit none
      integer n1,n2,n3,ntable
      double precision table(ntable,3)
c ntable should be 4*max(n1,n2,n3) +15


      call cffti(n1,table(1,1))
      call cffti(n2,table(1,2))
      call cffti(n3,table(1,3))

      return
      end
*****************************************************************************
      subroutine pubz3d(isign,n1,n2,n3,w,ld1,ld2,table,ntable,
     $    work,nwork)
      implicit none

      integer n1,n2,n3,ld1,ld2,isign,ntable,nwork
      double complex w(ld1,ld2,n3)
      double complex work( nwork)
      double precision table(ntable,3)

      integer i,j,k
c ntable should be 4*max(n1,n2,n3) +15
c nwork should be max(n1,n2,n3)
c
c   transform along X  first ...
c
      do 100 k = 1, n3
       do 90 j = 1, n2
        do 70 i = 1,n1
          work(i) = w(i,j,k)
70      continue
        if ( isign .eq. -1) call cfftf(n1,work,table(1,1))
        if ( isign .eq. 1) call cfftb(n1,work,table(1,1))
        do 80 i = 1,n1
          w(i,j,k) = work(i)
80      continue
90     continue
100   continue
c
c   transform along Y then ...
c
      do 200 k = 1,n3
       do 190 i = 1,n1
        do 170 j = 1,n2
          work(j) = w(i,j,k)
170     continue
        if ( isign .eq. -1) call cfftf(n2,work,table(1,2))
        if ( isign .eq. 1) call cfftb(n2,work,table(1,2))
        do 180 j = 1,n2
          w(i,j,k) = work(j)
180     continue
190    continue
200   continue
c
c   transform along Z finally ...
c
      do 300 i = 1, n1
       do 290 j = 1, n2
        do 270 k = 1,n3
          work(k) = w(i,j,k)
270     continue
        if ( isign .eq. -1) call cfftf(n3,work,table(1,3))
        if ( isign .eq. 1) call cfftb(n3,work,table(1,3))
        do 280 k = 1,n3
          w(i,j,k) = work(k)
280     continue
290    continue
300   continue

      return
      end
c----------------------------------------------------
      subroutine fill_bspline2(w,order,array,darray,d2array)
c---------- use standard B-spline recursions: see doc file
      implicit none
      integer order
      double precision w,array(order),darray(order),d2array(order)
      integer k
c do linear case
      call init_bspline(array,w,order)
c compute standard b-spline recursion
      do k = 3,order-2
       call one_pass_bspline(array,w,k)
      enddo
c perform standard b-spline 2nd differentiation
      call diff_bspline(array,darray,order)
      call diff_bspline(darray,d2array,order)
c perform standard b-spline differentiation
      call one_pass_bspline(array,w,order-1)
      call diff_bspline(array,darray,order)
c one more recursion
      call one_pass_bspline(array,w,order)

      return
      end
      subroutine diff_bspline(c,d,n)
      implicit none
      REAL*8 c(*),d(*)
      integer n
c Using notation from Essmann et al; w = u-[u] and
c array(j) = M_n(w + order - j)  where n is order
c DERIVATIVE:    d/dw M_n(w) = M_n-1(w) - M_n-1(w-1)
c i.e.   d/dw M_n(w+n-j) = M_n-1(w+n-j) - M_n-1(w+n-j-1)
c i.e.   new(j) = old(j-1) - old(j)
c where old is array before one_pass (thus n->n-1) and new is array afterwards

      integer j
      d(1) = -c(1)
      do j = 2,n
       d(j) = c(j-1) - c(j)
      enddo
      return
      end
c---------------------------------------------------------------------
      subroutine one_pass_bspline(c,w,n)
      implicit none
      REAL*8 c(*),w
      integer n
c Using notation from Essmann et al; w = u-[u] and
c array(j) = M_n(w + order - j)  where n is order
c RECURSION:  M_n(w) = (w/(n-1))*M_n-1(w)+((n-w)/(n-1))*M_n-1(w-1)
c i.e.   M_n(w+n-j) = ((w+n-j)/(n-1))*M_n-1(w+n-j)+((j-w)/(n-1))*M_n-1(w+n-j-1)
c i.e.   new(j) = ((w+n-j)/(n-1))*old(j-1) + ((j-w)/(n-1))*old(j)
c where old is array before one_pass (thus n->n-1) and new is array afterwards
c write backwards to do it with one array

      REAL*8 div
      integer j

      div = 1.d0 / (n-1)
      c(n) = div*w*c(n-1)
      do j = 1,n-2
       c(n-j) = div*((w+j)*c(n-j-1) + (n-j-w)*c(n-j))
      enddo
      c(1) = div*(1-w)*c(1)
      return
      end
c---------------------------------------------------------------------
      subroutine init_bspline(c,w,n)
      implicit none
c this gives back the spline array of order 2 (linear)
      integer n,i
      REAL*8 c(n),w
      c(2) = w
      c(1) = 1.d0 - w
      do i = 3,n
       c(i) = 0.d0
      enddo
      return
      end

      subroutine fill_bspline(w,order,array,darray)
c---------- use standard B-spline recursions: see doc file
      implicit none
      integer order
      double precision w,array(order),darray(order)

      integer k
c do linear case
      call init(array,w,order)
c compute standard b-spline recursion
      do 10 k = 3,order-1
       call one_pass(array,w,k)
10    continue
c perform standard b-spline differentiation
      call diff(array,darray,order)
c one more recursion
      call one_pass(array,w,order)
      return
      end
c---------------------------------------------------
      subroutine init(c,x,order)
      implicit none
      integer order
      double precision c(order),x
      c(order) = 0.d0
      c(2) = x
      c(1) = 1.d0 - x
      return
      end
c-------------------------------------
      subroutine one_pass(c,x,k)
      implicit none
      double precision c(*),x
      integer k

      double precision div
      integer j

      div = 1.d0 / (k-1)
      c(k) = div*x*c(k-1)
      do 100 j = 1,k-2
       c(k-j) = div*((x+j)*c(k-j-1) + (k-j-x)*c(k-j))
100   continue
      c(1) = div*(1-x)*c(1)
      return
      end
c-------------------------------------
      subroutine diff(c,d,order)
      implicit none
      double precision c(*),d(*)
      integer order

      integer j
      d(1) = -c(1)
      do 10 j = 2,order
       d(j) = c(j-1) - c(j)
10    continue
      return
      end
c-------------------------------------


C   FFT CALLS

      subroutine get_fftdims(nfft1,nfft2,nfft3,
     $       nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,
     $       sizfftab,sizffwrk)
      implicit none
      integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,
     $       nfftable,nffwork,sizfftab,sizffwrk
      integer n,nfftmax

      nfftmax = max(nfft1,nfft2,nfft3)
      nfftdim1 = nfft1
      n = nfft1/2
      if ( nfft1 .eq. 2*n )nfftdim1 = nfft1+1
      nfftdim2 = nfft2
      n = nfft2/2
      if ( nfft2 .eq. 2*n )nfftdim2 = nfft2+1
      nfftdim3 = nfft3
      n = nfft3/2
      if ( nfft3 .eq. 2*n )nfftdim3 = nfft3+1
#ifdef SGIFFT
      nfftable = 2*(nfftdim1+nfftdim2+nfftdim3+50)
      nffwork = 0
      sizfftab = nfftable
      sizffwrk  = nffwork
#endif
#if defined _FFT_CRAY_
      nfftable = 2*(nfftdim1+nfftdim2+nfftdim3+50)
      nffwork = 4*nfftdim1*nfftdim2*nfftdim3
      sizfftab = nfftable
      sizffwrk  = nffwork
#elif defined _FFT_T3E_ & defined PARALLEL
      nfftable = 2*(nfftdim1+nfftdim2+nfftdim3+50)
      nffwork = 4*nfftdim1*nfftdim2*nfftdim3
      sizfftab = nfftable
      sizffwrk  = nffwork
#elif defined _FFT_T3E_ & !defined  PARALLEL
      nfftable = 12*2*(nfftdim1+nfftdim2+nfftdim3)
      nffwork = 2*nfftdim1*nfftdim2*nfftdim3
      sizfftab = nfftable
      sizffwrk  = nffwork
#elif defined _GPFFT_
      nfftable = 4*1024
      nffwork = 2*nfftdim1*nfftdim2*nfftdim3
      sizfftab = nfftable
      sizffwrk  = nffwork
#else
      nfftable = 4*nfftmax + 15
      nffwork = nfftmax
      sizfftab = 3*nfftable
      sizffwrk  = 2*nfftmax
#endif
      return
      end
