subroutine load_bsp_moduli(bsp_mod1,bsp_mod2,bsp_mod3,nfft1,nfft2,nfft3,order)
  implicit none
  integer nfft1,nfft2,nfft3,order
  double precision bsp_mod1(nfft1),bsp_mod2(nfft2),bsp_mod3(nfft3)
  integer MAXORDER
  parameter (MAXORDER=25)
  integer MAXNFFT
  parameter (MAXNFFT=1000)
  double precision array(MAXORDER),darray(MAXORDER),w
  double precision bsp_arr(MAXNFFT)
  integer i,maxn,opt_infl
!!$c this routine loads the moduli of the inverse DFT of the B splines
!!$c bsp_mod1-3 hold these values, nfft1-3 are the grid dimensions,
!!$c Order is the order of the B spline approx.
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
  do i = 1,maxn
     bsp_arr(i) = 0.d0
  END do
  do  i = 2,order+1
     bsp_arr(i) = array(i-1)
  END do
  call DFTMOD(bsp_mod1,bsp_arr,nfft1)
  call factor_lambda(bsp_mod1,nfft1,order,bsp_mod1,opt_infl)
  call DFTMOD(bsp_mod2,bsp_arr,nfft2)
  call factor_lambda(bsp_mod2,nfft2,order,bsp_mod2,opt_infl)
  call DFTMOD(bsp_mod3,bsp_arr,nfft3)
  call factor_lambda(bsp_mod3,nfft3,order,bsp_mod3,opt_infl)
  return
end subroutine load_bsp_moduli
subroutine factor_lambda(bsp_mod,nfft,order,prefac,opt_infl)
  implicit none
  integer nfft,order,opt_infl
  REAL*8 bsp_mod(nfft),prefac(nfft)
  integer KCUT
  parameter (KCUT=50)
!!$C factor in optimal lambda coefficient for bspline approximant
!!$c of complex exponentials, thus modify influence function
  integer k,nf,m,order2
  REAL*8 lambda,ratio,x,pi,gsum,gsum2
  nf = nfft / 2
  pi= 3.14159265358979323846
  order2 = 2*order
!!$c something wrong with get_tarray for large order.
!!$c use clunky but reliable gamma_sum
  do k = 1,nfft
     if ( opt_infl .eq. 0 )then
        lambda = 1.d0
     else
        m = k-1
        if ( k .gt. nf )m = k - 1 - nfft
        x = pi*dfloat(m)/nfft
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
end subroutine factor_lambda
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
  frac = dfloat(m)/nfft 
  x = pi*frac
  gsum = 1.d0
  do k = 1,kcut
     gsum = gsum + (x/(x + pi*k))**order
  enddo
  do k = 1,kcut
     gsum = gsum + (x/(x - pi*k))**order
  enddo
  return
end subroutine gamma_sum

!!$c------------------------------------------------------------------------
subroutine DFTMOD(bsp_mod,bsp_arr,nfft)
  implicit none
  integer nfft
  double precision bsp_mod(nfft),bsp_arr(nfft)
!!$c Computes the modulus of the discrete fourier transform of bsp_arr,
!!$c  storing it into bsp_mod
  
  integer j,k
  double precision sum1,sum2,twopi,arg,tiny
  twopi = 2.d0*3.14159265358979323846
  tiny = 1.d-7
  do k = 1,nfft
     sum1 = 0.d0
     sum2 = 0.d0
     do j = 1,nfft
        arg = twopi*(k-1)*(j-1)/nfft
        sum1 = sum1 + bsp_arr(j)*dcos(arg)
        sum2 = sum2 + bsp_arr(j)*dsin(arg)
     end do
     bsp_mod(k) = sum1**2 + sum2**2
  end do
  do k = 1,nfft
     if ( bsp_mod(k) .lt. tiny ) bsp_mod(k) = 0.5d0*(bsp_mod(k-1) + bsp_mod(k+1))
  END do
  return
end subroutine DFTMOD
subroutine fill_bspline2(w,order,array,darray,d2array)
!!$c---------- use standard B-spline recursions: see doc file
  implicit none
  integer order
  double precision w,array(order),darray(order),d2array(order)
  integer k
!!$c do linear case
  call init_bspline(array,w,order)
!!$c compute standard b-spline recursion
  do k = 3,order-2
     call one_pass_bspline(array,w,k)
  enddo
!!$c perform standard b-spline 2nd differentiation
  call diff_bspline(array,darray,order)
  call diff_bspline(darray,d2array,order)
!!$c perform standard b-spline differentiation
  call one_pass_bspline(array,w,order-1)
  call diff_bspline(array,darray,order)
!!$c one more recursion
  call one_pass_bspline(array,w,order)
  
  return
end subroutine fill_bspline2
subroutine diff_bspline(c,d,n)
  implicit none
  REAL*8 c(*),d(*)
  integer n
!!$c Using notation from Essmann et al; w = u-[u] and
!!$c array(j) = M_n(w + order - j)  where n is order
!!$c DERIVATIVE:    d/dw M_n(w) = M_n-1(w) - M_n-1(w-1)
!!$c i.e.   d/dw M_n(w+n-j) = M_n-1(w+n-j) - M_n-1(w+n-j-1)
!!$c i.e.   new(j) = old(j-1) - old(j)
!!$c where old is array before one_pass (thus n->n-1) and new is array afterwards
  
  integer j
  d(1) = -c(1)
  do j = 2,n
     d(j) = c(j-1) - c(j)
  enddo
  return
end subroutine diff_bspline
!!$c---------------------------------------------------------------------
subroutine one_pass_bspline(c,w,n)
  implicit none
  REAL*8 c(*),w
  integer n
!!$c Using notation from Essmann et al; w = u-[u] and
!!$c array(j) = M_n(w + order - j)  where n is order
!!$c RECURSION:  M_n(w) = (w/(n-1))*M_n-1(w)+((n-w)/(n-1))*M_n-1(w-1)
!!$c i.e.   M_n(w+n-j) = ((w+n-j)/(n-1))*M_n-1(w+n-j)+((j-w)/(n-1))*M_n-1(w+n-j-1)
!!$c i.e.   new(j) = ((w+n-j)/(n-1))*old(j-1) + ((j-w)/(n-1))*old(j)
!!$c where old is array before one_pass (thus n->n-1) and new is array afterwards
!!$c write backwards to do it with one array
  
  REAL*8 div
  integer j
  
  div = 1.d0 / (n-1)
  c(n) = div*w*c(n-1)
  do j = 1,n-2
     c(n-j) = div*((w+j)*c(n-j-1) + (n-j-w)*c(n-j))
  enddo
  c(1) = div*(1-w)*c(1)
  return
end subroutine one_pass_bspline
!!$c---------------------------------------------------------------------
subroutine init_bspline(c,w,n)
  implicit none
!!$c this gives back the spline array of order 2 (linear)
  integer n,i
  REAL*8 c(n),w
  c(2) = w
  c(1) = 1.d0 - w
  do i = 3,n
     c(i) = 0.d0
  enddo
  return
end subroutine init_bspline

subroutine fill_bspline(w,order,array,darray)
!!$c---------- use standard B-spline recursions: see doc file
  implicit none
  integer order
  double precision w,array(order),darray(order)
  
  integer k
!!$c do linear case
  call init(array,w,order)
!!$c compute standard b-spline recursion
  do k = 3,order-1
     call one_pass(array,w,k)
  end do
!!$c perform standard b-spline differentiation
  call diff(array,darray,order)
!!$c one more recursion
  call one_pass(array,w,order)
contains
!!$c---------------------------------------------------
subroutine init(c,x,order)
  implicit none
  integer order
  double precision c(order),x
  c(order) = 0.d0
  c(2) = x
  c(1) = 1.d0 - x
  return
end subroutine init
!!$c-------------------------------------
subroutine one_pass(c,x,k)
!DEC$ ATTRIBUTES FORCEINLINE :: one_pass
  implicit none
  double precision c(*),x
  integer k
  
  double precision div
  integer j  
  div = 1.d0 / (k-1)
  c(k) = div*x*c(k-1)
  do 100 j = 1,k-2
     c(k-j) = div*((x+j)*c(k-j-1) + (k-j-x)*c(k-j))
100  continue
     c(1) = div*(1-x)*c(1)
     return
end subroutine one_pass
!!$c-------------------------------------
  subroutine diff(c,d,order)
!DEC$ ATTRIBUTES FORCEINLINE :: diff
    implicit none
    double precision c(*),d(*)
    integer order
    
    integer j
    d(1) = -c(1)
    do j = 2,order
       d(j) = c(j-1) - c(j)
    end do
    return
  end subroutine diff
!!$c-------------------------------------
end subroutine fill_bspline


