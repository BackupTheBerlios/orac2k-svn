c---------------------------------------------------------------------
      subroutine fft4rc(x,y,n)
c---------------------------------------------------------------------
c  Interface to essl/dxml/default  real to complex fourier transform 
c  library routines. Sequence to be transformed must be
c  doubled padding zero from n/2+1 to n. Default routine 
c  is taken from the FFTPACK  by P.A. Swarztrauber
c---------------------------------------------------------------------
#ifdef lessl
      integer  nx1,nx2,nx3
      real*8   aux1(30000),aux2(30000),aux3
      real*4   scl 
#endif
#ifdef ldxml
      INCLUDE 'DXMLDEF.FOR'
#endif
#ifdef nolib
      real*4 temp(30000)
      dimension wsave(30000)
#endif
      integer*4   n
      real*4      x(*)
      complex*8   y(*)

#ifdef ldxml     
      status = sfft('R','C','F',x,y,n,1)
#endif
#ifdef lessl
      scl=1.0
      nx1=30000
      nx2=30000
      nx3=0.0
      call srcft(1,x,1,y,1,n,1,1,scl,aux1,nx1,aux2,nx2,aux3,nx3)
      call srcft(0,x,1,y,1,n,1,1,scl,aux1,nx1,aux2,nx2,aux3,nx3)
#endif
#ifdef nolib
      do i=1,n
         temp(i)=x(i)
      end do   
      call rffti(n,wsave)
      call rfftf(n,temp,wsave)
      do i=1,n
         y(i)=cmplx(temp(2*i-2),temp(2*i-1))
      end do   
#endif      
      return
      end
