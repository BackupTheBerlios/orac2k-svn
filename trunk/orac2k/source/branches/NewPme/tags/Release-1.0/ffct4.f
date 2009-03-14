c---------------------------------------------------------------------
      subroutine ffct4(x,y,n)
c---------------------------------------------------------------------
c  Interface to essl/dxml/fftpack cosine fourier transform 
c  library routines 
c---------------------------------------------------------------------
      implicit none

#ifdef nolib
      real*4 wsave(60000),temp(0:30000)
#endif
#ifdef lessl
      integer nx1,nx2
      real*4  scl 
      real*8  aux1(30000),aux2(30000)
#endif
#ifdef ldxml
      INCLUDE 'DXMLDEF.FOR'
#endif
      integer*4   n,i
      real*4      x(0:*)
      real*4      y(0:*)
#ifdef ldxml     
      status=sfct('F',x,y,n,1,1)
#endif
#ifdef lessl
      scl=1.0
      nx1=30000
      nx2=30000
      call scosft(1,x,1,1,y,1,1,n,1,scl,aux1,nx1,aux2,nx2)
      call scosft(0,x,1,1,y,1,1,n,1,scl,aux1,nx1,aux2,nx2)
#endif
#ifdef nolib
      do i=0,n
         temp(i)=x(i)
      end do   
      call costi(n+1,wsave)
      call cost(n+1,temp,wsave)
      do i=0,n
         y(i)=temp(i)
      end do   
#endif      
      return
      end

