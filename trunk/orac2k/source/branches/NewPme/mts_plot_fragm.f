      subroutine mts_plot_fragm(i1,i2,beta,charge,xp0,yp0,zp0,ntap)
************************************************************************
*                                                                      *
*     Print snapshots of a given protein fragment                      *
*                                                                      *
************************************************************************

      implicit none

      real*8  xp0(*),yp0(*),zp0(*),co(3,3),charge(*),x,y,z,qq
      CHARACTER*7 beta(*)
      integer ntap,i,i1,i2

      INCLUDE 'unit.h'

*==================== EXECUTABLE STATEMENTS ============================

      do i=i1,i2
         x=xp0(i)
         y=yp0(i)
         z=zp0(i)
         qq=charge(i)*DSQRT(unitc)
         WRITE(kplot_fragm,1) beta(i),x,y,z,qq
      end do
1     FORMAT(5x,a7,4f16.9)
      
*================= END OF EXECUTABLE STATEMENTS ========================

      return
      end
