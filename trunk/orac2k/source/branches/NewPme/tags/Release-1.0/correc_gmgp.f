      subroutine correc_gmgp(vpx,vpy,vpz,gmgp,nmol)

************************************************************************
*                                                                      *
*     CORREC advances or velocities for half time step or corrects     *
*     velocities al full time step.                                    *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      REAL*8       gmgp(3,3)
      REAL*8       vpx(*),vpy(*),vpz(*)
      integer      nmol


*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER      i
      real*8       xc,yc,zc

*==================== EXECUTABLE STATEMENTS ============================

*=======================================================================
*  advance atomic velocities
*=======================================================================

      do i=1,nmol
         xc=gmgp(1,1)*vpx(i)+gmgp(1,2)*vpy(i)+gmgp(1,3)*vpz(i)
         yc=gmgp(2,1)*vpx(i)+gmgp(2,2)*vpy(i)+gmgp(2,3)*vpz(i)
         zc=gmgp(3,1)*vpx(i)+gmgp(3,2)*vpy(i)+gmgp(3,3)*vpz(i)
         vpx(i)=xc
         vpy(i)=yc
         vpz(i)=zc
      end do  

      RETURN
      END
