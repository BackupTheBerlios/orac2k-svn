      SUBROUTINE verlet_free(nstart,nend,xp0,yp0,zp0,vpx,vpy,vpz,dt)

************************************************************************
*                                                                      *
*  VERLET advance coordinates for one time step and velocities for     *
*  half time step. Coordinates and velocities arrays are overwritten.  *
*  This program is NOT a black box. Velocities are passed in common    *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      integer      nstart,nend
      real*8       dt
      REAL*8       xp0(*),yp0(*),zp0(*)
      REAL*8       vpx(*),vpy(*),vpz(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8       tfact
      INTEGER      i

*==================== EXECUTABLE STATEMENTS ============================

      do i=nstart,nend
         xp0(i)=xp0(i)+dt*vpx(i)
         yp0(i)=yp0(i)+dt*vpy(i)
         zp0(i)=zp0(i)+dt*vpz(i)
      end do   


      RETURN
      END

