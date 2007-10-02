      SUBROUTINE verlet(mass,nstart,nend,xp0,yp0,zp0,vpx,vpy,vpz,fpxn
     &     ,fpyn,fpzn,dt)

************************************************************************
*                                                                      *
*  VERLET advance coordinates for one time step and velocities for     *
*  half time step. Coordinates and velocities arrays are overwritten.  *
*  This program is NOT a black box. Velocities are passed in common    *
*                                                                      *
************************************************************************

!======================= DECLARATIONS ==================================

      IMPLICIT none

!----------------------- ARGUMENTS -------------------------------------

      integer      nstart,nend
      real*8       dt,mass(*)
      REAL*8       xp0(*),yp0(*),zp0(*)
      REAL*8       vpx(*),vpy(*),vpz(*)
      REAL*8       fpxn(*),fpyn(*),fpzn(*)

!-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8       ts2,tfact
      INTEGER      i

!==================== EXECUTABLE STATEMENTS ============================

      ts2=dt*dt

      do i=nstart,nend
         tfact=0.5*ts2/mass(i)
         xp0(i)=xp0(i)+dt*vpx(i)+tfact*fpxn(i)
         yp0(i)=yp0(i)+dt*vpy(i)+tfact*fpyn(i)
         zp0(i)=zp0(i)+dt*vpz(i)+tfact*fpzn(i)
         tfact=0.5*dt/mass(i)
         vpx(i)=vpx(i) + tfact*fpxn(i)
         vpy(i)=vpy(i) + tfact*fpyn(i)
         vpz(i)=vpz(i) + tfact*fpzn(i)
      end do   


      RETURN
      END

