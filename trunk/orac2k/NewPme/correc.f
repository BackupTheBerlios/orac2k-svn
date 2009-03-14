      subroutine correc(vpx,vpy,vpz,fpx,fpy,fpz,mass,nstart,nend,dt)

************************************************************************
*                                                                      *
*     CORREC advances or velocities for half time step or corrects     *
*     velocities al full time step.                                    *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      REAL*8       fpx(*),fpy(*),fpz(*)
      REAL*8       vpx(*),vpy(*),vpz(*)
      integer      nstart,nend
      real*8       dt,mass(*)


*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER      i
      real*8       tfact

*==================== EXECUTABLE STATEMENTS ============================

*=======================================================================
*  advance atomic velocities
*=======================================================================

      do i=nstart,nend
         tfact=0.5*dt/mass(i)
         vpx(i)=vpx(i) + fpx(i)*tfact
         vpy(i)=vpy(i) + fpy(i)*tfact
         vpz(i)=vpz(i) + fpz(i)*tfact
      end do  

      RETURN
      END
