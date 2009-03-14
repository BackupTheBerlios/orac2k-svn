#ifdef PARALLEL
      SUBROUTINE timer(vfcp,tfcp,elapse)
      REAL*8 t2(2),elapse,vfcp,tfcp,tsecnd
      INCLUDE 'mpif.h'

      tfcp=MPI_WTIME()
      elapse=0.0D0
      RETURN
      END
#else
      SUBROUTINE timer(vfcp,tfcp,elapse)
      IMPLICIT none
      REAL*8 t2(2),elapse,vfcp,tfcp,h,m,s,yy,mm,dd
      integer value(8)
      CHARACTER*10 date,time,zone
      CALL date_and_time(date,time,zone,value)
      READ(date(3:4),'(f2.0)') yy
      READ(date(5:6),'(f2.0)') mm
      READ(date(7:8),'(f2.0)') dd
      READ(time(1:2),'(f2.0)') h
      READ(time(3:4),'(f2.0)') m
      READ(time(5:10),'(f6.3)') s
      tfcp=dd*3600*24+h*3600+m*60+s
      elapse=0.0D0
      RETURN
      END
#endif
