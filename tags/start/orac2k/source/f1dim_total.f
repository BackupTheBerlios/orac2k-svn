      DOUBLE PRECISION FUNCTION f1dim_total(mapnl,mapdn,nmapdn,tag_bndg
     &     ,fudgec,x,xp0,yp0,zp0,fpx,fpy,fpz)


************************************************************************
*                                                                      *
*                                                                      *
************************************************************************


*======================= DECLARATIONS ==================================

      IMPLICIT none

      include 'parst.h'

*----------------------- ARGUMENTS -------------------------------------

      INTEGER mapnl(*),mapdn(2,*),nmapdn(*),tag_bndg(*)
      REAL*8  fudgec
      REAL*8  xp0(*),yp0(*),zp0(*),fpx(*),fpy(*),fpz(*),x

*-------------------- VARIABLES IN COMMONS -----------------------------

      INCLUDE 'cpropar.h'

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i
      REAL*8  utotal

*----------- LOCAL WORK ARRAYS FOR THE RUN -----------------------------
  

*==================== EXECUTABLE STATEMENTS ============================


      DO i=1,ntap
         xp0(i)=xp0(i)+x*fpx(i)
         yp0(i)=yp0(i)+x*fpy(i)
         zp0(i)=zp0(i)+x*fpz(i)
      END DO
#if defined DYNAMIC_MEM
      CALL comp_total_energy(mapnl,mapdn,nmapdn,tag_bndg,fudgec,x,xp0
     &     ,yp0,zp0,utotal)
#else
      CALL comp_total_energy(mapnl,mapdn,nmapdn,tag_bndg,fudgec,x,xp0
     &     ,yp0,zp0,utotal)
#endif
      DO i=1,ntap
         xp0(i)=xp0(i)-x*fpx(i)
         yp0(i)=yp0(i)-x*fpy(i)
         zp0(i)=zp0(i)-x*fpz(i)
      END DO
      f1dim_total=utotal
      
      RETURN
      END

