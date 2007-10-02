      subroutine mts_intra_n0(xp0,yp0,zp0,xcm,ycm,zcm,fpx,fpy,fpz
     &     ,ubond_slt,ubond_slv,ubend_slt,ubend_slv,uitors_slt
     &     ,uitors_slv,stressd,tag,nstart,nend,nlocal,ntot,node,nprocs
     &     ,ncube)

************************************************************************
*     FINTRAP stretching and bending routines for solute.
*     (used only with r-RESPA on)
************************************************************************

      
      USE Module_Extra_Forces, ONLY: Extra_force, Extra_Forces=>Force
      IMPLICIT none

      include 'parst.h'

*----------------------- ARGUMENTS -------------------------------------

      INTEGER nstart,nend,nlocal,ntot,tag(*)
      REAL*8  xp0(*),yp0(*),zp0(*),xcm(*),ycm(*),zcm(*),ubend_slt
     &     ,ubond_slt,uitors_slt,ubend_slv,ubond_slv,uitors_slv
     &     ,stressd(3,3)
      REAL*8  fpx(*),fpy(*),fpz(*)
      INTEGER node,nprocs,ncube

*-------------------- VARIABLES IN COMMONS -----------------------------

      INCLUDE 'cpropar.h'
      REAL*8 U_extra
      REAL*8  fppx(m1,2),fppy(m1,2),fppz(m1,2),work(m1)
      COMMON /rag1/ fppx,fppy,fppz,work

*-------------------- LOCAL VARIABLES ----------------------------------

#if defined T3E | defined _CRAY_
#define _AXPY_  saxpy
      INTEGER*8 n_loc,one
#else
#define _AXPY_  daxpy
      INTEGER n_loc,one
#endif
      REAL*8  xc,yc,zc,xd,yd,zd,sumx,sumy,sumz,st(3,3)
      INTEGER i,j,i1,count,m
      DATA one/1/

*=======================================================================
*--- Zero stress tensor ------------------------------------------------
*=======================================================================

      n_loc=nlocal
      IF(cpress .OR. pressure) THEN
         CALL zero0(stressd,9)
         CALL zero0(st,9)
      END IF

*=======================================================================
*----- Bonded interactions: stretching ------------------------------------
*=======================================================================

      IF(stretch) THEN
         CALL fpbond(ss_index,lstrtch,lstretch,lbnd_x,xp0,yp0,zp0
     &        ,potbo(1,2),potbo(1,1),ubond_slt,ubond_slv,fpx,fpy,fpz)
      END IF

*=======================================================================
*----- Bonded interactions: stretching ------------------------------------
*=======================================================================

      IF(Extra_Force) THEN
         CALL Extra_Forces(xp0,yp0,zp0,fpx,fpy,fpz,U_extra)
         IF(ubond_slt .NE. 0.0D0) THEN
            ubond_slt=ubond_slt+U_extra
         ELSE
            ubond_slt=U_extra
         END IF
      END IF


*=======================================================================
*----- Bonded interactions: bending ------------------------------------
*=======================================================================

      IF(bending) THEN
         CALL zeroa(fppx(1,1),fppy(1,1),fppz(1,1),ntot,1)
         CALL zeroa(fppx(1,2),fppy(1,2),fppz(1,2),ntot,1)
         CALL fpbend(ss_index,lbndg,lbend,tag,lbndg_x,xp0,yp0,zp0
     &        ,potbe(1,2),potbe(1,1),potbe(1,4),potbe(1,3),ubend_slt
     &        ,ubend_slv,fppx,fppy,fppz,m1)
         CALL _AXPY_(n_loc,1.0D0,fppx(nstart,2),one,fpx(nstart),one)
         CALL _AXPY_(n_loc,1.0D0,fppy(nstart,2),one,fpy(nstart),one)
         CALL _AXPY_(n_loc,1.0D0,fppz(nstart,2),one,fpz(nstart),one)
      END IF

*=======================================================================
*----- Bonded interactions: improper torsions --------------------------
*=======================================================================

      CALL fitors(ss_index,litr,litor,litr_x,xp0,yp0,zp0,potit(1,2),
     x     potit(1,1),potit(1,3),uitors_slt,uitors_slv,fpx,fpy,fpz)

*=======================================================================
*---- Compute stress tensor if coupling is by group --------------------
*=======================================================================

      IF((cpress .OR. pressure) .AND. coupl_grp) THEN
         DO j=nstart,nend
            i=atomp(j)
            xc=xcm(i)
            yc=ycm(i)
            zc=zcm(i)
            st(1,1)=st(1,1)+xc*fpx(j)
            st(1,2)=st(1,2)+yc*fpx(j)
            st(1,3)=st(1,3)+zc*fpx(j)
            st(2,1)=st(2,1)+xc*fpy(j)
            st(2,2)=st(2,2)+yc*fpy(j)
            st(2,3)=st(2,3)+zc*fpy(j)
            st(3,1)=st(3,1)+xc*fpz(j)
            st(3,2)=st(3,2)+yc*fpz(j)
            st(3,3)=st(3,3)+zc*fpz(j)
         END DO
         CALL DGEMM('N','N',3,3,3,1.0D0,co,3,st,3,0.0D0,stressd,3)
      END IF

*=======================================================================
*---- Add contribution to force from bendings among same group ---------
*=======================================================================

      IF(bending) THEN
         CALL _AXPY_(n_loc,1.0D0,fppx(nstart,1),one,fpx(nstart),one)
         CALL _AXPY_(n_loc,1.0D0,fppy(nstart,1),one,fpy(nstart),one)
         CALL _AXPY_(n_loc,1.0D0,fppz(nstart,1),one,fpz(nstart),one)
      END IF

#if defined PARALLEL
      IF(nprocs .GT. 1) THEN
         CALL P_merge_vecr8(stressd,9)
         CALL P_merge_r8(ubond_slt)
         CALL P_merge_r8(ubond_slv)
         CALL P_merge_r8(ubend_slt)
         CALL P_merge_r8(ubend_slv)
         CALL P_merge_r8(uitors_slt)
         CALL P_merge_r8(uitors_slv)
      END IF
#endif         

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
