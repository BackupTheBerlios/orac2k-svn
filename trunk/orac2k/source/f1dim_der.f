      DOUBLE PRECISION FUNCTION f1dim_der(iatom,icoord,mapnl,mapdn
     &     ,nmapdn,tag_bndg,fudgec,h,xp0,yp0,zp0,fpx
     &     ,fpy,fpz)

************************************************************************
*                                                                      *
*                                                                      *
************************************************************************


*======================= DECLARATIONS ==================================

      IMPLICIT none

      include 'parst.h'

*----------------------- ARGUMENTS -------------------------------------

      INTEGER iatom,icoord,mapnl(*),mapdn(2,*),nmapdn(*),tag_bndg(*)
      REAL*8  fudgec
      REAL*8  h,xp0(*),yp0(*),zp0(*),fpx(*),fpy(*),fpz(*)

*-------------------- VARIABLES IN COMMONS -----------------------------

      INCLUDE 'cpropar.h'

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i
      REAL*8  stressd(3,3),stressr(3,3),utotal,ucns,ucos,urcs
     &     ,coul_bnd_slv,conf_bnd_slv_n1,coul_bnd_slv_n1,self_slv
     &     ,uslvbon,uslvben,uslvtor,uslvitor,uumb,uptors,uitors,ubond
     &     ,ubend,ucnp,ucop,urcp,conf_bnd_slt_n1,coul_bnd_slt
     &     ,coul_bnd_slt_n1,self_slt,ucnsp,ucosp,urcsp,eer,fsbond,fsbend
     &     ,fsin14,unb14,cnb14

*----------- LOCAL WORK ARRAYS FOR THE RUN -----------------------------
  

*==================== EXECUTABLE STATEMENTS ============================

      IF(icoord .EQ. 1) THEN
         xp0(iatom)=xp0(iatom)+h
      ELSE IF(icoord .EQ. 2) THEN
         yp0(iatom)=yp0(iatom)+h
      ELSE IF(icoord .EQ. 3) THEN
         zp0(iatom)=zp0(iatom)+h
      END IF

      CALL get_total_energy(.TRUE.,mapnl,mapdn,nmapdn,tag_bndg,1,fudgec
     &     ,xp0,yp0,zp0,fpx,fpy,fpz,stressd,stressr,utotal,ucns,ucos
     &     ,urcs,coul_bnd_slv,conf_bnd_slv_n1,coul_bnd_slv_n1,self_slv
     &     ,fsbond,fsbend,fsin14,unb14,cnb14,uslvbon,uslvben,uslvtor
     &     ,uslvitor,uumb,uptors,uitors,ubond,ubend,ucnp,ucop,urcp
     &     ,conf_bnd_slt_n1,coul_bnd_slt,coul_bnd_slt_n1,self_slt,ucnsp
     &     ,ucosp,urcsp,eer)

      IF(icoord .EQ. 1) THEN
         xp0(iatom)=xp0(iatom)-h
      ELSE IF(icoord .EQ. 2) THEN
         yp0(iatom)=yp0(iatom)-h
      ELSE IF(icoord .EQ. 3) THEN
         zp0(iatom)=zp0(iatom)-h
      END IF
      f1dim_der=utotal
      
      RETURN
      END

