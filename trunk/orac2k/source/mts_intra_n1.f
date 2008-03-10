      SUBROUTINE mts_intra_n1(xp0,yp0,zp0,xcm,ycm,zcm,fpx,fpy,fpz,phi
     &     ,fudge_qq,fudge_lj,fold_dir,puhyd,conf_bnd_slt,conf_bnd_slv
     &     ,coul_bnd_slt,coul_bnd_slv,unb14,cnb14,ungrp,cngrp,uptors
     &     ,uslvtor,stressd,mapdn,nmapdn,uumb,gr,nstart,nend,node,nprocs
     &     ,ncube)

************************************************************************
*     FINTRAP call slow "bonded" force routines for the solute.   
************************************************************************

      IMPLICIT none

      include 'parst.h'

*----------------------- ARGUMENTS -------------------------------------

      REAL*8  xp0(*),yp0(*),zp0(*),xcm(*),ycm(*),zcm(*),phi(*),unb14
     &     ,cnb14
      REAL*8  fpx(*),fpy(*),fpz(*),stressd(3,3),puhyd,conf_bnd_slt
     &     ,conf_bnd_slv,coul_bnd_slt,coul_bnd_slv,uslvtor
      INTEGER fold_dir,nmapdn(*),mapdn(2,*)
      REAL*8 uptors,fudge_qq,fudge_lj,ungrp,cngrp,uumb,gr
      INTEGER node,nprocs,ncube,nstart,nend

*-------------------- VARIABLES IN COMMONS -----------------------------

      include 'cpropar.h'
      REAL*8  work(m1)
      COMMON /rag1/ work

*-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8  unb14_slt,cnb14_slt,ungrp_slt,cngrp_slt,uptors_slt
     &       ,unb14_slv,cnb14_slv,ungrp_slv,cngrp_slv,uptors_slv
      REAL*8  xc,yc,zc,st(3,3)
      INTEGER i,j,k,i1,count,m

*=======================================================================
*--- Zero stress tensor ------------------------------------------------
*=======================================================================

      puhyd=0.0D0
      IF(cpress .OR. pressure) THEN
         CALL zero0(stressd,9)
         CALL zero0(st,9)
      END IF

*=======================================================================
*----- Third neighbour interactions : Sum on the direct lattice --------
*=======================================================================

      CALL fnb14(ss_index,xp0,yp0,zp0,chrge,ntap,ecc1412,ecc146
     &     ,cutoff,clewld,alphal,int14,int14p,type14,int14_x,fudge
     &     ,lj_fudge,fpx,fpy,fpz,phi,unb14_slt,unb14_slv,cnb14_slt
     &     ,cnb14_slv)
      unb14=unb14_slt
      cnb14=cnb14_slt

      conf_bnd_slt=unb14_slt
      coul_bnd_slt=cnb14_slt
      conf_bnd_slv=unb14_slv
      coul_bnd_slv=cnb14_slv

*=======================================================================
*----- Interactions between atoms contained in the same group ----------
*=======================================================================

      CALL fnbgrp(ss_index,xp0,yp0,zp0,chrge,nbtype,ecc12,ecc6,clewld
     &     ,alphal,ingrp,ingrpp,ingrp_x,fpx,fpy,fpz,phi,ungrp_slt
     &     ,ungrp_slv,cngrp_slt,cngrp_slv)
      ungrp=ungrp_slt
      cngrp=cngrp_slt

      conf_bnd_slt=conf_bnd_slt+ungrp_slt
      coul_bnd_slt=coul_bnd_slt+cngrp_slt
      conf_bnd_slv=conf_bnd_slv+ungrp_slv
      coul_bnd_slv=coul_bnd_slv+cngrp_slv

*=======================================================================
*----- Bonded interactions: proper torsions ----------------------------
*=======================================================================

      CALL fptors(ss_index,ltor,ltors,ltor_x,xp0,yp0,zp0,potto(1,2)
     &     ,potto(1,1),uptors_slt,uptors_slv,fpx,fpy,fpz)
      uptors =uptors_slt
      uslvtor=uptors_slv

*=======================================================================
*---- Calculate forces if umbrella on Rg is to be done -----------------
*=======================================================================

      IF(.NOT. fold) THEN
         IF(wrtgyr) THEN
            CALL mts_umbrl(abmd_unbias,fold,fold_dir,xp0,yp0,zp0,uumb
     &           ,nbone,mback,ntap,rspset,omega,spring,gr,fpx,fpy,fpz)
         END IF
      ELSE
         CALL mts_umbrl(abmd_unbias,fold,fold_dir,xp0,yp0,zp0,uumb,nbone
     &        ,mback,ntap,rspset,omega,spring,gr,fpx,fpy,fpz)
      END IF

*=======================================================================
*---- Compute stress tensor if coupling is by group --------------------
*=======================================================================

      IF(cpress .OR. pressure) THEN
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
         CALL DGEMM('N','T',3,3,3,1.0D0,st,3,co,3,0.0D0,stressd,3)
      END IF
#if defined PARALLEL
      IF(nprocs .GT. 1) THEN
         CALL P_merge_r8(unb14)
         CALL P_merge_r8(cngrp)
         CALL P_merge_r8(uptors)
         CALL P_merge_r8(uslvtor)
         CALL P_merge_r8(conf_bnd_slt)
         CALL P_merge_r8(conf_bnd_slv)
         CALL P_merge_r8(coul_bnd_slt)
         CALL P_merge_r8(coul_bnd_slv)
         IF(cpress .OR. pressure .AND. coupl_grp) THEN
            CALL P_merge_vecr8(stressd,9)
         END IF
      END IF
#endif         

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
