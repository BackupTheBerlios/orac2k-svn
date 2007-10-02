      SUBROUTINE get_total_energy(i_forces,mapnl,mapdn,nmapdn,tag_bndg
     &     ,abmd_dir,fudgec,xp0,yp0,zp0,fpx,fpy,fpz,stressd,stressr
     &     ,utotal,ucns,ucos,urcs,coul_bnd_slv,conf_bnd_slv_n1
     &     ,coul_bnd_slv_n1,self_slv,fsbond,fsbend,fsin14,unb14,cnb14
     &     ,uslvbon,uslvben,uslvtor,uslvitor,uumb,uptors,uitors,ubond
     &     ,ubend,ucnp,ucop,urcp,conf_bnd_slt_n1,coul_bnd_slt
     &     ,coul_bnd_slt_n1,self_slt,ucnsp,ucosp,urcsp,eer)

************************************************************************
*   Time-stamp: <98/02/21 11:05:29 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Feb 14 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER abmd_dir,mapnl(*),mapdn(2,*),nmapdn(*),tag_bndg(*)
      REAL*8  xp0(*),yp0(*),zp0(*),fpx(*),fpy(*),fpz(*)
     &     ,stressd(3,3),stressr(3,3),utotal,fsbend,fsbond
      LOGICAL i_forces

*----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'parst.h'
      INCLUDE 'cpropar.h'
      INCLUDE 'unit.h'
      INCLUDE 'parallel.h'

*----------------------- SCRATCH COMMON -------------------------------*

      INTEGER worka(m1)
      REAL*8 xpa(m1),ypa(m1),zpa(m1),xpcma(npm),ypcma(npm),zpcma(npm)
     &     ,xpcm(npm),ypcm(npm),zpcm(npm),xpga(m11),ypga(m11),zpga(m11)
     &     ,xpg(m11),ypg(m11),zpg(m11),fpx1(m1),fpy1(m1),fpz1(m1)

      COMMON /rag2/ fpx1,fpy1,fpz1,xpa,ypa,zpa,xpcm,ypcm,zpcm,xpcma
     &     ,ypcma,zpcma,xpga,ypga,zpga,xpg,ypg,zpg,worka

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8  puconf,pucoul,puhyd,pubnd,ubend,uptors,uitors,uconf,ucoul
     &     ,ureal,urecp,urcs,urcp,urcsp,urcsp_m,eer,ucns,ucos,ucnsp
     &     ,ucosp,ucnp,ucop,conf_bnd_slt,coul_bnd_slt,conf_bnd_slv
     &     ,coul_bnd_slv,self_slt,self_slv,uslvtor,uslvitor,fscnstr_slt
     &     ,fscnstr_slv,conf_bnd_slt_n1,coul_bnd_slt_n1,conf_bnd_slv_n1
     &     ,coul_bnd_slv_n1,virs,virsp,virp,gsin14,gsbend
     &     ,gsbond,unb14,cnb14,fsin14,ungrp,cngrp,uumb,gr,ubond,uslvbon
     &     ,uslvben,unbond,cnbond,purecp,upconf,upcoul,fudgec,press
     &     ,st(3,3),st_n1(3,3),st_n0(3,3),pv
     &     ,ustot,uptot,upstot
      INTEGER iret,npp,npp_u
      CHARACTER*80 errmsg

      CHARACTER*1 rshell,rshk
      DATA   rshell/'h'/rshk/'h'/

*----------------------- EXECUTABLE STATEMENTS ------------------------*


*=======================================================================
*---- length of neighbor list                                       ----
*=======================================================================

#if defined DYNAMIC_MEM
      npp_u=mpp/tgroup
#else
      npp_u=mpp
#endif
*=======================================================================
*----- Calculate solute center of mass coordinates and velocities ------
*=======================================================================

      CALL inicmp(ss_index,xp0,yp0,zp0,xpcm,ypcm,zpcm,mass,nprot
     &     ,protl)

*=======================================================================
*---- Calculate group position  ----------------------------------------
*=======================================================================

      CALL appbou(xp0,yp0,zp0,xpg,ypg,zpg,pmass,1,ngrp,grppt)
      
*=======================================================================
*-------- Find out the first and last group of each protein ------------
*=======================================================================

      CALL fndgrp(nprot,ngrp,protl,grppt,atomg,protg,groupp,atomp,npm
     &     ,iret,errmsg)
      IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)

*=======================================================================
*--- Change frame to get xpa, ypa, zpa etc in box fractions ------------
*=======================================================================

      CALL change_frame(co,oc,-1,ntap,xp0,yp0,zp0,xpa,ypa,zpa)
      CALL change_frame(co,oc,-1,ngrp,xpg,ypg,zpg,xpga,ypga,zpga)
      CALL change_frame(co,oc,-1,nprot,xpcm,ypcm,zpcm,xpcma,ypcma
     &     ,zpcma)

      urcs=0.0D0
      urcp=0.0D0
      urcsp=0.0D0
      eer=0.0D0
      self_slt=0.0D0
      self_slv=0.0D0
      fscnstr_slt=0.0D0
      fscnstr_slv=0.0D0
      fsin14=0.0D0
      unb14=0.0D0
      cnb14=0.0D0
      coul_bnd_slt=0.0D0
      coul_bnd_slv=0.0D0
      IF(clewld) THEN
         CALL cself(ss_index,ntap,alphal,rkcut,chrge,self_slt,self_slv)
         CALL ferrf(ss_index,alphal,chrge,1.0D0,xp0,yp0,zp0,0,lcnstr
     &        ,lconstr,lconstr_x,fscnstr_slt,fscnstr_slv,fpx1,fpy1,fpz1
     &        ,erf_corr,erf_arr_corr,delew,rlew)
#if defined PARALLEL
         IF(nprocs .GT. 1) THEN
            CALL P_merge_r8(fscnstr_slt)
            CALL P_merge_r8(fscnstr_slv)
         END IF
#endif
      END IF
      CALL zeroa(fpx1,fpy1,fpz1,ntap,1)

      CALL mts_forces(rshell,xpa,ypa,zpa,xpga,ypga,zpga,xpcma
     &     ,ypcma,zpcma,mapnl,mapdn,nmapdn,ucns,ucos,virs,virsp,ucnp
     &     ,ucop,ucnsp,ucosp,fpx1,fpy1,fpz1,stressd,worka,cpu_h,ncpu_h
     &     ,nstart_h,nend_h,nstart_ah,nend_ah,nlocal_ah,node,nprocs
     &     ,ncube,P_dyn_update_shell)
      IF(clewld) THEN
         CALL mts_furier(node,nodex,nodey,nodez,ictxt,npy,npz,descQ
     &        ,nprocs,ncube,nstart_1,nend_1,nlocal_1,nstart_2,nend_2
     &        ,nlocal_2,xp0,yp0,zp0,xpa,ypa,zpa,xpcma,ypcma,zpcma,urcsp
     &        ,urcs,urcp,virsp,virs,virp,fpx1,fpy1,fpz1,fsin14,fsbend
     &        ,fsbond,fscnstr_slt,fscnstr_slv,coul_bnd_slt,coul_bnd_slv
     &        ,rshell,rshk,eer,stressr,fudgec,tag_bndg)
      END IF
      CALL mts_intra_n1(xp0,yp0,zp0,xpcma,ypcma,zpcma,fpx1,fpy1
     &     ,fpz1,fudge,lj_fudge,abmd_dir,puhyd,conf_bnd_slt_n1
     &     ,conf_bnd_slv_n1,coul_bnd_slt_n1,coul_bnd_slv_n1,unb14,cnb14
     &     ,ungrp,cngrp,uptors,uslvtor,st_n1,mapdn,nmapdn,uumb,gr
     &     ,nstart_1,nend_1,node,nprocs,ncube)
      CALL mts_intra_n0(xp0,yp0,zp0,xpcma,ypcma,zpcma,fpx1,fpy1
     &     ,fpz1,ubond,uslvbon,ubend,uslvben,uitors,uslvitor
     &     ,st_n0,tag_bndg,nstart_1,nend_1,nlocal_1,ntot_1,node
     &     ,nprocs,ncube)

      IF(slv_exist) THEN
         uconf=ucns+ucos+urcs+coul_bnd_slv+conf_bnd_slv_n1
     &        +coul_bnd_slv_n1+self_slv + uumb
         ucoul=ucos+urcs+coul_bnd_slv+coul_bnd_slv_n1+self_slv
         ureal=ucos+coul_bnd_slv_n1
         urecp=urcs+self_slv+coul_bnd_slv
      END IF
      
*---        Energy terms for the protein         -----------------------
      
      IF(slt_exist) THEN 
         pubnd =uptors+uitors+ubond+ubend
         unbond=ucnp
         cnbond=ucop
         purecp=urcp
         puconf= unbond + conf_bnd_slt_n1 + uumb
         pucoul= cnbond + purecp + coul_bnd_slt + coul_bnd_slt_n1
     &        +self_slt
      END IF   
      
      
*---        Mixed terms                          -----------------------
      
      IF(slv_exist .AND. slt_exist) THEN 
         upconf= ucnsp+ucosp+urcsp
         upcoul= ucosp+urcsp
      END IF
      
      IF(pressure) THEN
         CALL daxpy(9,1.0D0,st_n1,1,stressd,1)
         CALL daxpy(9,1.0D0,st_n0,1,stressd,1)
      END IF
      
*----       If pme=T no solute-solvent separation ----------------------
      
      IF(pme) THEN 
         IF(slv_exist.AND.(.NOT.slt_exist)) THEN
            uconf=uconf + eer
            ucoul=ucoul + eer
            urecp= eer 
         END IF
         IF(slt_exist.AND.(.NOT.slv_exist)) THEN
            purecp=eer 
            pucoul=pucoul+eer 
         END IF
         IF(slt_exist.AND.slv_exist) THEN 
            upconf=upconf + eer 
            upcoul=upcoul + eer 
         END IF
      END IF
 
      ustot=0.d0
      ustot=(ucns+ucos+urcs+uslvbon+uslvben+uslvtor+uslvitor
     &     +conf_bnd_slv_n1+coul_bnd_slv+coul_bnd_slv_n1+self_slv)*efact
     &     /1000.0D0
      uptot=0.d0
      uptot=(urcp+ucop+ucnp+ubond+ubend+uptors+uitors
     &     +conf_bnd_slt_n1+coul_bnd_slt+coul_bnd_slt_n1+self_slt)*efact
     &     /1000.d0

      upstot=0.d0
      upstot = (urcsp+ucosp
     &     +ucnsp)*efact/1000.d0

      utotal= ustot+uptot+upstot+eer*efact/1000.d0
      IF(i_forces) THEN
         CALL dcopy(ntap,fpx1,1,fpx,1)
         CALL dcopy(ntap,fpy1,1,fpy,1)
         CALL dcopy(ntap,fpz1,1,fpz,1)
      END IF 

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
