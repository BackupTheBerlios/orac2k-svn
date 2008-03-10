************************************************************************
*   Time-stamp: <2006-02-05 16:00:29 marchi>                             *
*                                                                      *
*   Dump restart file and do tests at time H                           *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Jan 31 1998 -                                     *
*                                                                      *
************************************************************************

*=======================================================================
*---- Dump the restart file when required ------------------------------
*=======================================================================

      IF(nsave .NE. 0) THEN
         IF(mstep .NE. 0 .AND. MOD(mstep,nsave).EQ.0)THEN
#if  defined PARALLEL            
*=======================================================================
*---------  Expand group velocities in parallel runs ------------------
*=======================================================================

            CALL P_expand_r8x3(vcax,vcay,vcaz,nstart_cm,nend_cm
     &           ,nlocal_cm,node,nprocs,1)
            CALL P_expand_r8x3(vpx,vpy,vpz,nstart_2,nend_2,nlocal_2,node
     &           ,nprocs,1)

#endif
            IF(node .EQ. 0) THEN
               grflag= igr .OR. gofr
               WRITE(kprint,70200)
               CALL comp_vel_labframe(vpx,vpy,vpz,vpx1,vpy1,vpz1,co
     &              ,nprot,protl,vcax,vcay,vcaz)
               IF(restart_write) CALL dumprs(kdump_out,restart_out,nstep
     &              ,temp,ntap,ngrp,nprot,xp0,yp0,zp0,vpx1,vpy1,vpz1,eta
     &              ,etap,3,thermos,sumarray,ssmarray,navg,navg,co,oc
     &              ,cpress,vco,grflag,krdf,igrn)
            END IF
         END IF
      END IF

*=======================================================================
*--- Stop the run smoothly if a file called STOP is found --------------
*=======================================================================

      INQUIRE(FILE='./STOP',EXIST=exist)
#ifdef OSF1
      TimeToGo=TRemain()
      IF(TimeToGo .LT. 500) THEN
         exist=.TRUE.
         WRITE(kprint,70500)
      END IF
#endif
      IF(exist) THEN
#if  defined PARALLEL            
*=======================================================================
*---------  Expand group velocities in parallel runs ------------------
*=======================================================================

         CALL P_expand_r8x3(vcax,vcay,vcaz,nstart_cm,nend_cm,nlocal_cm
     &        ,node,nprocs,1)
         CALL P_expand_r8x3(vpx,vpy,vpz,nstart_2,nend_2,nlocal_2,node
     &        ,nprocs,1)
#endif
         IF(node .EQ. 0) THEN
            grflag= igr .OR. gofr   
            WRITE(kprint,70200)
            CALL matinv(3,3,co,oc,volume)
            volume=volume*boxl**3
            CALL comp_vel_labframe(vpx,vpy,vpz,vpx1,vpy1,vpz1,co,nprot
     &           ,protl,vcax,vcay,vcaz)
            IF(restart_write) CALL dumprs(kdump_out,restart_out,nstep
     &           ,temp,ntap,ngrp,nprot,xp0,yp0,zp0,vpx1,vpy1,vpz1,eta
     &           ,etap,3,thermos,sumarray,ssmarray,navg,navg,co,oc
     &           ,cpress,vco,grflag,krdf,igrn)
            IF(wrtgyr) THEN
               CLOSE(kgyr)
            END IF
            WRITE(kprint,70700)
         END IF
         STOP
      END IF
 
*=======================================================================
*---- Print out all energies and forces if (ltest_times)            ----
*=======================================================================

      IF(ltest_times)  THEN
         CALL GetKineticTest(nstart_2,nend_2,nstart_cm,nend_cm,ucek
     &        ,vpx,vpy,vpz,vcax,vcay,vcaz,tmass,mass,co)
         ucek=ucek*efact
         pucek=0.0D0
#if  defined PARALLEL
         CALL P_merge_r8(ucek)
#endif
         IF(thermos) THEN
            CALL comp_thermos_energy(neta,ndf_thermos,t,qmass,eta
     &           ,etap,uceh,hpot,temph)
         END IF
         CALL mts_test(node,efact,ntmtss,fpx_h,fpx_l,fpy_m,fpx_n1,ntmtsp
     &        ,time,urcs_h,urcs_l,urcs_m,ucos_h,ucos_l,ucos_m,ucns_h
     &        ,ucns_l,ucns_m,uslvbon,uslvben,uslvtor,uslvitor
     &        ,conf_bnd_slv_n1,coul_bnd_slv,coul_bnd_slv_n1,self_slv
     &        ,ucek,urcp_h,urcp_l,urcp_m,ucop_h,ucop_l,ucop_m,ucnp_h
     &        ,ucnp_l,ucnp_m,ubond,ubend,uitors,uptors,fsbond,fsbend
     &        ,fsin14,cnb14,unb14,cngrp,ungrp,conf_bnd_slt_n1
     &        ,coul_bnd_slt,coul_bnd_slt_n1,self_slt,pucek,urcsp_h
     &        ,urcsp_l,urcsp_m,ucosp_h,ucosp_l,ucosp_m,ucnsp_h,ucnsp_l
     &        ,ucnsp_m,eer_m,eer_l,eer_h,nmol,nato_slv,ntap,nstep,ktest
     &        ,lfirst,maxstp,cpress,volume,pext,ucepr,thermos,uceh,hpot)
      END IF
