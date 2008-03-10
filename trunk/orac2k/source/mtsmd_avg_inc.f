************************************************************************
*   Time-stamp: <2006-04-03 12:05:32 marchi>                             *
*                                                                      *
*   Compute averages and do some analysis at time step M               *
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Jan 31 1998 -                                     *
*                                                                      *
************************************************************************


*=======================================================================
*---        Computes all K.E.s, temperatures 
*=======================================================================

            IF(coupl_grp) THEN
               ptr=get_pointer_protl(nstart_cm-1,protl)
               CALL comp_vel_labframe(vpx,vpy,vpz,vpx1,vpy1,vpz1,co
     &              ,nlocal_cm,protl(ptr),vcax(nstart_cm),vcay(nstart_cm
     &              ),vcaz(nstart_cm))
               IF(nprocs .EQ. 1) CALL comp_vcm(vpx1,vpy1,vpz1,oc,nprotb
     &              ,protlb,mass,tmassb,vcbx,vcby,vcbz)
               CALL kinetic(nstart_2,nend_2,ss_index,co,nato_slv
     &              ,nmol,cnstpp_slv,ntap,cnstpp,xp0,yp0,zp0,vpx1,vpy1
     &              ,vpz1,vcbx,vcby,vcbz,nprotb,protlb,mass,tmassb,wmtp
     &              ,ucek,pucek,temp,tempt,tempr,temppr,tcm,rcm,stresstk
     &              ,cpress,isostress,vco,masspp,ucepr,temppra)
#if  defined PARALLEL
               CALL P_merge_r8(ucek)
               CALL P_merge_r8(pucek)
               CALL P_merge_r8(temp)
               CALL P_merge_r8(temppr)
#endif
            ELSE
               CALL kinetic(1,ntap,ss_index,co,nato_slv,nmol
     &              ,cnstpp_slv,ntap,cnstpp,xp0,yp0,zp0,vpx,vpy,vpz,vcax
     &              ,vcay,vcaz,nprot,protl,mass,tmass,wmtp,ucek,pucek
     &              ,temp,tempt,tempr,temppr,tcm,rcm,stresstk,cpress
     &              ,isostress,vco,masspp,ucepr,temppra)
            END IF            
            
*=======================================================================
*---        Gather all energy terms 
*=======================================================================
            
            
*--         Energy terms for the solvent         -----------------------
            
1010        CONTINUE
            fsrtalc=0.0D0
            IF(slv_exist) THEN
               uconf=ucns_h+ucos_h+ucns_l+ucos_l+ucns_m+ucos_m+urcs_h
     &              +urcs_l+urcs_m+coul_bnd_slv+conf_bnd_slv_n1
     &              +coul_bnd_slv_n1+self_slv + uumb
               ucoul=ucos_h+ucos_l+ucos_m+urcs_h+urcs_l+urcs_m
     &              +coul_bnd_slv+coul_bnd_slv_n1+self_slv
               ureal=ucos_h+ucos_l+ucos_m+coul_bnd_slv_n1
               urecp=urcs_h+urcs_l+urcs_m+self_slv+coul_bnd_slv
            END IF
            
*---        Energy terms for the protein         -----------------------
            
            IF(slt_exist) THEN 
               pubnd=uptors+uitors+ubond+ubend
               unbond=ucnp_m + ucnp_l + ucnp_h + virtual_energy
               cnbond=ucop_m + ucop_l + ucop_h + virtual_energy
               purecp=urcp_h+urcp_l+urcp_m
               puconf= unbond + conf_bnd_slt_n1 + uumb
               pucoul= cnbond + purecp + coul_bnd_slt + coul_bnd_slt_n1
     &              +self_slt
            END IF   
            
*---        Mixed terms                          -----------------------
            
            IF(slv_exist .AND. slt_exist) THEN 
               upconf= ucnsp_h + ucnsp_l + ucnsp_m + 
     &              ucosp_h + ucosp_l + ucosp_m + 
     &              urcsp_h + urcsp_l + urcsp_m
               upcoul= ucosp_h + ucosp_l + ucosp_m + 
     &              urcsp_h + urcsp_l + urcsp_m
            END IF
            IF(pressure) THEN
               DO i=1,3
                  DO j=1,3
                     stressd(i,j)=stressd_h(i,j)+stressd_l(i,j)
     &                    +stressd_m(i,j)+stressd_n1(i,j)
     &                    +stressd_n0(i,j)
                     stressr(i,j)=stressr_h(i,j)+stressr_l(i,j)
     &                    +stressr_m(i,j)
                  END DO
               END DO

               CALL comp_stress_conf(stressd,stressr,stress_conf,oc
     &              ,volume,unitp,press_conf)
               CALL comp_stress_kinetic(vcax,vcay,vcaz,tmass,co
     &              ,nstart_cm,nend_cm,volume,unitp,stress_kin,press_kin
     &              )
#if  defined PARALLEL            
               CALL P_merge_r8(press_kin)
               CALL P_merge_vecr8(stress_kin,9)
#endif
               CALL dcopy(9,stress_kin,1,stress_tot,1)
               CALL daxpy(9,1.0D0,stress_conf,1,stress_tot,1)
            END IF
            
*----       If pme=T no solute-solvent separation ----------------------
            
            IF(pme) THEN 
               IF(polar) THEN
                  eer_m=eer_m+uself_dip+uself
               END IF
               IF(slv_exist.AND.(.NOT.slt_exist)) THEN
                  uconf=uconf + eer_h + eer_l + eer_m
                  ucoul=ucoul + eer_h + eer_l + eer_m
                  urecp= eer_h + eer_l + eer_m
               END IF
               IF(slt_exist.AND.(.NOT.slv_exist)) THEN
                  purecp=eer_h + eer_l + eer_m
                  pucoul=pucoul+eer_h + eer_l + eer_m
               END IF
               IF(slt_exist.AND.slv_exist) THEN 
                  upconf=upconf + eer_h + eer_l + eer_m
                  upcoul=upcoul + eer_h + eer_l + eer_m
               END IF
               eer=eer_h + eer_l + eer_m
            END IF
            IF(energy_then_die) go to 1020
            
*=======================================================================
*---   Now for each timestep of shell M do inline analysis, scale    ---
*---     temperatures, compute averages and dump hystory files       ---
*=======================================================================
            

            IF(node .EQ. 0) THEN
*=======================================================================
*---- Compute instantaneous X-rms --------------------------------------
*=======================================================================
            
               IF(anxrms) CALL CalcXrmsSecStruct(anxca,anxbc,anxhe,anxal
     &              ,SecStructure,SecStructTotal,SecPointer
     &              ,grppt,mres,protlb,wca,whe,wbc,xp0,yp0,zp0
     &              ,xpt0,ypt0,zpt0,nato_slt,errca,errhe,errbc
     &              ,erral,drpca,drpbc,drphe,drpal)

*=======================================================================
*---- Compute averaged structure ---------------------------------------
*=======================================================================
            
               IF(avg_str) THEN
                  IF(avg_ca) CALL calc_avg_str(anprot,annpro
     &                ,anpoint,protlb,wca,xpt0,ypt0,zpt0,xp_avg
     &                ,yp_avg,zp_avg
     &                ,qt,xp0,yp0,zp0,nato_slt,iter_avg)
                  IF(avg_he) CALL calc_avg_str(anprot,annpro
     &                ,anpoint,protlb,whe,xpt0,ypt0,zpt0,xp_avg
     &                ,yp_avg,zp_avg
     &                ,qt,xp0,yp0,zp0,nato_slt,iter_avg)
                  IF(avg_bc) CALL calc_avg_str(anprot,annpro
     &                ,anpoint,protlb,wbc,xpt0,ypt0,zpt0,xp_avg
     &                ,yp_avg,zp_avg
     &                ,qt,xp0,yp0,zp0,nato_slt,iter_avg)
               END IF
            
               IF(navg_str .NE. 0) THEN
                  IF(MOD(ninner,navg_str) .EQ.0) THEN
                     CALL dcopy(nato_slt,xp_avg,1,xpo,1)
                     CALL dcopy(nato_slt,yp_avg,1,ypo,1)
                     CALL dcopy(nato_slt,zp_avg,1,zpo,1)
                     fact=1.0D0/DFLOAT(iter_avg)
                     CALL dscal(nato_slt,fact,xpo,1)
                     CALL dscal(nato_slt,fact,ypo,1)
                     CALL dscal(nato_slt,fact,zpo,1)
                     fstep=time*DFLOAT(ninner)/dfloat(mrespa*lrespa)
                     IF(avg_ca) WRITE(kavg,80000) 
                     IF(avg_he) WRITE(kavg,80100)
                     IF(avg_bc) WRITE(kavg,80101) 
                     CALL plotd(fstep,kavg,beta,xpt0,ypt0,zpt0,xpo,ypo
     &                    ,zpo,nato_slt,nres,m1,prsymb)
                  END IF
               END IF
               IF(navg_str_xrms .NE. 0) THEN
                  IF(MOD(ninner,navg_str_xrms) .EQ.0) THEN
                     CALL dcopy(nato_slt,xp_avg,1,xpo,1)
                     CALL dcopy(nato_slt,yp_avg,1,ypo,1)
                     CALL dcopy(nato_slt,zp_avg,1,zpo,1)
                     fact=1.0D0/DFLOAT(iter_avg)
                     CALL dscal(nato_slt,fact,xpo,1)
                     CALL dscal(nato_slt,fact,ypo,1)
                     CALL dscal(nato_slt,fact,zpo,1)
                     fstep=time*DFLOAT(ninner)/dfloat(mrespa*lrespa)
                     CALL calc_avg_xrms(avg_ca,avg_he,avg_bc,fstep
     &                    ,kavg_xrms,xpt0,ypt0,zpt0,xpo,ypo,zpo
     &                    ,wca,whe,wbc,protlb
     &                    ,anprot,annpro,anpoint,nato_slt)
                  END IF
               END IF
               
               IF(debug) THEN
                  IF(.NOT. ASSOCIATED(fpx)) THEN
                     ALLOCATE(fpx(ntap),fpy(ntap),fpz(ntap))
                  END IF
                  DO i=1,ntap
                     fpx(i)=fpx_n0(i)+fpx_n1(i)+fpx_m(i)+fpx_l(i)
     &                    +fpx_h(i)
                     fpy(i)=fpy_n0(i)+fpy_n1(i)+fpy_m(i)+fpy_l(i)
     &                    +fpy_h(i)
                     fpz(i)=fpz_n0(i)+fpz_n1(i)+fpz_m(i)+fpz_l(i)
     &                    +fpz_h(i)
                  END DO
                  CALL prtfrc(kprint,ngrp,grppt,nres,M1,prsymb,beta,xp0
     &                 ,yp0,zp0,fpx,fpy,fpz)
               END IF
            END IF
               
*=========== Write instantaneous results, averages and fluctuations ====
c ==
            

1020        CONTINUE 
            CALL prtacc(node,pucek,puhyd,puconf,pueng,pucoul,self_slt
     &           ,fsbond,purecp,fsbend,fsin14,unb14,cnb14,ubend,ubond
     &           ,uitors,uptors,pubnd,uceh,hpot,ucoul,uconf,urecp,ureal
     &           ,fsrtal,ucek,upconf,upcoul,uslvbon,uslvben,uslvtor
     &           ,uslvitor,eer,Uind,U_Thole,U_solv,uumb,temp,temph,tcm
     &           ,rcm,tempt,tempr,temppr,gr,gra,ucepr,stress_tot
     &           ,press_conf,pressc,press_kin,temppra,errca,errhe,errbc
     &           ,erral,drpca,drpbc,drphe,drpal,sum_econf,sum_ecoul
     &           ,sum_enbnd,sum_etotpot,sum_eslvint,sum_eebond
     &           ,sum_eebend,sum_eeptors,sum_eeitors,sum_tote,sum_ucek
     &           ,sum_temp,sum_tempt,sum_tempr,sum_temppr,sum_temph
     &           ,sum_pecek,sum_pehyd,sum_peconf,sum_pecoul,sum_percip
     &           ,sum_enb14,sum_ebend,sum_ebond,sum_eitor,sum_eptor
     &           ,sum_pnbd,sum_pebnd,sum_pepot,sum_ptote,sum_gr
     &           ,sum_epcoul,sum_epconf,sum_co,sum_st,sum_presst
     &           ,sum_press,sum_pressc,sum_pressk,sum_volume,sum_pv
     &           ,sum_temppra,ssm_econf,ssm_ecoul,ssm_enbnd,ssm_etotpot
     &           ,ssm_eslvint,ssm_eebond,ssm_eebend,ssm_eeptors
     &           ,ssm_eeitors,ssm_tote,ssm_ucek,ssm_temp,ssm_tempt
     &           ,ssm_tempr,ssm_temppr,ssm_temph,ssm_pecek,ssm_pehyd
     &           ,ssm_peconf,ssm_pecoul,ssm_percip,ssm_enb14,ssm_ebend
     &           ,ssm_ebond,ssm_eitor,ssm_eptor,ssm_pnbd,ssm_pebnd
     &           ,ssm_pepot,ssm_ptote,ssm_gr,ssm_epcoul,ssm_epconf
     &           ,ssm_co,ssm_st,ssm_presst,ssm_press,ssm_pressc
     &           ,ssm_pressk,ssm_volume,ssm_pv,ssm_temppra,energy,ninner
     &           ,nstep)
            IF(energy_then_die) STOP
            
            
*=========== Dump files ================================================
c =
            IF(node .EQ. 0) THEN
               IF(ndipole.gt.0) THEN
                  IF(MOD(ninner,ndipole).EQ.0)THEN
                     CALL comp_dip(co,xpga,ypga,zpga,xpa,ypa,zpa
     &                    ,chrge,dips,ntap,ngrp,grppt)
                     fstep=time*DFLOAT(ninner)/dfloat(mrespa*lrespa)
                     aux=elechg*unitl/3.336D-30
                     WRITE(kdipole,106) fstep, (dips(j)*aux,j=1,3)
106                  FORMAT(' Dip. ',f12.3,' t ', 3e15.5)
                  END IF
               END IF
               IF(nascii .NE. 0) THEN
                  IF(MOD(ninner,nascii) .EQ. 0) THEN
                     IF(ascii_nocell) THEN
                        CALL change_frame(co,oc,1,ntap,xpa,ypa,zpa,xpo
     &                       ,ypo,zpo)
                     ELSE IF(ascii_wsc) THEN
                        CALL tr_wsc(co,xpa,ypa,zpa,xpo,ypo,zpo,mass
     &                       ,nprotb,protlb)
                        CALL change_frame(co,oc,1,ntap,xpo,ypo,zpo,xpo
     &                       ,ypo,zpo)
                     ELSE
                        CALL tr_inbox(xpa,ypa,zpa,xpo,ypo,zpo,mass
     &                       ,nprotb,protlb)
                        CALL change_frame(co,oc,1,ntap,xpo,ypo,zpo,xpo
     &                       ,ypo,zpo)
                     END IF
                     fstep=time*DFLOAT(ninner)/dfloat(mrespa*lrespa)
                     CALL plotc(co,abmd,gr,gra,fstep,beta,xpo,ypo,zpo
     &                    ,ntap,nres,m1,prsymb,chrge)
                  END IF
               END IF
               
*=========== Dump configurations to a .PDB file ========================
c =
               
               IF(nplot.GT.0) THEN
                  IF(MOD(ninner,nplot).EQ.0)THEN
                     CALL mts_plotp(beta,mback,nbone,xp0,yp0,zp0,ntap
     &                    ,nres,m1,prsymb)
                  END IF
               END IF
               
               IF(nplot_fragm .GT. 0) THEN
                  IF(MOD(ninner,nplot_fragm).EQ.0 ) THEN
                     fstep=time*DFLOAT(ninner)/dfloat(mrespa*lrespa)
                     write(kplot_fragm,'(''Fstep '',f20.5,2i8)') fstep
     &                    ,ninner
                     do i=1,nfragm
                        write(kplot_fragm,'(''Fragm No. '',i10)') i
                        fragm_1=fragm(1,i)
                        fragm_2=fragm(2,i)
                        CALL mts_plot_fragm(fragm_1,fragm_2,beta,chrge
     &                       ,xp0,yp0,zp0,ntap)
                     end do
                  END IF
               END IF
               
               IF(nplot_center .NE. 0) THEN
                  IF(MOD(ninner,nplot_center) .EQ. 0) THEN
                     CALL tr_inbox(xpa,ypa,zpa,xpo,ypo,zpo,mass,nprotb
     &                    ,protlb)
                     CALL change_frame(co,oc,1,ntap,xpo,ypo,zpo,xpo,ypo
     &                    ,zpo)
                     fstep=time*DFLOAT(ninner)/dfloat(mrespa*lrespa)
                     CALL plot_center(abmd,gr,gra,fstep,beta,xpo,ypo,zpo
     &                    ,ntap,nres,m1,prsymb,chrge)
                  END IF
               END IF

               IF(gofr) THEN
                  IF(MOD(ninner,gofr_ncomp) .EQ. 0 .AND. nstep .GT.
     &                 nrject)THEN
                     l2=maxint
                     CALL calc_gofr(nato_slt,nato_slv,type_slv,ss_index
     &                    ,atomp,gofr_neighbor,gofr_intra,co,xpa,ypa,zpa
     &                    ,xpga,ypga,zpga,wca,whe,delrg,nstart_h,nend_h
     &                    ,nstart_ah,nend_ah,ntap,ngrp,grppt,l2,krdf
     &                    ,ngrdon,gofr_cut)
                     IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
                  END IF
               END IF
               
               IF(gofr) THEN
                  IF(MOD(ninner,gofr_nprint) .EQ. 0 .AND. nstep .GT.
     &                 nrject)THEN
                     fstep=time*DFLOAT(ninner)/dfloat(mrespa*lrespa)
                     vol_gofr=sum_volume/DFLOAT(ninner)
                     IF(slt_exist) THEN
                        offset=0
                        CALL write_gofrp(.NOT.gofr_avg,fstep,krdf,maxint
     &                       ,offset,wca,whe,nato_slt,delrg,gofr_cut
     &                       ,ngrdon,iret,errmsg)
                        IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
                     END IF
                     IF(slv_exist .AND. slt_exist) THEN
                        CALL write_gofrw(.FALSE.,fstep,krdf,maxint,g1
     &                       ,delrg,gofr_cut,itype_slv,betab_slv
     &                       ,vol_gofr,ntype_slv,nmol,ngrdon,3,iret
     &                       ,errmsg)
                        IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
                     ELSE IF(slv_exist .AND. (.NOT. slt_exist)) THEN
                        CALL write_gofrw(.TRUE.,fstep,krdf,maxint,g1
     &                       ,delrg,gofr_cut,itype_slv,betab_slv
     &                       ,vol_gofr,ntype_slv,nmol,ngrdon,3,iret
     &                       ,errmsg)
                        IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
                     END IF
                  END IF
               END IF
               
*=======================================================================
*---- Initialize G of Rs when required ---------------------------------
*=======================================================================
               
               IF(gofr .AND. gofr_avg  .AND. nstep .GT. nrject) THEN
                  IF(MOD(ninner,gofr_navg) .EQ. 0) THEN
                     CALL zero_gofr(maxint,krdf,ngrdon,offset_slv)
                  END IF
               END IF
               
*=======================================================================
*----- Write topology --------------------------------------------------
*=======================================================================
               
               IF(prttopl) THEN
                  IF(MOD(ninner,ntop_print) .EQ. 0) THEN
                     fstep=time*DFLOAT(ninner)/dfloat(mrespa*lrespa)
                     IF(top_bonds(1) .GT. 0) CALL write_bonds(ktopol
     &                    ,fstep,top_bonds,lbnd,lbond,xp0,yp0,zp0)
                     IF(top_bendings(1) .GT. 0) CALL write_bends(ktopol
     &                    ,fstep,top_bendings,lbndg,lbend,xp0,yp0,zp0)
                     IF(top_ptors(1) .GT. 0) CALL write_tors('P',ktopol
     &                    ,fstep,top_ptors,ltor,ltors,xp0,yp0,zp0)
                     IF(top_itors(1) .GT. 0) CALL write_tors('I',ktopol
     &                    ,fstep,top_itors,litr,litor,xp0,yp0,zp0)
                  END IF
               END IF
               
*=========== Dump the hystory of the run when required =================
c =
               IF(nconf.GT.0 .AND. node .EQ. 0) THEN
                  IF(MOD(ninner,nconf).EQ.0 .AND. ninner.GT.mrject)THEN
                     fstep=time*DFLOAT(ninner)/dfloat(mrespa*lrespa)
                     CALL write_confc(co,xp0,yp0,zp0,ntap,fstep,ninner
     &                    ,nconf,divide_records,atom_record)
                  END IF
               END IF
            END IF
            
*=========== Scale the temperature when required =======================
               
            IF(ninner-mrject.EQ.0)THEN
               
               IF(coupl_grp) THEN
                  ptr=get_pointer_protl(nstart_cm-1,protl)
                  CALL comp_vel_labframe(vpx,vpy,vpz,vpx1,vpy1,vpz1,co
     &                 ,nlocal_cm,protl(ptr),vcax(nstart_cm)
     &                 ,vcay(nstart_cm),vcaz(nstart_cm))
                  IF(nprocs .EQ. 1) CALL comp_vcm(vpx1,vpy1,vpz1,oc
     &                 ,nprotb,protlb,mass,tmassb,vcbx,vcby,vcbz)
                  CALL kinetic(nstart_2,nend_2,ss_index,co,nato_slv
     &                 ,nmol,cnstpp_slv,ntap,cnstpp,xp0,yp0,zp0,vpx1
     &                 ,vpy1,vpz1,vcbx,vcby,vcbz,nprotb,protlb,mass
     &                 ,tmassb,wmtp,ucek,pucek,temp,tempt,tempr,temppr
     &                 ,tcm,rcm,stresstk,cpress,isostress,vco,masspp
     &                 ,ucepr,temppra)
#if  defined PARALLEL            
                  CALL P_merge_r8(ucek)
                  CALL P_merge_r8(pucek)
                  CALL P_merge_r8(temp)
                  CALL P_merge_r8(temppr)
#endif
               ELSE 
                  CALL kinetic(1,ntap,ss_index,co,nato_slv,nmol
     &                 ,cnstpp_slv,ntap,cnstpp,xp0,yp0,zp0,vpx,vpy,vpz
     &                 ,vcax,vcay,vcaz,nprot,protl,mass,tmass,wmtp,ucek
     &                 ,pucek,temp,tempt,tempr,temppr,tcm,rcm,stresstk
     &                 ,cpress,isostress,vco,masspp,ucepr,temppra)
               END IF
               
               IF(slv_exist) ustot=(ucns_h+ucos_h+ucns_l+ucos_l+ucns_m
     &              +ucos_m+urcs_h+urcs_l+urcs_m+uslvbon+uslvben
     &              +uslvtor+uslvitor+coul_bnd_slv+self_slv
     &              )*efact/1000.0D0
               
               uptot=0.d0
               IF(slt_exist) uptot=(urcp_h+urcp_l+urcp_m+ucop_h+ucop_l
     &              +ucop_m+ucnp_h+ucnp_l+virtual_energy+ucnp_m+ubond
     &              +ubend+uptors+uitors+coul_bnd_slt+self_slt)*efact
     &              /1000.d0
               
               ektot = (ucek+pucek)/1000.d0
               nscal=0
               
               upstot=0.d0
               IF(slv_exist .AND. slt_exist)
     &              upstot = (urcsp_h+urcsp_l+urcsp_m+ucosp_h+ucosp_l
     &              +ucosp_m+ucnsp_h+ucnsp_l+ucnsp_m+eer_m+eer_l+eer_h)
     &              *efact/1000.d0
               E0 = ektot+ustot+uptot+upstot
               nrject=-1
               mrject=-1
               ninner = 0
               nstep=0
               WRITE(kprint,13000) nscale
            END IF
            IF(annealing .AND. ninner-mrject .GT. 0) THEN
               CALL anneal(annealing_fact,vpx,vpy,vpz,ntap)
               IF(cnstpp .NE. 0) THEN
                  CALL rattle_correc(nstart_2,nend_2,tm,xp0,yp0,zp0
     &                 ,vpx,vpy,vpz,ntap,cnstp,dssp,coeffp,cnstpp,mass
     &                 ,dnit,cnst_protp,cnst_protl,mim_lim,gcpu_rt
     &                 ,iret,errmsg)
                  IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
               END IF
            END IF
            
            IF(ninner-mrject .LT. 0) THEN
               
               IF(DABS(t-temp).GT.dtemp .AND. (.NOT. heating)) THEN
                  nscale=nscale+1
                  CALL change_origin(-1,nprot,protl,vpx,vpy,vpz,vpx,vpy
     &                 ,vpz,vcax,vcay,vcaz,co)
                  CALL ranvel(t,mass,ntap,vpx,vpy,vpz,xp0,yp0,zp0,co,
     &                 .TRUE.)
                  IF(cnstpp .NE. 0) THEN
                     CALL rattle_correc(nstart_2,nend_2,tm,xp0,yp0
     &                    ,zp0,vpx,vpy,vpz,ntap,cnstp,dssp,coeffp,cnstpp
     &                    ,mass,dnit,cnst_protp,cnst_protl,mim_lim
     &                    ,gcpu_rt,iret,errmsg)
                     IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
                  END IF
                  CALL comp_vcm(vpx,vpy,vpz,oc,nprot,protl,mass,tmass
     &                 ,vcax,vcay,vcaz)
                  IF(.not.start_conf) WRITE(kprint,70100)
               END IF
               IF(heating .AND. MOD(ninner,mrespa*lrespa) .EQ. 0) THEN
                  nscale=nscale+1
                  CALL change_origin(-1,nprot,protl,vpx,vpy,vpz,vpx,vpy
     &                 ,vpz,vcax,vcay,vcaz,co)

                  CALL HEAT_Set_slv(ucek,vpx,vpy,vpz)
                  CALL HEAT_Set_slt(temp_heat,pucek,vpx,vpy,vpz)

                  IF(cnstpp .NE. 0) THEN
                     CALL rattle_correc(nstart_2,nend_2,tm,xp0,yp0
     &                    ,zp0,vpx,vpy,vpz,ntap,cnstp,dssp,coeffp,cnstpp
     &                    ,mass,dnit,cnst_protp,cnst_protl,mim_lim
     &                    ,gcpu_rt,iret,errmsg)
                     IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
                  END IF
                  CALL comp_vcm(vpx,vpy,vpz,oc,nprot,protl,mass,tmass
     &                 ,vcax,vcay,vcaz)

                  IF(temp_heat .LE. t) temp_heat=temp_heat+dtemp_heat
               END IF

               IF(cpress) THEN
                  IF(DABS(t-temppra) .GT. dtemppr .OR. MOD(ninner,scale)
     &                 .EQ. 0) THEN
                     CALL set_tempp(masspr,vco,temppra,t)
                     IF(isostress   .OR. FixedAngles_Stress) THEN
                        CALL rattle_correc_co(co,dssco,cnstco,vco,masspp
     &                       ,iret,errmsg)
                        IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
                     END IF
                     
                     IF(.not.start_conf) WRITE(kprint,70120)
                  END IF
               END IF
               IF(thermos) THEN
                  IF(DABS(t-temph) .GT. dtemph .OR. MOD(ninner,scale)
     &                 .EQ. 0) THEN
                     CALL set_tempt(neta,qmass,etap,temph,t)
                     WRITE(kprint,70130)
                  END IF
               END IF
               
               ig=0
            END IF

