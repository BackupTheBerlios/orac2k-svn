
*-------------------- VARIABLES IN COMMON ------------------------------

      INTEGER maxstp,nsave,nrject,nconf,nprint,nplot,nprop,nform,icl,icm
     &     ,icn,nmol,nhskip,nupdte,nupdte_index,nunit,llacc,lldon,nbone
     &     ,ntap,lbond,lbend,ltors,litor,int14p,nspw,maxrun,nflag,scale
     &     ,nheating,kgyr,ngyr,int13p,dunitgi,dunfrci,dunitgr,dunfrcr
     &     ,duitgi,dufrci,ngrdon,ntop,lstretch,lstrtch,lconstr,lcnstr
     &     ,nxrms,nvi,ndipole,nnative,nrespa,n1respa,n0respa,mrespa
     &     ,lrespa,ntmtss,ntmtsp,ncofactor,nfragm,fragm(2,m1)
     &     ,gofr_nprint,gofr_navg,gofr_ncomp,sofk_nprint,sofk_navg
     &     ,sofk_ncomp,nmol_create,mim_lim,ndf_thermos,neta,navg_str
     &     ,nrms,navg_str_xrms,atoms_diss(2,2),mol_diss(2),diss_list
     &     ,ncavities,ncontact_slt,contact_slt(2,pcontact_slt)
     &     ,start_anl,stop_anl,update_anl,top_bonds(ntopol+1)
     &     ,top_bendings(ntopol+1),top_ptors(ntopol+1),top_itors(ntopol
     &     +1),ntop_print,nfragm_dist,nhbonds,nhhisto,nabmd,corr_atoms
     &     ,nbinew,nstep_freq,nx_cav,ny_cav,nz_cav,nfreq_ef,nmol_ef
     &     ,pmol_ef,nfreq_dp,nxslt,ncalc_cofm,nlda_mol,nlda_atm
     &     ,nlda_zero,ninst_fit,ninst_lda,natms,patms,nlda_flu,nlda_hyd
     &     ,nprot_hyd,min_hyd,max_hyd,min_lda,max_lda,nlda_grp,nprot_lda
     &     ,Virtual_types,virtual_atoms


      INTEGER lbnd_x,lconstr_x,lbndg_x,ltor_x,litr_x,int14_x,int13_x
     &     ,ingrp_x

      INTEGER protl,nprot,protlb,nprotb,nprot_charges,prot_charges
     &     ,clsthl,nclsth,nprot_fix,prot_fix

      CHARACTER*1 shell_pme,P_dyn_update_shell

      INTEGER lacc,ldon,nbtype,nres,nbun,mend,mback,lbnd,concta,concth
     &     ,lbndg,ltor,int14,type14,litr,nhtype,lpnbd,lphyd,pea,peb,pec
     &     ,nec,nascii,atres,ngrp,grppt,ingrp,ingrpp,mres,resg,atomg
     &     ,int13,bsitp,asitp,nbsitp,nasitp,nrigg,nrigat,wres,totvaf
     &     ,prtvaf,novaf,divide_spline,iseed,protg,groupp,atomp,iz,nv1
     &     ,nv2,nv3,nv4,nospt,ns1,ns2,ns3,ns4,niop,nios,ngprint
     &     ,itor_ptype,nplot_fragm,nplot_center,ncx,ncy,ncz,type_table
     &     ,abmd_tors_num,divide_records,atom_record,abmd_natp
     &     ,abmd_nat,atres_map1,atres_map2,Boltz_Group
      INTEGER nplmap,plmap,ntime_q
      INTEGER rgmap,rgmapc,nsct,rgnato,nfield,ffield,fexcl,nfatm,nexatm,
     &     mfield,nesad,neead,maxsim,rggrup,rggmap,menerg,
     &     kenerg,nhoov,anntom,ancorr,annslv,annpro,anpoint,
     &     nvol,nfft1,nfft2,nfft3,nfft3_start,nfft3_local,nfft2_start
     &     ,nfft2_local,pme_order,bfgs_m,nfreq_polar,NSecStruct
     &     ,SecPointer(m1*10),SecStructTotal,NFreqPotential

      REAL*8  time,time_q,t,dtemp,radius,aaxis,baxis,caxis,volume
     &     ,volumepr,alphal,rkcut,hrson,hrsoff,hanon,hanoff,rspon,mhoov
     &     ,whoov,rspoff,rspcut,rcut_LJ,rcut_EL,hrcut,hrcut_update,hacut
     &     ,wmtp,alf,bet,gam,phyd2,phyd1,dtemph,delrg,gofr_cut,sofk_cut
     &     ,sofk_delta,cut_field,alphaf,qmass,dtemppr,rcut_hb,acut_hb
     &     ,a2cut_hb,annealing_fact,TimeLimit,TimeStartRun,SigmaIon
     &     ,ChargeIon,SmoothFactor,r_min_virt(100),mesos_rho,mass_pfix
      REAL*8  rmol,co,oc,chrge,mass,potbe,potto,ptorj,potit,pnbd1,pnbd2
     &     ,pnbd3,pnbd4,ecc6,ecc12,ecc146,ecc1412,c6pw,c12pw,deps,potbo
     &     ,potbo_cnst,ett,emass,echrg,ecut,xecm,yecm,zecm,lmbdt,tlmbdb
     &     ,csecp,stret,epcut,omega,rspset,spring,pmass,uealfa,uemin
     &     ,uemax,xet,yet,zet,sjorg,ejorg,c6jorg,c12jorg,fudge,lj_fudge
     &     ,lj_fudgeb,antstar,anfact,radii,eps_energy,xpt0,ypt0,zpt0
     &     ,erfc_bin,hhisto_bin,erfc_arr,rpass,nutime,erf_arr_corr,delew
     &     ,rlew,ruew,hstep_freq,size_atom_cav,bin_size_cav
     &     ,rmax_size_cav,coeff_hyd,plrzbij,coeff_lda,polar_scale
     &     ,polar_counter,Ext_ef
      REAL*8  densurfpt,probrad
      REAL*8  xsct,ysct,zsct,rdsct,hdmass,wcoef,duitgr,dufrcr
      REAL*8  pext,masspr,wpr
      REAL*8  taut,taup,mu,compressibility
      REAL*8  rcuth,rcutl,rcutm,rtolh,rtoll,rtolm,rneil,rneim,rneih,kl
     &     ,klt,km,kmt,distmax,abmd_tors_dir,abmd_cryst_dir
     &     ,abmd_kvect,P_nei_fact,abmd_nat_theta,abmd_nat_dir
     &     ,native_dist,native_theta,bfgs_factr
      LOGICAL linser,protei,solven,hydbnd,cutoff,mixrww,mixrwp,clewld,
     &     stinit,debug,debug_rs,debug_ct,debug_st,debug_bt,debug_it
     &     ,debug_pt,tpgfil,prttpg,electr,elinit,sction,hmass,tpgwbn
     &     ,grpcut,rescm,lenerg,uenerg,amphi,bending,slvpdb,stoprun,pme
     &     ,erfc_spline,hhisto,erfc_spline_corr,ltest_times,md_respa
     &     ,linked_cell,avg_str,abmd,abmd_tors,abmd_cryst,abmd_unbias
     &     ,abmd_native,dissociate,diss_mol,diss_atoms,associate
     &     ,remove_momentum,avg_rms,diffusion,avg_ca,avg_he,avg_bc
     &     ,write_ff_pars,annealing,heating,change_cell,resize_cell
     &     ,erf_corr,write_grad,equi,P_dyn_update,native,check_native
     &     ,efield,lmol_ef,print_ef,print_dp,inst_fit,calc_cofm,latms
     &     ,sel_lda,avg_lda,lda_flu,lda_hyd,prot_hyd,res_time,prot_lda
     &     ,polar,EPotential,DoFreeEnergy,Virtual_residue,skip_step
     &     ,scan_traj,Old_dipoles,mesos,old_tpg
      LOGICAL lfield,prtatl,prthbl,prtnbl,prtcnl,prtbal,prtptl,prtitl
     &     ,prtbndl,prtseq,prtipw,prtn14,fold,scharge,UnCharge,pfix
     &     ,nonbnd,noneq,wrtgyr,iobdmp,dmprnd_i,dmprnd_o,hoover,thermos
     &     ,widom,duitgl,dufrcl,vacf,time_corr,not_time_corr,analys
     &     ,anxrms,anxrms_cell,anrms,anslv,anxca,anxbc,anxhe,anxal
     &     ,anprot,adihed,stretch,cpress,pressure,isostress,coupl_mol
     &     ,coupl_grp,coupl_atm,lberendsen,protvold,protvolc,pvexclud
     &     ,sgroup,recstrc,mdsim,minimize,frequencies,conj_grad,steepest
     &     ,l_bfgs_b,igr,prttopl,prttor,prtben,prtstr,stretch_heavy
     &     ,lslvint,template,spectra,clean,start_conf,energy_then_die
     &     ,lgr,vacf_mts,replicate,gofr_avg,gofr_neighbor,gofr_intra
     &     ,gofr,s_of_k,restart_old,read_co,landersen,adjust_cnstr
     &     ,readjust_cnstr,tpgprm_read,cavities
     &     ,occupy_space,fragm_dist,hbonds_res,hbonds_tot,hbonds_vor
     &     ,ascii_nocell,ascii_wsc,restart_read,restart_write
     &     ,SecStructure,AddTime,zero_charges
      CHARACTER*80 titld0,gerald,cgroup,restart_in,restart_out
     &     ,cavities_file,dmpfil_i,dmpfil_o,What_to_Do_Pol
      CHARACTER*7 beta,betb,alnbd,types_virt(100)
      CHARACTER*8 anatom
      CHARACTER*8 prsymb,Residue_virt
      CHARACTER*80 title0

*-------  New definition of the solvent array --------------------------

      INTEGER nbun_slv,nato_slv,ngrp_slv,lbond_slv,lbend_slv,ltors_slv
     &     ,litor_slv,nbone_slv,int14p_slv
     &     ,int13p_slv,llacc_slv,lldon_slv,mback_slv,nres_slv,mend_slv
     &     ,grppt_slv,int14_slv,type14_slv,int13_slv,lacc_slv,ldon_slv
     &     ,nhtype_slv,lbnd_slv,lbndg_slv,ltor_slv,litr_slv,atres_slv
     &     ,concta_slv,nbtype_slv,bsitp_slv,asitp_slv,nbsitp_slv
     &     ,nasitp_slv,nrigg_slv,nrigat_slv,icl_slv,icm_slv,icn_slv
     &     ,ss_point,ss_index

      REAL*8  chrge_slv,mass_slv,potbo_slv,potbe_slv,potto_slv
     &     ,potit_slv

      CHARACTER*7 beta_slv,betb_slv,alnbd_slv
      CHARACTER*8 slv_residue,slv_redef_type
      LOGICAL slv_exist,slt_exist,slv_generate,slv_create,slv_redef
     &     ,slv_read,adihed_slv,pdb_read,slt_create,slv_randomize
     &     ,slv_add

*-------------------- VARIABLES IN COMMON ------------------------------

      COMMON /ipropr/ nhskip,llacc,lldon,nbone,ntap,lbond,lstretch
     &     ,lconstr,lbend,ltors,litor,lpnbd,lphyd,int14p,nspw,nunit
     &     ,nplmap,nbun,ngrp,ingrpp,int13p,nrigg,nrigat,itor_ptype,iz
     &     ,mim_lim,ntime_q,ndf_thermos(3),neta,plmap(40),wres
     &     ,lacc(2,m1),ldon(2,m1),nbtype(m1),mend(nores),atres(2,nores)
     &     ,atres_map1(m1+nores),atres_map2(m1+nores),nres(m1,2)
     &     ,mback(m1),grppt(2,m11),boltz_Group(m11),lbnd(2,m9),lstrtch(2
     &     ,m9),lcnstr(2,m9),concta(m1,m10),concth(m1,m10),lbndg(3,m2)
     &     ,ltor(4,m3),int14(2,m3),type14(m3),litr(4,m4),nhtype(m6,m6)
     &     ,ingrp(2,m12),mres(2,nores),int13(2,m2),atomg(m1),bsitp(m13
     &     ,m14),asitp(m13,m14),nbsitp(m14),nasitp(m14),protg(m15)
     &     ,groupp(m11),atomp(m1),type_table(m6,m6),lbnd_x(m9)
     &     ,lconstr_x(m9),lbndg_x(m2),ltor_x(m3),litr_x(m4),int14_x(m3)
     &     ,int13_x(m2),ingrp_x(m12),virtual_atoms(m1),Virtual_types

      COMMON /rrntim/ time,time_q(10),t,dtemp,radius,aaxis,baxis,caxis
     &     ,alf,bet,dtemph,dtemppr,gam,volume,volumepr,rmol(3,8),co(3,3)
     &     ,oc(3,3),whoov,mhoov(hoov),mass_pfix,qmass(3)
     &     ,annealing_fact,TimeLimit,TimeStartRun,SigmaIon,ChargeIon
     &     ,SmoothFactor,r_min_virt,mesos_rho,P_nei_fact,nutime,nhoov
     &     ,titld0,gerald,dmpfil_i,dmpfil_o,What_to_Do_Pol,restart_in
     &     ,restart_out,cavities_file,linser,AddTime,replicate,protei
     &     ,solven,hydbnd,cutoff,mixrww,mixrwp,clewld,pme
     &     ,remove_momentum,erfc_spline,erfc_spline_corr,erf_corr,tpgfil
     &     ,hhisto,rescm,tpgwbn,prttpg,grpcut,iobdmp,dmprnd_i,dmprnd_o
     &     ,slvpdb,hoover,thermos,adihed,stretch,stretch_heavy,recstrc
     &     ,mdsim,minimize,frequencies,conj_grad,steepest,write_grad
     &     ,l_bfgs_b,md_respa,stoprun,ltest_times,clean,start_conf
     &     ,annealing,heating,Old_Dipoles,mesos,old_tpg,change_cell
     &     ,resize_cell,energy_then_die,title0,linked_cell,restart_old
     &     ,restart_read,restart_write,read_co,landersen,adjust_cnstr
     &     ,readjust_cnstr,tpgprm_read,equi,P_dyn_update,native
     &     ,check_native,SecStructure
      COMMON /rmolpr/ alphal,rkcut,delrg,rpass,rcuth,rcutl
     &     ,rcutm,rtolh,rtoll,rtolm,rneil,rneim,rneih,kl,klt
     &     ,km,kmt,distmax,widom,shell_pme,P_dyn_update_shell

      COMMON /rpropr/ hrson,hrsoff,hanon,hanoff,rspon,rspoff,rspcut
     &     ,rcut_LJ,rcut_EL,hrcut,hrcut_update,hacut,wmtp,eps_energy
     &     ,bfgs_factr,hstep_freq,size_atom_cav,bin_size_cav
     &     ,rmax_size_cav,erfc_bin,hhisto_bin,erfc_arr(4,mspline)
     &     ,erf_arr_corr(4,mspline),delew,rlew,ruew,chrge(m1),mass(m1)
     &     ,pmass(m1),wcoef,potbe(m2,7),potbo(m9,2),potbo_cnst(m9,2)
     &     ,potto(m3,2),potit(m4,3),pnbd1(m6),pnbd2(m6),pnbd3(m6)
     &     ,pnbd4(m6),ecc6(m5),ecc12(m5),ecc146(m5),ecc1412(m5),c6pw(m7)
     &     ,c12pw(m7),deps(m7),phyd2(m5),phyd1(m5),csecp(4,m13,m14)
     &     ,radii(m5),prsymb(n5),beta(m1),betb(m1),alnbd(m6),types_virt
     &     ,Residue_virt
      COMMON /irntim/ maxstp,nsave,nrject,nconf,nprint,nplot,nprop,
     &     nascii,nupdte,nupdte_index,maxrun,scale,nheating,nflag(10)
     &     ,nrespa,mrespa,n0respa,n1respa,lrespa,ntmtss(3),ntmtsp(5)
     &     ,ngprint,nplot_fragm,nplot_center,ncx,ncy,ncz,nstep_freq
     &     ,bfgs_m
      COMMON /imolpr/ nform,icl,icm,icn,nmol,nmol_create,debug,debug_rs
     &     ,debug_ct,debug_st,debug_bt,debug_it,debug_pt,stinit,ngrdon
     &     ,iseed,igr,ntop,nfft1,nfft2,nfft3,nfft3_start,nfft3_local
     &     ,nfft2_start,nfft2_local,pme_order,protl(m1*2),nprot
     &     ,protlb(m1*2),nprotb,clsthl(m1*2),nclsth,nprot_charges
     &     ,prot_charges(npm),nprot_fix,prot_fix(npm+1),nx_cav,ny_cav
     &     ,nz_cav
      COMMON /quantum/pea,peb,pec,nec,ett,emass,echrg,stret,ecut,xecm,
     &                yecm,zecm,xet(2),yet(2),zet(2),lmbdt,tlmbdb,
     &                epcut,electr,elinit
      COMMON /region/ hdmass,xsct(10),ysct(10),
     &                zsct(10),rdsct(10),hmass,sction,nsct
      COMMON /regmap/ rgmap(m1),rgmapc(m1),rggmap(m11),rggrup,rgnato
      COMMON /prttop/ prtatl,prthbl,prtnbl,prtcnl,prtbal,prtptl,prtbndl
     &     ,prtseq,prtitl,prtipw,prtn14,prttor,prtben,prtstr,prttopl
      COMMON /field/  ffield(f1,f2),fexcl(f3,f2),nfatm(f2),nexatm(f2),
     &                nfield,mfield,lfield
      COMMON /fold/   abmd_tors_dir,abmd_cryst_dir,abmd_nat_dir
     &     ,abmd_nat_theta,native_dist,native_theta,abmd_kvect,omega
     &     ,rspset,spring,atoms_diss,mol_diss,diss_list(2,m1)
     &     ,abmd_tors_num,nabmd,abmd_natp,abmd_nat(2,nores),kgyr,ngyr
     &     ,nnative,scharge,UnCharge,pfix,zero_charges,nonbnd,fold,abmd
     &     ,abmd_tors,abmd_cryst,dissociate,associate,diss_mol
     &     ,diss_atoms,abmd_unbias,wrtgyr,abmd_native
      COMMON /neq/    nesad(2),neead(2),maxsim,noneq
      COMMON /umeng/  uealfa,uemin,uemax,menerg,kenerg,lenerg,uenerg
      COMMON /jorgn/  fudge,lj_fudge,lj_fudgeb,sjorg(m7),ejorg(m7)
     &     ,c6jorg(m7),c12jorg(m7),ptorj(m3,4),resg(m11),amphi,bending
      COMMON /vafcom/  totvaf,prtvaf,novaf,vacf,vacf_mts,time_corr
     &     ,not_time_corr,nv1,nv2,nv3,nv4,nospt,ns1,ns2,ns3,ns4,nios
     &     ,niop,spectra,lgr
      COMMON /dummy/  duitgr(10),dufrcr(10),duitgi(10),dufrci(10),
     &                dunitgi,dunfrci,dunitgr,dunfrcr,duitgl,dufrcl
      COMMON /anlysis1/ antstar,anfact,xpt0(m1),ypt0(m1),zpt0(m1)
     &     ,rcut_hb,acut_hb,a2cut_hb,gofr_cut,sofk_cut,sofk_delta
     &     ,cut_field,alphaf,anntom,annslv,ancorr,annpro,start_anl
     &     ,stop_anl,update_anl,nxrms,nrms,navg_str,navg_str_xrms,nvi
     &     ,ncofactor,ntop_print,nfragm_dist,nbinew,nhbonds,nhhisto
     &     ,nfragm,fragm,ncavities,top_bonds,top_bendings,top_ptors
     &     ,top_itors,ncontact_slt,contact_slt,ndipole,gofr_nprint
     &     ,gofr_navg,gofr_ncomp,sofk_nprint,sofk_ncomp,divide_records
     &     ,atom_record,divide_spline,anpoint(2,ndan1),corr_atoms(1+m1)
     &     ,nxslt,anatom(ndan1),nfreq_ef,lmol_ef,nmol_ef,pmol_ef(f1
     &     ,ndan1),nfreq_dp,ncalc_cofm,nlda_mol,nlda_atm,nlda_zero
     &     ,ninst_fit,ninst_lda,natms,patms(2,ndan1),nlda_flu,nlda_hyd
     &     ,nprot_hyd,min_hyd,max_hyd,min_lda,max_lda,nlda_grp,nprot_lda
     &     ,coeff_lda,Ext_ef(3),plrzbij(m5),coeff_hyd,polar_scale
     &     ,polar_counter,nfreq_polar,NSecStruct,SecPointer
     &     ,SecStructTotal,NFreqPotential



      COMMON /anlysis2/ prot_lda,analys,gofr,s_of_k,gofr_avg
     &     ,gofr_neighbor,gofr_intra,anprot,anxrms,anxrms_cell,anrms
     &     ,anslv,anxca,anxbc,anxhe,anxal,avg_str,avg_ca,avg_he,avg_bc
     &     ,avg_rms,write_ff_pars,diffusion,template,cavities
     &     ,fragm_dist,hbonds_res
     &     ,hbonds_tot,hbonds_vor,ascii_nocell,ascii_wsc,polar,efield
     &     ,EPotential,print_ef,print_dp,latms,sel_lda,avg_lda,inst_fit
     &     ,calc_cofm,lda_flu,lda_hyd,prot_hyd,res_time,DoFreeEnergy
     &     ,Virtual_residue,skip_step,scan_traj

      COMMON /press/ pressure,cpress,pext,masspr,wpr,isostress,coupl_mol
     &     ,coupl_grp,coupl_atm
      COMMON /extbath/ taut,taup,mu,compressibility,lberendsen
      COMMON /surfvol/ protvold,protvolc,densurfpt,probrad,pvexclud,nvol
      COMMON /spaceg/  cgroup,sgroup

*-------  New definition of the solvent common -------------------------

      COMMON /slvarray1/ chrge_slv(slvatm),mass_slv(slvatm)
     &     ,potbo_slv(slv2,2),potbe_slv(slv3,4),potto_slv(slv4,4)
     &     ,potit_slv(slv5,4),ss_point(m1+1,2),ss_index(m1)

      COMMON /slvarray2/ nbun_slv,nato_slv,ngrp_slv,lbond_slv,lbend_slv
     &     ,ltors_slv,litor_slv,nbone_slv,int14p_slv
     &     ,int13p_slv,llacc_slv,lldon_slv,mback_slv(slvatm)
     &     ,nres_slv(slvatm,2),mend_slv(slvunit),grppt_slv(2,slv_group)
     &     ,int14_slv(2,slv4),type14_slv(slv4),int13_slv(2,slv3)
     &     ,lacc_slv(2,slvatm),ldon_slv(slvatm),nhtype_slv(types,types)
     &     ,lbnd_slv(2,slv2),lbndg_slv(3,slv3),ltor_slv(4,slv4)
     &     ,litr_slv(4,slv5),atres_slv(2,slvunit),concta_slv(slvatm,slv6
     &     ),nbtype_slv(slvatm),bsitp_slv(slv_cns1,slv_cns2)
     &     ,asitp_slv(slv_cns1,slv_cns2),nbsitp_slv(slv_cns2)
     &     ,nasitp_slv(slv_cns2),nrigg_slv,nrigat_slv,icl_slv,icm_slv
     &     ,icn_slv

      COMMON /slvarray3/ slv_exist,slt_exist,slv_generate,slv_create
     &     ,slv_redef,slv_read,slv_randomize,pdb_read,adihed_slv
     &     ,slt_create,slv_add,beta_slv(slvatm),betb_slv(slvatm)
     &     ,alnbd_slv(slvatm),slv_residue,slv_redef_type

*-------------------- VARIABLES IN COMMONMS END ------------------------

