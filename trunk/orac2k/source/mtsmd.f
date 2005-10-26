      SUBROUTINE mtsmd(mapnl,xp0,yp0,zp0,xpg,ypg,zpg,eta,xpcm,ypcm,zpcm)

************************************************************************
*                                                                      *
*     MTSMD is the driver of the MD run when multiple time scales      *
*     methods are used. The integration algorithm is r-RESPA           *
*     of order o(dt^3) (Tuckermann et al. JCP 97 1990 (1992))          *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     Written by Piero Procacci, CECAM-ENS Lyon 1995                   *
*     New constant pressure and temperature algorithms implemented     *
*     by Massimo Marchi                                                *
*                                                                      *
*______________________________________________________________________*
*                                                                      *
*  MTSMD externals:                                                    *
*     adjust_bonds analyse_voronoi anneal appbou asng_xrmsw            *  
*     bcnstp calc_avg_str calc_avg_xrms calc_gofr calc_xrms            *
*     change_coord_inner change_frame change_origin check_topology     *
*     chkewl collision comp_abmd_fdiss comp_dip comp_fcm               *
*     comp_forcep comp_molmass comp_neigh_vor comp_stress_conf         *
*     comp_stress_kinetic comp_thermos_energy comp_thermos_forces      *
*     comp_vcm comp_vel_labframe comp_voronoi copy_protl correc        *
*     correc_etap correc_exp_scale correc_matr correc_stress cov_thermo*
*     cself daxpy dcopy dscal dumprs                                   *
*     erfc_spline_init erf_corr_cutoff ferrf fft_pme_init              *
*     find_diss_mol find_igint_vor find_thermos fndgrp gen_abmd_kvect  *
*     gen_stress_tag_bnd get_prot_cnstr get_type_slv icopy             *
*     inicmp int_corr_erf_spline kinetic lc_index lc_list              *
*     matinv mts_forces mts_furier mts_intra_n0   mts_intra_n1         *
*     mts_plotp mts_plot_fragm mts_test plotc plotd                    *
*     plot_center prtacc prtat prtba prtbnd                            *
*     prtcn prtfrc prtit prtite prtpt prtsq                            *
*     ranvel rattle_correc rattle_correc_co rattle_verlet              *
*     rattle_verlet_co readrs setup_skin_shell set_const_co            *
*     set_tempp set_tempt starting_verlet timer trans_center           *
*     tr_inbox verlet verlet_free verlet_free_eta write_bends          *
*     write_bonds write_confc write_gofrp write_gofrw write_pot_bond   *
*     write_pot_nbond write_tors xerror zero zero0                     *
*     zero3x3 zeroa zero_gofr zero_voronoi                             *
* 								       *
************************************************************************


*======================= DECLARATIONS ==================================

      USE Module_Fourier, ONLY: Fourier_init=>Init
      USE Class_Gauss_Param, ONLY: Param_Init=>Init, Found_Boltz 
      USE Class_Gcharges, ONLY: Gauss_Charges, Gauss_Init=>Init
      USE Poisson_Boltzmann
      USE Module_Extra_Forces, ONLY: inp_1, inp_2, inp_3, pop1, pop2,
     &     pop3, Extra_Force, Extra_Init=> Stretch_Init

      USE Module_Stress, ONLY: FixedAngles_Stress
      IMPLICIT none
#if defined PARALLEL
      INTERFACE
         SUBROUTINE P_Reduce_Forces(x,y,z,nstart,nend,nlocal,node,nprocs
     &        )
         REAL(8) :: x(*),y(*),z(*)
         INTEGER, OPTIONAL :: nstart,nend,nlocal,node,nprocs
         END SUBROUTINE P_Reduce_Forces
         SUBROUTINE P_Comm_Intra(nstart,nend,node,nprocs,x,y,z,xc,yc,zc
     &        ,nato,atomp,lstrtch,lstretch,lbnd_x,lbndg,lbend,lbndg_x
     &        ,ltor,ltors,ltor_x,litr,litor,litr_x,int14,int14p,int14_x
     &        ,int13,int13p,int13_x,ingrp,ingrpp,ingrp_x)


         INTEGER :: nstart,nend,node,nprocs
         REAL(8) :: x(*),y(*),z(*),xc(*),yc(*),zc(*)
         INTEGER, OPTIONAL :: nato,atomp(*),lstrtch(2,*),lstretch
     &        ,lbnd_x(*),lbndg(3,*),lbend,lbndg_x(*),ltor(4,*),ltors
     &        ,ltor_x(*),litr(4,*),litor,litr_x(*),int14(2,*),int14p
     &        ,int14_x(*),int13(2,*),int13p,int13_x(*),ingrp(2,*),ingrpp
     &        ,ingrp_x(*)
         END SUBROUTINE P_Comm_Intra
         SUBROUTINE P_Change_Decomposition(Decmp_name,nato,vpx,vpy,vpz
     &        ,nstart,nend,nstart1,nend1,node,nprocs)
         REAL(8) :: vpx(*),vpy(*),vpz(*)
         INTEGER :: nstart,nend,nstart1,nend1
         INTEGER :: node,nprocs,nato
         CHARACTER(*) :: Decmp_name
         END SUBROUTINE P_Change_Decomposition
         
      END INTERFACE
#endif
      include 'parst.h'

      INTEGER ma,mb,mh,ms,mf,memp,t1,mbs,numpr
      PARAMETER (ma=tsites,mb=tsitep,mh=elsiz,ms=secsiz,memp=1
     &     ,mf=mb*4,t1=types*(types+1)/2)

      PARAMETER (mbs=2*mb, numpr=npm)

*----------------------- ARGUMENTS -------------------------------------

      INTEGER mapnl(*)
      REAL*8  xp0(*),yp0(*),zp0(*),xpg(*),ypg(*),zpg(*),eta(*),xpcm(*)
     &     ,ypcm(*),zpcm(*)

*-------------------- VARIABLES IN COMMONS -----------------------------

      INCLUDE 'cpropar.h'
      INCLUDE 'fourier.h'
      INCLUDE 'pme.h'
      INCLUDE 'lc_list.h'
      INCLUDE 'unit.h'
      INCLUDE 'parallel.h'

*-------------------- DEFINITION OF AN EXTERNAL FUNCTION ---------------

      EXTERNAL  near0
      LOGICAL near0

*-------------------- LOCAL VARIABLES ----------------------------------

      CHARACTER*80 errmsg
      CHARACTER*7  eora
      character*1  rshell,rshk

      INTEGER    l2
      INTEGER    kp,kt,hours,min,nato_slt,iter_avg
      REAL*8  uptot,upstot,ustot,ektot,E0
      REAL*8  dentr,a_free,ts_free,dfree

      INTEGER ninner,nscal,ntot_fragm,fragm_1,fragm_2,ninn1,ninn0
      INTEGER i,j,nsstt,nscale,ncconf,tstep,nstep,mstep,iret,ka
      INTEGER ig,flag,il,im,in0,in1
      LOGICAL elflag,linit,lfirst,lprint
     &     ,lreturn,lupdate,lfalse,lrject,Poisson_Boltz
      REAL(8) :: virsp_h=0.0D0,virs_h=0.0D0,virs_m=0.0D0,virp_h=0.0D0
     &     ,ucns_h=0.0D0,ucns_l=0.0D0,ucns_m=0.0D0,ucos_h=0.0D0,ucos_l
     &     =0.0D0,ucos_m=0.0D0,ucnsp_h=0.0D0,ucnsp_l=0.0D0,ucnsp_m
     &     =0.0D0,ucosp_h=0.0D0,ucosp_l=0.0D0,ucosp_m=0.0D0,virsp_l
     &     =0.0D0,virs_l=0.0D0,virp_l=0.0D0,ucnp_h=0.0D0,ucnp_l=0.0D0
     &     ,ucnp_m=0.0D0,ucop_h=0.0D0,ucop_l=0.0D0,ucop_m=0.0D0,tl,tm
     &     ,tn0,tn1,tm2,tm4,tl2,time2,tn12,tn02,virsp_m,de_remove,fudgec
     &     ,fstep,dnit,time_q1,time_q2
      REAl(8) ::  urcs_h=0.0D0,urcs_l=0.0D0,urcs_m=0.0D0,urcp_h=0.0D0
     &     ,urcp_l=0.0D0,urcp_m=0.0D0,urcsp_h=0.0D0,urcsp_l=0.0D0
     &     ,urcsp_m=0.0D0,virp_m=0.0D0,eer=0.0D0,eer_h=0.0D0,eer_l=0.0D0
     &     ,eer_m=0.0D0,stressr_h(3,3),stressr_l(3,3),stressr_m(3,3)
     &     ,stressr_n1(3,3),stressr_n0(3,3),stressd_h(3,3),stressd_l(3,3
     &     ),stressd_m(3,3),stressd_n1(3,3),stressd_n0(3,3),stressd_p(3
     &     ,3),stressd(3,3),stressr(3,3),stress_conf(3,3),stress_kin(3,3
     &     ),stress_tot(3,3),prt_m(3,3),prt_l(3,3),prt_h(3,3),prt_n1(3,3
     &     ),prt_n0(3,3),st_m(3,3),gmgp(3,3),press_m,press_l,press_h
     &     ,press_n1,press_n0,press_conf,press_kin
      REAL(8) :: stressc(3,3)

      REAL*8  temp,fcpu,elapse,puconf,pucoul,puhyd,pubnd
     &     ,ubend,uptors,uitors,uconf,ucoul,ureal,urecp,ucek,pucek,tempt
     &     ,fsrtal,fsrtalc,tempr,temppr,tcm,rcm,elaps,fnstep,unbond
     &     ,cnbond,unb14,cnb14,timeq,xl,xlcut,sftalp,fsbond,fsbend
     &     ,fsin14,gsbond,gsbend,purecp,vfcp,tfcp,elps,etime,etimeq
     &     ,ttime,ttimeq,uumb,ungrp,cngrp,upconf,upcoul,timef,timefq
     &     ,gcpu,temph,uceh,hpot,ubond,aux,uslvbon,uslvben,gcpu_hd
     &     ,gcpu_ld,gcpu_md,gcpu_1,gcpu_2,gcpu_3,theoric_speed_up,ucepr
     &     ,gcpu_u,gcpu_0,gsin14,dips(3),vol_gofr,qt(4),gcpu_inn
     &     ,vfcp_inn,tfcp_inn,tdelta_inn,gcpu_hh,vfcp_hh,tfcp_hh
     &     ,tdelta_hh,gcpu_ll,vfcp_ll,tfcp_ll,tdelta_ll,gcpu_mm,vfcp_mm
     &     ,tfcp_mm,tdelta_mm,gcpu_rt,gcpu_sh,ngcpu_rt,ngcpu_sh
      REAL*8 CTime_h,GTime_h,CTime_l,GTime_l,CTime_m,GTime_m
      REAL*8  gr,pueng,gra
      INTEGER navg
      PARAMETER (navg = 60)

      REAL*8  sumarray(navg),ssmarray(navg)
      REAL*8 sum_econf,sum_ecoul,sum_enbnd,sum_etotpot,sum_tote,sum_ucek
     &     ,sum_temp,sum_tempt,sum_tempr,sum_temppr,sum_temph
     &     ,sum_pecek,sum_pehyd,sum_peconf,sum_pecoul,sum_percip
     &     ,sum_enb14,sum_ebend,sum_ebond,sum_eitor,sum_eptor,sum_pnbd
     &     ,sum_pebnd,sum_pepot,sum_ptote,sum_gr,sum_epcoul,sum_epconf
     &     ,sum_presst,sum_press,sum_pressc,sum_pressk,sum_st(3,3)
     &     ,sum_co(3,3),sum_volume,sum_pv,sum_temppra,sum_eslvint
     &     ,sum_eebond,sum_eebend,sum_eeptors,sum_eeitors

      REAL*8 ssm_econf,ssm_ecoul,ssm_enbnd,ssm_etotpot,ssm_tote,ssm_ucek
     &     ,ssm_temp,ssm_tempt,ssm_tempr,ssm_temppr,ssm_temph
     &     ,ssm_pecek,ssm_pehyd,ssm_peconf,ssm_pecoul,ssm_percip
     &     ,ssm_enb14,ssm_ebend,ssm_ebond,ssm_eitor,ssm_eptor,ssm_pnbd
     &     ,ssm_pebnd,ssm_pepot,ssm_ptote,ssm_gr,ssm_epcoul,ssm_epconf
     &     ,ssm_presst,ssm_press,ssm_pressc,ssm_pressk,ssm_co(3,3)
     &     ,ssm_st(3,3),ssm_volume,ssm_pv,ssm_temppra,ssm_eslvint
     &     ,ssm_eebond,ssm_eebend,ssm_eeptors,ssm_eeitors

      EQUIVALENCE      
     &     (sum_econf, sumarray(1)), (sum_ecoul, sumarray(2)),
     &     (sum_enbnd, sumarray(3)), (sum_etotpot, sumarray(4)),
     &     (sum_tote, sumarray(5)), (sum_ucek, sumarray(6)),
     &     (sum_temp, sumarray(7)), (sum_tempt, sumarray(8)),
     &     (sum_tempr, sumarray(9)), (sum_temppr, sumarray(10)),
     &     (sum_temph, sumarray(11)),
     &     (sum_pecek, sumarray(13)), (sum_pehyd, sumarray(14)),
     &     (sum_peconf, sumarray(15)), (sum_pecoul, sumarray(16)),
     &     (sum_enb14, sumarray(17)), (sum_ebend, sumarray(18)),
     &     (sum_ebond, sumarray(19)), (sum_eitor, sumarray(20)),
     &     (sum_eptor, sumarray(21)), (sum_pnbd, sumarray(22)),
     &     (sum_pebnd, sumarray(23)), (sum_pepot, sumarray(24)),
     &     (sum_ptote, sumarray(25)), (sum_gr, sumarray(26)),
     &     (sum_epcoul , sumarray(27)), (sum_epconf, sumarray(28)),
     &     (sum_st(1,1), sumarray(29)), (sum_co(1,1), sumarray(38)),
     &     (sum_presst, sumarray(47)), (sum_press, sumarray(48)),
     &     (sum_pressc, sumarray(49)), (sum_pressk, sumarray(50)),
     &     (sum_volume, sumarray(51)), (sum_pv, sumarray(52)),
     &     (sum_temppra,sumarray(53)), (sum_eslvint,sumarray(54)),
     &     (sum_percip,sumarray(55)),  (sum_eebond, sumarray(56)),
     &     (sum_eebend, sumarray(57)),  (sum_eeptors, sumarray(58)),
     &     (sum_eeitors, sumarray(59))
      EQUIVALENCE      
     &     (ssm_econf, ssmarray(1)), (ssm_ecoul, ssmarray(2)),
     &     (ssm_enbnd, ssmarray(3)), (ssm_etotpot, ssmarray(4)),
     &     (ssm_tote, ssmarray(5)), (ssm_ucek, ssmarray(6)),
     &     (ssm_temp, ssmarray(7)), (ssm_tempt, ssmarray(8)),
     &     (ssm_tempr, ssmarray(9)), (ssm_temppr, ssmarray(10)),
     &     (ssm_temph, ssmarray(11)), 
     &     (ssm_pecek, ssmarray(13)), (ssm_pehyd, ssmarray(14)),
     &     (ssm_peconf, ssmarray(15)), (ssm_pecoul, ssmarray(16)),
     &     (ssm_enb14, ssmarray(17)), (ssm_ebend, ssmarray(18)),
     &     (ssm_ebond, ssmarray(19)), (ssm_eitor, ssmarray(20)),
     &     (ssm_eptor, ssmarray(21)), (ssm_pnbd, ssmarray(22)),
     &     (ssm_pebnd, ssmarray(23)), (ssm_pepot, ssmarray(24)),
     &     (ssm_ptote, ssmarray(25)), (ssm_gr, ssmarray(26)),
     &     (ssm_epcoul , ssmarray(27)), (ssm_epconf, ssmarray(28)),
     &     (ssm_st(1,1), ssmarray(29)),(ssm_co(1,1),ssmarray(38)),
     &     (ssm_presst, ssmarray(47)), (ssm_press, ssmarray(48)),
     &     (ssm_pressc, ssmarray(49)), (ssm_pressk, ssmarray(50)),
     &     (ssm_volume, ssmarray(51)), (ssm_pv, ssmarray(52)),
     &     (ssm_temppra,ssmarray(53)), (ssm_eslvint,ssmarray(54)),
     &     (ssm_percip,ssmarray(55)),  (ssm_eebond, ssmarray(56)),
     &     (ssm_eebend, ssmarray(57)),  (ssm_eeptors, ssmarray(58)),
     &     (ssm_eeitors, ssmarray(59))

*=======================================================================
*    Variables needed to compute the stress tensor and implement       =
*    constant-stress simulations                                       =
*=======================================================================

      REAL*8  pressc
      REAL*8  stresstk(3,3)
      REAL*8  temppra,energy,time_fs

*----------- LOCAL WORK ARRAYS FOR THE RUN -----------------------------

      INTEGER mrject,ngrp_old,nprot_old,nind,indxyz,ind_a
      INTEGER, DIMENSION (:), ALLOCATABLE :: indxi,indxj,indxk

      REAL*8  max_dist(3)
      REAL(8), DIMENSION (:), POINTER ::  fpx_n0,fpy_n0,fpz_n0,fcax_n0
     &     ,fcay_n0,fcaz_n0
      REAL(8), DIMENSION (:), POINTER ::  fpx_n1,fpy_n1,fpz_n1,fcax_n1
     &     ,fcay_n1,fcaz_n1
      REAL(8), DIMENSION (:), POINTER ::  fpx_m,fpy_m,fpz_m,fcax_m
     &     ,fcay_m,fcaz_m
      REAL(8), DIMENSION (:), POINTER ::  fpx_l,fpy_l,fpz_l,fcax_l
     &     ,fcay_l,fcaz_l
      REAL(8), DIMENSION (:), POINTER ::  fpx_h,fpy_h,fpz_h,fcax_h
     &     ,fcay_h,fcaz_h
      REAL(8), DIMENSION (:), POINTER ::  phi

      REAL(8), DIMENSION (:), POINTER ::  fpx,fpy,fpz
      REAL*8 grad_max
      
      REAL*8  fth(3)
      INTEGER mapdn(2,mf),nmapdn(mf)
      REAL(8), DIMENSION (:), POINTER ::  vpx,vpy,vpz,vpx1,vpy1,vpz1
     &     ,vcax,vcay,vcaz,vcbx,vcby,vcbz

      TYPE(Gauss_Charges) :: Boltzmann

      REAL*8  etap(hoov),vh1(hoov),vco(3,3)

      INTEGER get_pointer_protl
      
*--- DYNAM is a scratch common block: here to save storage; not passed
*    to any of the external

      COMMON /dynam/ fpx_h,fpy_h,fpz_h,fpx_l,fpy_l,fpz_l,fpx_m,fpy_m
     &     ,fpz_m,fpx_n0,fpy_n0,fpz_n0,fpx_n1,fpy_n1,fpz_n1,fcax_m
     &     ,fcay_m,fcaz_m,fcax_l,fcay_l,fcaz_l,fcax_h,fcay_h,fcaz_h
     &     ,fcax_n0,fcay_n0,fcaz_n0,fcax_n1,fcay_n1,fcaz_n1,vpx,vpy,vpz
     &     ,vpx1,vpy1,vpz1,etap,vh1,vcax,vcay,vcaz,vcbx,vcby,vcbz,vco
     &     ,mapdn,tag_bndg,nmapdn

*     Phony forces and energies for neighbor list

      REAL*8  ucns_p,ucos_p,virs_p,virsp_p,ucnp_p,ucop_p,ucnsp_p,ucosp_p
     &     ,fpx_p(1),fpy_p(1),fpz_p(1),rcuth_save,rtolh_save,rcutl_save
     &     ,rtoll_save,conf_bnd_slt,coul_bnd_slt,conf_bnd_slv
     &     ,coul_bnd_slv,self_slt,self_slv,uslvtor,uslvitor,fscnstr_slt
     &     ,fscnstr_slv,conf_bnd_slt_n1,coul_bnd_slt_n1,conf_bnd_slv_n1
     &     ,coul_bnd_slv_n1
      
      LOGICAL dpress,dhoover,grflag,exist,mask(m1),h_skin,l_skin,m_skin
     &     ,n1_skin,n0_skin,Param_out
      INTEGER dthree,n_q
      REAL*8 zz1(3,3),lzz,lzz1
      INTEGER igrn,krdf(maxint*g1),worka(m1,m10),type_slv(slvatm)
      INTEGER numatoms,ntype_slv,offset_slv,nbetab_slv,abmd_cryst_nvect
      REAL*8 pol_type(m5),Ext(3)
      INTEGER tag_bndg(m2)
      PARAMETER(nbetab_slv=90,abmd_cryst_nvect=1)
      INTEGER itype_slv(nbetab_slv)
      CHARACTER*1 betab_slv(nbetab_slv)
      CHARACTER*80 tag,file_dbg,Polarization_Model(3)
      INTEGER mapnl0(m1),ingrpp0,ingrp0(2,3*m11),ingrp0_x(3*m11)
      REAL*8  dssp(mbs),coeffp(mbs)
      REAL*8  xpo(mb),ypo(mb),zpo(mb),coo(3,3),co2(3,3),oc2(3,3)
      REAL*8  xpa(mb),ypa(mb),zpa(mb),xpga(m11),ypga(m11),zpga(m11)
     &     ,xpcma(npm),ypcma(npm),zpcma(npm),lx(mb),ly(mb),lz(mb)
      REAL*8  work(mspline),wca(m1),whe(m1),wbc(m1),errca(npm),errhe(npm
     &     ),errbc(npm),erral(npm),drpca(m1),drpbc(m1),drphe(m1)
     &     ,drpal(m1),xp_avg(m1),yp_avg(m1),zp_avg(m1),tmass(numpr)
     &     ,tmassb(numpr),masspp(3),dssco(5),cnstco(2,5)
     &     ,abmd_cryst_vect(3,abmd_cryst_nvect),rtollo,yy1,yy2,yy
     &     ,ffwork(2),fact,temp_heat,dtemp_heat,U_conf,U_ele,Udirect
     &     ,Urecip,Uind,uself_dip,uself,Ugrp,U_Thole,Fixt_Dipoles(3,m1)
     &     ,U_solv

      REAL*8 dummy1,dummy3,TimeCurrent,virtual_energy,n_plus,n_minus
      INTEGER TimeToGo,TRemain
      INTEGER cnstp(2,mbs),cnstpp,offset,abmd_dir,cnstpp_slv
     &     ,cnst_protl(2*m1),cnst_protl_1(2*m1),cnst_protl_ex(2*m1)
     &     ,mapnl_save(m8),nmapnl,cnst_protp,cnst_protp_1,cnst_protp_ex
     &     ,count,tot_protl,tot_cnst,nmol_slt,nmol_slv,nmol_cm,natom_slt
     &     ,natom_slv,ntotal,npp,len,ierr,npp_u,npp_h
     &     ,npp_l,npp_m,M_get_length,ihplen,nbun_slt,counter


      COMMON /rag2/ temp,fcpu,elapse,puconf,pucoul,puhyd,pubnd,ubend
     &     ,uptors,uitors,uconf,ucoul,ureal,urecp,ucek,pucek,tempt
     &     ,fsrtal,tempr,temppr,tcm,rcm,fnstep,unbond,cnbond,unb14
     &     ,cnb14,timeq,xl,xlcut,sftalp,fsbond,fsbend,fsin14,gsbond
     &     ,gsbend,purecp,vfcp,tfcp,elps,etime,etimeq,ttime,ttimeq,uumb
     &     ,ungrp,cngrp,upconf,upcoul,timef,timefq,gcpu,temph,uceh,hpot
     &     ,ubond,aux,ucepr,temppra,dssp,coeffp,xpo,ypo,zpo,xpa,ypa,zpa
     &     ,xpga,ypga,zpga,xpcma,ypcma,zpcma,wca,whe,wbc,errca,errhe
     &     ,errbc,erral,drpca,drpbc,drphe,drpal,xp_avg,yp_avg,zp_avg
     &     ,tmass,tmassb,lx,ly,lz,work,krdf,cnstp,cnst_protl
     &     ,cnst_protl_1,cnst_protl_ex,mapnl_save
      DATA eora/'E.O.R. '/

*==================== EXECUTABLE STATEMENTS ============================



      lfalse = .FALSE.
      IF(.NOT.clean) THEN
         WRITE(kprint,20900)
      END IF   
*===  Check if the dimension of the work array are sufficient 


      IF(mb.LT.ntap) THEN
         errmsg=' While running MD: PARAMETER MB dimensions the work'
     &        //' arrays in sufficiently. Abort. '
         CALL xerror(errmsg,80,1,2)
         STOP
      END IF

      if (lrespa.EQ.0) lrespa=1
      if (mrespa.EQ.0) mrespa=1
      if (n1respa.Eq.0)n1respa=1
   
*=======================================================================
*----- Initialize some stuff -------------------------------------------
*=======================================================================

*===  set few variable to zero

      lreturn=.false.
      lprint=.false.
      flag=nflag(1)
      IF(flag .EQ. 1) THEN
         mrject=0
         nrject=0
      END IF
      IF(heating) THEN
         temp_heat=0.0D0
         dtemp_heat=t/DBLE(nheating)
      END IF
      mrject=nrject*lrespa*mrespa
      time_fs=time*DBLE(maxstp+nrject)
      nscale=0
      ncconf=0
      mstep=0
      temph=0.0D0
      temppra=0.0D0
      urcs_h=0.0D0
      urcs_l=0.0D0
      urcs_m=0.0D0
      urcp_h=0.0D0
      urcp_l=0.0D0
      urcp_m=0.0D0
      urcsp_h=0.0D0
      urcsp_l=0.0D0
      urcsp_m=0.0D0
      eer_h=0.0D0
      eer_l=0.0D0
      eer_m=0.0D0
      CALL zero3x3(stressr_h)
      CALL zero3x3(stressr_l)
      CALL zero3x3(stressr_m)
      CALL zero3x3(stressr_n1)
      CALL zero3x3(stressr_n0)
      CALL zero3x3(stressd_h)
      CALL zero3x3(stressd_l)
      CALL zero3x3(stressd_m)
      CALL zero3x3(stressd_n1)
      CALL zero3x3(stressd_n0)
      uumb=0.0D0
      conf_bnd_slv=0.0D0
      conf_bnd_slt=0.0D0
      tstep=0
      cnstpp=0
      cnstpp_slv=0
      cnst_protp=0
      fudgec=1.d0-fudge
      nsstt=0
      iret=0
      errmsg=' '
      nato_slt=ntap-nmol*nato_slv
      iter_avg=0
      gcpu_inn=0.0D0
      gcpu_mm=0.0D0
      gcpu_ll=0.0D0
      gcpu_hh=0.0D0
      gcpu_rt=0.0D0
      gcpu_sh=0.0D0
      ngcpu_rt=0
      ngcpu_sh=0
      CALL zero0(cpu_h,nprocs)
      ncpu_h=0
      npp_u=mpp
      npp_h=mpp
      npp_l=mpp
      npp_m=1
      indxyz=indmax

*=======================================================================
*--  Check if any of the particle is a poisson_boltzman charge      ----
*=======================================================================

      Poisson_Boltz=Found_Boltz(ntap,betb)
      IF(Poisson_Boltz) THEN
         CALL reset_nbond(ntap,betb,atomg,ss_index,Boltz_Group)
         n_plus=2.0D0/DSQRT(unitc)
         n_minus=2.0D0/DSQRT(unitc)
         Param_out=Param_Init(n_plus,n_minus,efact,300.0D0)
         Boltzmann=Gauss_Init(ntap,betb)
      END IF

*=======================================================================
*--  Set up NPT masses                                              ----
*=======================================================================

      IF(cpress) THEN
         DO i=1,3
            masspp(i)=masspr
         END DO
         IF(isostress .OR. FixedAngles_Stress) CALL set_const_co(co
     &        ,dssco,cnstco)
      END IF
      CALL comp_molmass(nprot,protl,mass,tmass)

*=======================================================================
*--  Set up some ABMD parameters                                    ----
*=======================================================================

      IF(abmd) THEN
         IF(fold) THEN
            IF(.NOT. abmd_unbias) THEN
               IF(rspset .LT. 0.0D0) THEN
                  abmd_dir=-1
               ELSE
                  abmd_dir=1
               END IF
               rspset=0.0D0
            END IF
         ELSE IF(dissociate .OR. associate) THEN
            IF(dissociate) THEN
               abmd_dir=1
            ELSE IF(associate) THEN
               abmd_dir=-1
            END IF
            rspset=0.0D0
            CALL find_diss_mol(diss_mol,diss_atoms,atoms_diss,mol_diss
     &           ,diss_list,nprot,protl)
         ELSE IF(abmd_cryst) THEN
            rspset=0.0D0
            abmd_dir=abmd_cryst_dir
            CALL gen_abmd_kvect(abmd_kvect,abmd_cryst_nvect
     &           ,abmd_cryst_vect,oc)
         ELSE IF(abmd_native) THEN
            IF(.NOT. abmd_unbias) rspset=0.0D0
            abmd_dir=abmd_nat_dir
            CALL get_atres(atres,nres(1,1),nato_slt,nbun_slt)
            CALL get_native(knative_tpl,nat_listp,nat_list,nores,atres
     &           ,beta,atres_map1,atres_map2,iret,errmsg)
            IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
         END IF
      END IF

*=======================================================================
*--- Initialize the velocities if at the beginning of the run        ---
*=======================================================================

      ALLOCATE(vpx(ntap),vpy(ntap),vpz(ntap))
      ALLOCATE(vpx1(ntap),vpy1(ntap),vpz1(ntap))
      ALLOCATE(vcax(ntap),vcay(ntap),vcaz(ntap))
      ALLOCATE(vcbx(ntap),vcby(ntap),vcbz(ntap))

      CALL zeroa(vpx,vpy,vpz,ntap,1)
      CALL zeroa(vpx1,vpy1,vpz1,ntap,1)
      CALL zeroa(vcax,vcay,vcaz,nprot,1)
      IF(avg_str) THEN
         CALL zeroa(xp_avg,yp_avg,zp_avg,ntap,1)
      END IF

*=======================================================================
*----- Initialize field arrays if needed 
*=======================================================================

      if(start_conf) then 
         max_dist(1) = distmax
         max_dist(2) = distmax
         max_dist(3) = distmax
      end if

      IF(nplot_fragm .GT. 0) THEN
         do i=1,nfragm
            IF ( fragm(1,i).gt.ntap.or.fragm(2,i).gt.ntap) THEN 
               WRITE(kprint,10260) i
10260          FORMAT(//  
     &              '**ERROR: fragments limits for',i5,' DEF_FRAGM',
     &              '(&READ_SOLUTE) exceeds protein atoms' )
               STOP
            END IF 
         end do
         ntot_fragm=0
         do i=1,nfragm
            fragm_1=fragm(1,i)
            fragm_2=fragm(2,i)
            ntot_fragm = ntot_fragm + fragm_2-fragm_1+1 
         end do    
      END IF

*=======================================================================
*----- Copy old CO and OC matrix to temporary arrays -------------------
*=======================================================================

      IF(change_cell .OR. resize_cell) THEN
         CALL dcopy(9,co,1,co2,1)
         CALL dcopy(9,oc,1,oc2,1)
      END IF

*=======================================================================
*----- Read input from restart file if needed --------------------------
*=======================================================================

      OPEN(unit=99,file='.WRITE_JOBTAG')
      WRITE(99,'(2i10)') 1000,0
      kp=0
      kt=-1
      IF(flag.GT.0) THEN
         elflag=electr.AND.(.NOT.elinit)
         WRITE(kprint,70300)
         IF(restart_old) THEN
            CALL readrs_old(kdump_in,restart_in,nstep,temp,ntap,ngrp_old
     &           ,nprot_old,xp0,yp0,zp0,vpx,vpy,vpz,xpg,ypg,zpg,xpcma
     &           ,ypcma,zpcma,vpx1,vpy1,vpz1,vcax,vcay,vcaz,eta,etap,vh1
     &           ,dthree,dhoover,sumarray,ssmarray,navg,navg,co,oc
     &           ,dpress,vco,zz1,lzz,lzz1,grflag,krdf,igrn,maxint*g1)
         ELSE
            CALL readrs(kdump_in,restart_in,nstep,temp,ntap,ngrp_old
     &           ,nprot_old,xp0,yp0,zp0,vpx,vpy,vpz,eta,etap,dthree
     &           ,dhoover,sumarray,ssmarray,navg,navg,co,oc,dpress,vco
     &           ,grflag,krdf,igrn,maxint*g1)
         END IF
            
         WRITE(kprint,70400)

         ninner=nstep*lrespa*mrespa
         ninn1=nstep*lrespa*mrespa*n1respa

         IF(flag .EQ. 1) THEN
            time_fs=time*DBLE(maxstp)
         END IF
      END IF
      
*=======================================================================
*---  set velocities to zero for fixed molecules
*=======================================================================

      IF(pfix) CALL Init_Fixed_Molecule(mass,ntap,vpx,vpy,vpz)
      
*=======================================================================
*---  set fake neighbor lists
*=======================================================================

      IF(AddTime) THEN
         maxstp=nstep+maxstp
         IF(flag .EQ. 1 .AND. nstep .GE. maxstp ) THEN
            WRITE(errmsg,'(''In MTSMD: Step number in restart = '',
     &           i6,'' while in input = '',i6,'' Increase TIME.'')')
     &           nstep,maxstp
            CALL xerror(errmsg,80,1,2)
         END IF
      ELSE
         IF(flag .EQ. 1 .AND. nstep .GE. maxstp ) THEN
            WRITE(errmsg,'(''In MTSMD: Step number in restart = '',
     &           i6,'' while in input = '',i6,'' Increase TIME.'')')
     &           nstep,maxstp
            CALL xerror(errmsg,80,1,2)
         END IF
      END IF
      REWIND(99)
      WRITE(99,'(2i10)') maxstp,maxrun
      CLOSE(unit=99)

*=======================================================================
*---  Write a Banner Page                                            ---
*=======================================================================

      WRITE(kprint,'(//)')
      WRITE(kprint,1200)
      WRITE(kprint,1300)
      WRITE(kprint,1100)
      WRITE(kprint,1100)
      WRITE(kprint,1300)
      WRITE(kprint,1200)
      WRITE(kprint,'(//)')

*=======================================================================
*---  Reset the averages when NFLAG = 0 or 2                         ---
*=======================================================================

      IF(flag.EQ.2. .OR. flag .EQ. 3 .OR.flag.EQ.0 
     &     .OR. flag .EQ. 4) THEN
          CALL zero(sumarray,ssmarray,navg,navg)
         nstep=0
         ninner=0
         ninn1=0
         ninn0=0
         a_free=0.0D0
         ts_free=0.0D0
      END IF

      IF(anxrms .OR. gofr .OR. avg_str) CALL asng_xrmsw(ss_point,m1+1
     &     ,wca,whe,wbc,beta,mback,nbone)

      IF(anxrms_cell) THEN
         anprot=.TRUE.
         annpro=annpro+1
         anpoint(1,annpro)=1
         anpoint(2,annpro)=ntap
      END IF


      IF(gofr) CALL get_type_slv(nato_slt,nato_slv,beta,betab_slv
     &     ,nbetab_slv,ntype_slv,type_slv,itype_slv,offset_slv
     &     ,types_gofr,iret,errmsg)
      IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)

*=======================================================================
*---- Initialize G of Rs -----------------------------------------------
*=======================================================================

      IF(gofr) THEN
         CALL zero_gofr(maxint,krdf,ngrdon,offset_slv)
      END IF

*=======================================================================
*-------- Compute two vectors used by the constraint subroutine --------
*=======================================================================

      IF(lconstr .NE. 0) THEN
         IF(adjust_cnstr) THEN
            WRITE(kprint,14000) 
            CALL adjust_bonds(ss_index,ntap,lcnstr,lconstr
     &        ,xp0,yp0,zp0,potbo_cnst,M9,iret,errmsg)
            IF(iret .EQ. 2) CALL xerror(errmsg,80,1,20)
         END IF
         cnstpp=0
         CALL bcnstp(ss_index,lcnstr,lconstr,potbo_cnst(1,2),mass,
     &        dssp,coeffp,cnstp,mbs,cnstpp,cnstpp_slv,iret,errmsg)
         IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
         IF(stretch_heavy) THEN
            CALL get_prot_cnstr(ntap,clsthl,nclsth,cnstp,cnstpp
     &           ,cnst_protl,cnst_protp,tot_protl,tot_cnst,mask,worka
     &           ,m1)
         ELSE
            CALL get_prot_cnstr(ntap,protl,nprot,cnstp,cnstpp,cnst_protl
     &           ,cnst_protp,tot_protl,tot_cnst,mask,worka,m1)
         END IF
      ELSE
      END IF
      virtual_energy=0.0D0
      IF(virtual_residue) THEN
         CALL get_atres(atres,nres(1,1),ntap,nbun)
         CALL virtual_extract(nbun,atres,beta,xp0,yp0,zp0,nres,m1
     &        ,prsymb,chrge,virtual_atoms)

         CALL Get_virtual_selfenergy(xp0,yp0,zp0,virtual_atoms,ecc6
     &        ,ecc12,nbtype,type_table,m6,virtual_energy)
      END IF

*=======================================================================
*--  Set up prot and protl vectors                                  ----
*=======================================================================

      CALL copy_protl(protl,protlb,nprot,nprotb)
      IF(coupl_atm .OR. coupl_grp) THEN
         CALL comp_molmass(nprot,protl,mass,tmassb)
         IF(lconstr .NE. 0) THEN
            nprot=nclsth
            CALL icopy(tot_protl,clsthl,1,protl,1)
         ELSE
            nprot=ntap
            count=0
            DO i=1,ntap
               protl(1+count)=1
               protl(2+count)=i
               count=count+2
            END DO
         END IF
         IF(nprot .GT. npm) THEN
            errmsg=
     &           'Number of molecules exceeds physical dimensions. '
     &           / /'Increase _N_PROT_MAX_ to '
            WRITE(errmsg(75:80),'(i6)') nprot
            CALL xerror(errmsg,80,1,2)            
         END IF
         CALL comp_molmass(nprot,protl,mass,tmass)
         CALL comp_nmol(ss_index,nprot,protl,nmol_slt,nmol_slv,natom_slt
     &        ,natom_slv)
      ELSE
         CALL comp_nmol(ss_index,nprot,protl,nmol_slt,nmol_slv,natom_slt
     &        ,natom_slv)
      END IF

*=======================================================================
*----- Setup Extra Stretching ------------------------------- ----------
*=======================================================================

      IF(Extra_Force) THEN
         IF(inp_1 % K /= 0.0D0) CALL Extra_Init(beta, nprot,protl,ntap
     &        ,inp_1,pop1)
         IF(inp_2 % K /= 0.0D0) CALL Extra_Init(beta, nprot,protl,ntap
     &        ,inp_2,pop2)
         IF(inp_3 % K /= 0.0D0) CALL Extra_Init(beta, nprot,protl,ntap
     &        ,inp_3,pop3)
      END IF

*=======================================================================
*----- Choose on which shell to carry out the gmgp correction ----------
*=======================================================================

      IF(cpress) THEN
         CALL setup_skin_shell(coupl_grp,coupl_mol,h_skin,l_skin,m_skin
     &        ,n1_skin,n0_skin)
      END IF

*=======================================================================
*----- Calculate solute center of mass coordinates and velocities ------
*=======================================================================

      IF(resize_cell) THEN
         CALL dcopy(9,co2,1,co,1)
         CALL dcopy(9,oc2,1,oc,1)
      END IF

      CALL inicmp(ss_index,xp0,yp0,zp0,xpcm,ypcm,zpcm,mass,nprot,protl)
      CALL change_frame(co,oc,-1,nprot,xpcm,ypcm,zpcm,xpcma,ypcma,zpcma)
      CALL comp_vcm(vpx,vpy,vpz,oc,nprot,protl,mass,tmass,vcax,vcay
     &        ,vcaz)
*=======================================================================
*----- Change the cell according to input commands ---------------------
*=======================================================================

      IF(change_cell) THEN
         CALL change_origin(1,nprot,protl,xp0,yp0,zp0,lx,ly,lz,xpcma
     &        ,ypcma,zpcma,co)

         CALL dcopy(9,co2,1,co,1)
         CALL dcopy(9,oc2,1,oc,1)

         CALL change_origin(-1,nprot,protl,xp0,yp0,zp0,lx,ly,lz,xpcma
     &        ,ypcma,zpcma,co)
         CALL change_origin(1,nprot,protl,xp0,yp0,zp0,lx,ly,lz,xpcma
     &        ,ypcma,zpcma,co)
         CALL change_frame(co,oc,1,nprot,xpcma,ypcma,zpcma,xpcm,ypcm
     &        ,zpcm)
      END IF

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

      IF(nflag(1) .NE. 0 .AND. ngrp_old .NE. ngrp ) THEN
         WRITE(kprint,98000) ngrp_old,ngrp
         IF(cpress) THEN
            CALL set_tempp(masspr,vco,temppra,0.0D0)
         END IF
         IF(thermos) THEN
            CALL set_tempt(neta,qmass,etap,temph,0.0D0)
            CALL zero0(eta,neta)
         END IF
      END IF
      IF(nflag(1) .NE. 0 .AND. nprot_old .NE. nprot ) THEN
         WRITE(kprint,99000) nprot_old,nprot
         IF(cpress) THEN
            CALL set_tempp(masspr,vco,temppra,0.0D0)
         END IF
         IF(thermos) THEN
            CALL set_tempt(neta,qmass,etap,temph,0.0D0)
            CALL zero0(eta,neta)
         END IF
      END IF
      IF(.NOT. dhoover .AND. thermos) THEN
         WRITE(kprint,81000) 
         CALL set_tempt(neta,qmass,etap,temph,0.0D0)
         CALL zero0(eta,neta)
      END IF
      IF(.NOT. dpress .AND. cpress) THEN
         WRITE(kprint,82000) 
         CALL set_tempp(masspr,vco,temppra,0.0D0)
      END IF

*=======================================================================
*===== print potential and equilibirum value  (by matteo) ==============
*=======================================================================

      if(equi .and. slt_exist) THEN
        gcpu = unite*avogad/1000.d0/4.184d0
        write(kequi,*) lbond, ' ===== BOND TABLE  '
        CALL equi_bond(kequi,gcpu,betb,lstrtch,lstretch,xp0,yp0,zp0,
     &                 potbo(1,1),potbo(1,2))
        write(kequi,*) lbond, ' ===== CONSTR TABLE'
        CALL equi_bond(kequi,gcpu,betb,lcnstr,lconstr,xp0,yp0,zp0,
     &                 potbo_cnst(1,1),potbo_cnst(1,2))
        write(kequi,*) lbend, ' ===== BEND TABLE  '
        CALL equi_bend(kequi,gcpu,betb,lbndg,lbend,xp0,yp0,zp0,
     &                 potbe(1,1),potbe(1,2),potbe(1,3),potbe(1,4))
        write(kequi,*) ltors,
     &              ' ===== TORSION TABLE '
        CALL equi_tors(kequi,gcpu,betb,ltor,ltors,xp0,yp0,zp0,
     &                 potto(1,1),potto(1,2))
        write(kequi,*) litor,
     &              ' ===== IMPROPER TABLE '
        CALL equi_impro(kequi,gcpu,betb,litr,litor,xp0,yp0,zp0,
     &                 potit(1,1),potit(1,2),potit(1,3))
      END IF


*=======================================================================
*---- Allocate memory for neighbor list                             ----
*=======================================================================

      IF(linked_cell) THEN
         indxyz=(2*(ncx-1)+1)*(2*(ncy-1)+1)*(2*(ncz-1)+1)
         ALLOCATE(indxi(indxyz),indxj(indxyz),indxk(indxyz))
      END IF
      npp_u=mpp


*=======================================================================
*----- Get for each node its atomic decomposition ----------------------
*=======================================================================

      CALL change_frame(co,oc,-1,ngrp,xpg,ypg,zpg,xpga,ypga,zpga)
#if defined PARALLEL
      WRITE(kprint,93410) nprocs
      WRITE(kprint,14010)
#endif
      CALL P_setup_decomp(node,nprocs,ncube,rbyte,nbyte,nstart_h
     &     ,nend_h,nlocal_h,nstart_ah,nend_ah,nlocal_ah,nstart_1,nend_1
     &     ,nlocal_1,nstart_ex0,nend_ex0,nlocal_ex0,nstart_ex,nend_ex
     &     ,nlocal_ex,nstart_2,nend_2,nlocal_2,nstart_cmi,nend_cmi
     &     ,nlocal_cmi,nstart_cm,nend_cm,nlocal_cm,nstart_cme,nend_cme
     &     ,nlocal_cme,nstart_cme0,nend_cme0,nlocal_cme0,nstart_g1
     &     ,nend_g1,nlocal_g1,nstart_gex,nend_gex,nlocal_gex,nstart_gex0
     &     ,nend_gex0,nlocal_gex0,ptr_ex,ptr_cme,ptr_gex,ntot_1,ntot_cmi
     &     ,ntot_ex,ntot_cme,nprot,protl,cnst_protp_1,cnst_protl_1
     &     ,cnst_protp_ex,cnst_protl_ex,cnstp,cnst_protp,cnst_protl
     &     ,lstrtch,lstretch,lbndg,lbend,cpu_h,ncpu_h,atomp,atomg
     &     ,ntap,ngrp,grppt,worka,xpa,ypa,zpa,xpga,ypga,zpga)

#if defined PARALLEL
      CALL P_Comm_Intra(nstart_1,nend_1,node,nprocs,xp0,yp0,zp0,xp0,yp0
     &     ,zp0,ntap,atomp,lstrtch,lstretch,lbnd_x,lbndg,lbend,lbndg_x
     &     ,ltor,ltors,ltor_x,litr,litor,litr_x,int14,int14p,int14_x
     &     ,int13,int13p,int13_x,ingrp,ingrpp,ingrp_x)

      IF(nprocs > 1) THEN
         npp_u=4*mpp/nprocs
      ELSE
         npp_u=mpp
      END IF
         
#else
      npp_u=mpp
#endif

*=======================================================================
*---- Divide mapnl -----------------------------------------------------
*=======================================================================

      CALL min_pack(ntap,mapnl,nmapnl)
      CALL icopy(nmapnl,mapnl,1,mapnl_save,1)
      CALL mapnl_divide(node,nstart_h,nend_h,grppt,mapnl)

*=======================================================================
*--- Initialize stress tag tables for bendings -------------------------
*=======================================================================

      CALL gen_stress_tag_bnd(lbend,3,lbndg,lbndg_x,atomp,tag_bndg)

*=======================================================================
*----- If the system is at constant temperature initialize some --------
*----- variables and convert masses of the 3 thermostats ---------------
*=======================================================================

      IF(thermos) THEN
         CALL find_thermos(ntap,lconstr,natom_slt,natom_slv,nmol_slt
     &        ,nmol_slv,cnstpp_slv,nprot,cpress,isostress,ndf_thermos)
     &        
         CALL cov_thermos(slv_exist,slt_exist,qmass,ndf_thermos,t)
         IF(flag .EQ. 0) THEN
            CALL zero0(eta,neta)
         END IF
      END IF

*=======================================================================
*-------- Initialize velocities when starting without a restart --------
*=======================================================================

      IF(FLAG.EQ.0.AND.t.GT.1d-15 .OR. FLAG .EQ. 3) THEN
         LINIT=.TRUE.
         
         IF(heating) THEN
            CALL ranvel(temp_heat,mass,ntap,vpx,vpy,vpz,xp0,yp0,zp0,co
     &           ,linit)
         ELSE
            CALL ranvel(t,mass,ntap,vpx,vpy,vpz,xp0,yp0,zp0,co,linit)
         END IF
         IF(cnstpp .NE. 0) THEN
            aux=0.0D0
            CALL rattle_correc(nstart_2,nend_2,time,xp0,yp0,zp0
     &           ,vpx,vpy,vpz,ntap,cnstp,dssp,coeffp,cnstpp,mass,dnit
     &           ,cnst_protp,cnst_protl,mim_lim,aux,iret,errmsg)
            IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
#if defined PARALLEL            
            CALL P_expand_r8x3(vpx,vpy,vpz,nstart_2,nend_2,nlocal_2,node
     &           ,nprocs,1)
            CALL P_expand_r8x3(xp0,yp0,zp0,nstart_2,nend_2,nlocal_2,node
     &           ,nprocs,1)
#endif
         END IF
         CALL comp_vcm(vpx,vpy,vpz,oc,nprot,protl,mass,tmass,vcax,vcay
     &        ,vcaz)
         LINIT=.FALSE.

*=======================================================================
*----- Set to zero velocities of the barostat --------------------------
*=======================================================================

         CALL zero0(temppra,1)
         CALL zero0(temph,1)
         IF(cpress) THEN
            CALL set_tempp(masspr,vco,temppra,0.0D0)
         END IF
         IF(thermos .AND. FLAG .NE. 3) THEN
            CALL set_tempt(neta,qmass,etap,temph,0.0D0)
         ELSE IF(thermos .AND. FLAG .EQ. 3) THEN
            CALL set_tempt(neta,qmass,etap,temph,t)            
         END IF
      END IF

*=======================================================================
*--- Compute the volume of the system ----------------------------------
*=======================================================================

      CALL matinv(3,3,co,oc,volume)
      volume=volume*boxl**3
     
*========================================================================
*==== Calls init routine for conventional kspace Ewald or PME -----------
*========================================================================

      IF(clewld) THEN
         IF(pme) THEN
            numatoms=ntap
            CALL fft_pme_init(numatoms,nfft1,nfft2,nfft3,pme_order
     &           ,sizfftab,sizffwrk,siztheta,siz_Q,sizheap,sizstack
     &           ,bsp_mod1,bsp_mod2,bsp_mod3,fftable,ffwork)
            if ( siz_Q .GT. MAXT ) THEN
               write(kprint,78410)
               stop
            END IF
            rshk=shell_pme
            IF(erf_corr) THEN
               CALL erf_corr_cutoff(oc,delew,rkcut,nfft1,nfft2,nfft3)
               CALL int_corr_erf_spline(rlew,ruew,nbinew,alphal,rkcut
     &              ,erf_arr_corr,work)
            END IF
            CALL Pme_init(node,nprocs,nodex,nodey,nodez,npy,npz,ictxt
     &           ,descQ,fftable,nfft1,nfft2,nfft3,nfft3_start
     &           ,nfft3_local,nfft2_start,nfft2_local,iret,errmsg)
            IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
         ELSE
            aux=(rtolh+rcuth)**2
            CALL chkewl(oc,aux,alphal,rkcut,volume)
            rshk=shell_pme
            IF(erf_corr) THEN  
              CALL int_corr_erf_spline(rlew,ruew,nbinew,alphal,rkcut
     &              ,erf_arr_corr,work)
            END IF
         END IF

*=======================================================================
*----- Set up table for erfc spline ------------------------------------
*=======================================================================

         IF(erfc_spline) THEN
            aux=rcuth+rtolh
            aux=aux*aux
            CALL erfc_spline_init(aux,alphal,erfc_bin,mspline,erfc_arr
     &           ,work,iret,errmsg,rkcut,erfc_spline_corr)
            IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
         END IF
      END IF
#if defined PARALLEL
      CALL P_Reduce_Forces(xp0,yp0,zp0,nstart_2,nend_2,nlocal_2,node
     &     ,nprocs)
#endif
*=======================================================================
*---- Calculate ewald self and the intramolecular terms. ---------------
*---- Only the contribution to the latter coming from  -----------------
*---- bond constraints is computed at this point -----------------------
*=======================================================================

      if(polar) then
         fact = polar_scale
         CALL dscal(ntap,fact,chrge,1)
         mapnl(1:m1)=0
         CALL igmap(ngrp,grppt,ingrpp0,ingrp0,3*m11,mapnl0,errmsg,iret)
         IF(What_to_do_Pol .EQ. ' ') THEN
            What_to_do_Pol='Full'
         END IF
         CALL Make_Pol_Models(What_to_Do_Pol,Polarization_Model)
         IF(What_To_Do_Pol .EQ. 'Full') WRITE(*,88000)
         IF(What_To_Do_Pol .EQ. 'Gauss') WRITE(*,87000)
      end if

      IF(clewld) THEN
         IF(pme .AND. shell_pme .NE. '0') THEN 
            CALL cself(ss_index,ntap,alphal,rkcut,chrge,self_slt
     &           ,self_slv)
            CALL ferrf(ss_index,alphal,chrge,1.0D0,xp0,yp0,zp0,0
     &           ,lcnstr,lconstr,lconstr_x,fscnstr_slt,fscnstr_slv,fpx_p
     &           ,fpy_p,fpz_p,erf_corr,erf_arr_corr,delew,rlew)
         ELSE IF(.NOT. pme) THEN
            CALL cself(ss_index,ntap,alphal,rkcut,chrge,self_slt
     &           ,self_slv)
            CALL ferrf(ss_index,alphal,chrge,1.0D0,xp0,yp0,zp0,0
     &           ,lcnstr,lconstr,lconstr_x,fscnstr_slt,fscnstr_slv,fpx_p
     &           ,fpy_p,fpz_p,erf_corr,erf_arr_corr,delew,rlew)
         END IF
#if defined PARALLEL
         IF(nprocs .GT. 1) THEN
            CALL P_merge_r8(fscnstr_slt)
            CALL P_merge_r8(fscnstr_slv)
         END IF
#endif
         IF(polar) THEN
            self_slv=0.0D0
            self_slt=0.0D0
c$$$            fscnstr_slt=0.0D0
c$$$            fscnstr_slv=0.0D0
         END IF
      END IF


c$$$*=======================================================================
c$$$*----- EPotential ------------------------------------------------------
c$$$*=======================================================================
c$$$
c$$$      IF(EPotential) THEN
c$$$         CALL CompElecPotentialOnGrid(co,xp0,yp0,zp0,pmass,ss_index
c$$$     &        ,mass,nfft1,nfft2,nfft3,nprocs,node,ncube,rbyte,nbyte,ntap
c$$$     &        ,ngrp,grppt,chrge,nstart_h,nend_h,nlocal_h,nstart_ah
c$$$     &        ,nend_ah,nlocal_ah,nprot,protl,pme_order,alphal,bsp_mod1
c$$$     &        ,bsp_mod2,bsp_mod3,rkcut) 
c$$$         WRITE(*,*) 'Success!'
c$$$      END IF
c$$$
*=======================================================================
*---  Check if there are no error in the list of bonds, bendings -------
*---  and torsions -----------------------------------------------------
*=======================================================================

      IF(prttopl) THEN
         CALL check_topology('B',top_bonds,lbond,iret,errmsg)
         IF(iret .NE. 0) CALL xerror(errmsg,80,1,2)
         CALL check_topology('D',top_bendings,lbend,iret,errmsg)
         IF(iret .NE. 0) CALL xerror(errmsg,80,1,2)
         CALL check_topology('P',top_ptors,ltors,iret,errmsg)
         IF(iret .NE. 0) CALL xerror(errmsg,80,1,2)
         CALL check_topology('I',top_itors,litor,iret,errmsg)
         IF(iret .NE. 0) CALL xerror(errmsg,80,1,2)
      END IF


*=======================================================================
*----- Print titles for the run ----------------------------------------
*=======================================================================


      CALL prtite
      IF(prttpg) THEN
         WRITE(kprint,1000)
         WRITE(kprint,'(/ / / /)')
         IF(prtseq) CALL prtsq(nbun,mend,prsymb)
         IF(prtatl) CALL prtat(ss_index,ntap,beta,betb,nres,m1,prsymb
     &        ,chrge,mass)
         IF(prtbndl) CALL prtbnd(beta,nres(1,1),lstrtch,lstretch,potbo
     &        ,m9)
         IF(prtcnl) CALL prtcn(beta,nres(1,1),lcnstr,lconstr,potbo_cnst
     &        ,m9)
         IF(prtbal) CALL prtba(beta,nres(1,1),lbndg,lbend,potbe,m2)
         IF(prtptl) CALL prtpt(beta,nres(1,1),ltor,ltors,potto,m3)
         IF(prtitl) CALL prtit(beta,nres(1,1),litr,litor,potit,m4)
         WRITE(kprint,'(/ / / /77(''=''))')
      END IF

*======================================================================
*==== Write force field when required. Only the parameters used -------
*==== in the actual simulation will be printed ------------------------
*======================================================================

      IF(write_ff_pars) THEN
         WRITE(kprint,89000)
         WRITE(kprint,90000)
         CALL write_pot_bond(lstretch,2,lstrtch,m9,2,potbo,betb,1)
         WRITE(kprint,95000)
         CALL write_pot_bond(lconstr,2,lcnstr,m9,2,potbo_cnst,betb,-1)
         WRITE(kprint,94000)
         WRITE(kprint,91000)
         CALL write_pot_bond(lbend,3,lbndg,m2,4,potbe,betb,1)
         WRITE(kprint,94000)
         WRITE(kprint,92000)
         CALL write_pot_bond(ltors,4,ltor,m3,2,potto,betb,1)
         WRITE(kprint,94000)
         WRITE(kprint,93000)
         CALL write_pot_bond(litor,4,litr,m4,3,potit,betb,1)
         WRITE(kprint,94000)
         WRITE(kprint,96000)
         WRITE(kprint,97000)
         CALL write_pot_nbond(ecc6,ecc12,ecc146,ecc1412,nbtype,mass
     &        ,lj_fudgeb,betb,ntap)
         WRITE(kprint,94000)
         WRITE(kprint,89500)
      END IF

*======================================================================
*==== Write a Banner Before intermediate results ----------------------
*======================================================================

      WRITE(kprint,'(/ /)')
      WRITE(kprint,1200)
      WRITE(kprint,1300)
      WRITE(kprint,1110)
      WRITE(kprint,1110)
      WRITE(kprint,1300)
      WRITE(kprint,1200)
      WRITE(kprint,'(/ /)')

      IF(rneih.GT.0.d0) THEN
         lupdate=.true.
      else
         lupdate=.false.
      END IF   

*=======================================================================
*--- Change frame to get xpa, ypa, zpa etc in box fractions ------------
*=======================================================================

      CALL change_frame(co,oc,-1,ntap,xp0,yp0,zp0,xpa,ypa,zpa)
      CALL change_frame(co,oc,-1,ngrp,xpg,ypg,zpg,xpga,ypga,zpga)
      CALL change_frame(co,oc,-1,nprot,xpcm,ypcm,zpcm,xpcma,ypcma,zpcma)

*=======================================================================
*--- Compute lx, ly, lz by subtracting the coordinates of the c.of m.---
*=======================================================================

      CALL change_origin(1,nprot,protl,xp0,yp0,zp0,lx,ly,lz,xpcma,ypcma
     &     ,zpcma,co)

*==== Phony call to forces: Computes only neighbor lists (OLD UPDATE)          

      CALL timer(vfcp,tfcp,elapse)
      gcpu=tfcp
      IF(lupdate) THEN
         IF( .NOT. linked_cell) THEN
*--      update shell h neighbor list
            CALL mts_forces('u',xpa,ypa,zpa,xpga,ypga,zpga,xpcma,ypcma
     &           ,zpcma,mapnl,mapdn,nmapdn,ucns_p,ucos_p,virs_p,virsp_p
     &           ,ucnp_p,ucop_p,ucnsp_p,ucosp_p,fpx_p,fpy_p,fpz_p
     &           ,stressd_p,worka,cpu_h,ncpu_h
     &           ,nstart_h,nend_h,nstart_ah,nend_ah,nlocal_ah,node
     &           ,nprocs,ncube,P_dyn_update_shell)
         ELSE
            aux= rcuth+rtolh+rneih
            CALL lc_index(indxyz,ncx,ncy,ncz,nind,indxi,indxj,indxk,aux
     &           ,co)
            CALL timer(vfcp,tfcp,elapse)
            gcpu=-gcpu + tfcp
            write(kprint,15011) gcpu
            CALL timer(vfcp,tfcp,elapse)
            gcpu=tfcp

            CALL lc_list(ncx,ncy,ncz,nind,indxi,indxj,indxk,aux,co,xpga
     &           ,ypga,zpga,ngrp,nstart_h,nend_h,node,nprocs,ncube
     &           ,worka,kprint,.TRUE.)
         END IF
      END IF

      CALL timer(vfcp,tfcp,elapse)
      gcpu=-gcpu + tfcp
      gcpu_u=gcpu
      write(kprint,16011) gcpu

*=======================================================================
*---- Zeroes all forces ------------------------------------------------
*=======================================================================
      ALLOCATE(fpx_n0(ntap),fpy_n0(ntap),fpz_n0(ntap))
      ALLOCATE(fpx_n1(ntap),fpy_n1(ntap),fpz_n1(ntap))
      ALLOCATE(fpx_m(ntap),fpy_m(ntap),fpz_m(ntap))
      ALLOCATE(fpx_l(ntap),fpy_l(ntap),fpz_l(ntap))
      ALLOCATE(fpx_h(ntap),fpy_h(ntap),fpz_h(ntap))
      ALLOCATE(phi(ntap))

      CALL zeroa(fpx_h,fpy_h,fpz_h,ntap,1)
      CALL zeroa(fpx_l,fpy_l,fpz_l,ntap,1)
      CALL zeroa(fpx_m,fpy_m,fpz_m,ntap,1)
      CALL zeroa(fpx_n1,fpy_n1,fpz_n1,ntap,1)
      CALL zeroa(fpx_n0,fpy_n0,fpz_n0,ntap,1)

      write(kprint,15000) 

*=======================================================================
*===== Computes all forces at time = 0 ================================= 
*=======================================================================

      GTime_h=0.0D0
      GTime_l=0.0D0
      GTime_m=0.0D0

#ifdef PARALLEL

*=======================================================================
*===== Make estimates of communication for forces ====================== 
*=======================================================================

      CALL P_Reduce_Forces(fpx_h,fpy_h,fpz_h)
      CALL P_Reduce_Forces(fpx_h,fpy_h,fpz_h)

      CALL timer(vfcp,tfcp,elapse)
      CTime_h=tfcp
      DO i=1,10
         CALL P_Reduce_Forces(fpx_h,fpy_h,fpz_h)
      END DO
      CALL timer(vfcp,tfcp,elapse)
      IF(pme) THEN
         IF(rshk .NE. 'h') THEN
            GTime_h = (tfcp -CTime_h)/10.0D0
         END IF
         IF(rshk .NE. 'l') THEN
            GTime_l = (tfcp -CTime_h)/10.0D0
         END IF
         IF(rshk .NE. 'm') THEN
            GTime_m = (tfcp -CTime_h)/10.0D0
         END IF
      ELSE
         GTime_l = (tfcp -CTime_h)/10.0D0
         GTime_m = (tfcp -CTime_h)/10.0D0
      END IF

      CALL P_Change_Decomposition("N2-N1",ntap,fpx_h,fpy_h,fpz_h
     &     ,nstart_2,nend_2,nstart_1,nend_1,node
     &     ,nprocs)
      
      CALL P_Change_Decomposition("N2-EX",ntap,fpx_h,fpy_h,fpz_h
     &     ,nstart_2,nend_2,nstart_ex,nend_ex,node
     &     ,nprocs)
      
      CALL P_Change_Decomposition("CM-CME",nprot,fpx_h,fpy_h,fpz_h
     &     ,nstart_cm,nend_cm,nstart_cme
     &     ,nend_cme,node,nprocs)
      CALL P_Change_Decomposition("CM-CMI",nprot,fpx_h,fpy_h,fpz_h
     &     ,nstart_cm,nend_cm,nstart_cmi
     &     ,nend_cmi,node,nprocs)

      CALL timer(vfcp,tfcp,elapse)
      CTime_m=tfcp
      DO i=1,10
         CALL P_Change_Decomposition("N2-N1",ntap,fpx_h,fpy_h,fpz_h
     &        ,nstart_2,nend_2,nstart_1,nend_1,node
     &        ,nprocs)
         CALL P_Change_Decomposition("N2-EX",ntap,fpx_h,fpy_h,fpz_h
     &        ,nstart_2,nend_2,nstart_ex,nend_ex,node
     &        ,nprocs)
         
         CALL P_Change_Decomposition("CM-CME",nprot,fpx_h,fpy_h,fpz_h
     &        ,nstart_cm,nend_cm,nstart_cme
     &        ,nend_cme,node,nprocs)
         CALL P_Change_Decomposition("CM-CMI",nprot,fpx_h,fpy_h,fpz_h
     &        ,nstart_cm,nend_cm,nstart_cmi
     &        ,nend_cmi,node,nprocs)
      END DO
      CALL timer(vfcp,tfcp,elapse)
      GTime_m = GTime_m + 2.0D0*(tfcp -CTime_m)/10.0D0
      
      CALL P_expand_r8x3(fpx_h,fpy_h,fpz_h,nstart_1,nend_1,nlocal_1,node
     &     ,nprocs,1)
      IF(ptr_ex .LT. ntap) THEN
         CALL P_expand_r8x3(fpx_h,fpy_h,fpz_h,nstart_ex0
     &        ,nend_ex0,nlocal_ex0,node,nprocs,ptr_ex)
      END IF

      CALL timer(vfcp,tfcp,elapse)
      CTime_m=tfcp
      DO i=1,10
         CALL P_expand_r8x3(fpx_h,fpy_h,fpz_h,nstart_1,nend_1,nlocal_1
     &        ,node,nprocs,1)
         IF(ptr_ex .LT. ntap) THEN
            CALL P_expand_r8x3(fpx_h,fpy_h,fpz_h,nstart_ex0
     &           ,nend_ex0,nlocal_ex0,node,nprocs,ptr_ex)
         END IF
      END DO
      CALL timer(vfcp,tfcp,elapse)
      GTime_m = GTime_m + (tfcp -CTime_m)/10.0D0

#endif

************************************************************************
***                        Direct SPACE                              ***
***        keep this order to compute nested neighbor lists          ***
************************************************************************

*=======================================================================
*==== H-contribution in direct space
*=======================================================================

      CALL timer(vfcp,tfcp,elapse)
      gcpu=tfcp
      IF(.NOT. polar) THEN

*---  CALLs mts_forces

         rshell='h'
         CALL mts_forces(rshell,xpa,ypa,zpa,xpga,ypga,zpga,xpcma,ypcma
     &        ,zpcma,mapnl,mapdn,nmapdn,ucns_h,ucos_h,virs_h,virsp_h
     &        ,ucnp_h,ucop_h,ucnsp_h,ucosp_h,fpx_h,fpy_h,fpz_h,stressd_h
     &        ,worka,cpu_h,ncpu_h,nstart_h,nend_h
     &        ,nstart_ah,nend_ah,nlocal_ah,node,nprocs,ncube
     &        ,P_dyn_update_shell)
         CALL timer(vfcp,tfcp,elapse)
#if defined PARALLEL
         IF(pme) THEN
            IF(rshk .NE. rshell) THEN
               CALL P_Reduce_Forces(fpx_h,fpy_h,fpz_h)
            END IF
         ELSE
            CALL P_Reduce_Forces(fpx_h,fpy_h,fpz_h)
         END IF
#endif
         gcpu_hd=-gcpu + tfcp
         gcpu_1 = gcpu_hd + GTime_h

*=======================================================================
*==== L-contribution in direct space (long-ranged) ---------------------
*=======================================================================

         CALL timer(vfcp,tfcp,elapse)
         gcpu=tfcp
         
         rshell='l'
         CALL mts_forces(rshell,xpa,ypa,zpa,xpga,ypga,zpga,xpcma,ypcma
     &        ,zpcma,mapnl,mapdn,nmapdn,ucns_l,ucos_l,virs_l,virsp_l
     &        ,ucnp_l,ucop_l,ucnsp_l,ucosp_l,fpx_l,fpy_l,fpz_l,stressd_l
     &        ,worka,cpu_h,ncpu_h,nstart_h,nend_h
     &        ,nstart_ah,nend_ah,nlocal_ah,node,nprocs,ncube
     &        ,P_dyn_update_shell)
         CALL timer(vfcp,tfcp,elapse)
#if defined PARALLEL
         IF(pme) THEN
            IF(rshk .NE. rshell) THEN
               CALL P_Reduce_Forces(fpx_l,fpy_l,fpz_l)
            END IF
         ELSE
            CALL P_Reduce_Forces(fpx_l,fpy_l,fpz_l)
         END IF
#endif
         gcpu_ld=-gcpu + tfcp
         gcpu_2 = gcpu_ld + GTime_l
         
*=======================================================================
*------- M-contribution in direct space (short-ranged) -----------------
*-------          and protein torsion ----------------------------------
*=======================================================================

         CALL timer(vfcp,tfcp,elapse)
         gcpu=tfcp
         rshell='m'
         CALL mts_forces(rshell,xpa,ypa,zpa,xpga,ypga,zpga,xpcma,ypcma
     &        ,zpcma,mapnl,mapdn,nmapdn,ucns_m,ucos_m,virs_m,virsp_m
     &        ,ucnp_m,ucop_m,ucnsp_m,ucosp_m,fpx_m,fpy_m,fpz_m,stressd_m
     &        ,worka,cpu_h,ncpu_h,nstart_h
     &        ,nend_h,nstart_ah,nend_ah,nlocal_ah,node,nprocs,ncube
     &        ,P_dyn_update_shell)
         CALL timer(vfcp,tfcp,elapse)
#if defined PARALLEL
         IF(pme) THEN
            IF(rshk .NE. rshell) THEN
               CALL P_Reduce_Forces(fpx_m,fpy_m,fpz_m)
            ENDIF
         ELSE
            CALL P_Reduce_Forces(fpx_m,fpy_m,fpz_m)
         END IF
#endif

         IF(abmd .AND. (dissociate .OR. associate)) CALL
     &        comp_abmd_fdiss(abmd_dir,rspset,diss_list,spring,xpa,ypa
     &        ,zpa,co,fpx_m,fpy_m,fpz_m,uumb,gr)
         
         CALL timer(vfcp,tfcp,elapse)
         gcpu_md=-gcpu + tfcp
         gcpu_3 = gcpu_md + GTime_m
      END IF
************************************************************************
***                    Reciprocal SPACE                              ***
***     order is reverted to avoid to compute twice pme corrections  ***
************************************************************************

*=======================================================================
*==== M-contribution in reciprocal Space -------------------------------
*=======================================================================

      IF(polar) THEN

         CALL Polarization_Forces(fscnstr_slt,fscnstr_slv,ntap,chrge
     &        ,ss_index,plrzbij,kpol_inp,lbohr,lpnbd,alnbd,betb,nbtype
     &        ,protl,anprot,annpro,anpoint,polar,pol_type,Ext_ef,xp0,yp0
     &        ,zp0,xpa,ypa,zpa,xpga,ypga,zpga,xpcma,ypcma,zpcma,co,oc
     &        ,pmass,mapnl,fpx_m,fpy_m,fpz_m,worka,cpu_h
     &        ,ncpu_h,nstart_h,nend_h,nlocal_h,nstart_ah,nend_ah
     &        ,nlocal_ah,node,nodex,nodey,nodez,ictxt,npy,npz,descQ
     &        ,nprocs,ncube,fstep,tag_bndg,ngrp,grppt,ingrpp,ingrp
     &        ,ingrp_x,errmsg,iret,nstep,unitc,efact,pi,alphal,U_conf
     &        ,U_ele,Udirect,Urecip,Uind,uself_dip,uself,Ugrp,U_Thole
     &        ,ucop_m,ucos_m,ucosp_m,ucnp_m,ucns_m,ucnsp_m,U_solv
     &        ,Old_dipoles,Fixt_Dipoles,mesos,Mesos_Rho,volume
     &        ,Polarization_Model)
         
         eer_m=Urecip
         CALL timer(vfcp,tfcp,elapse)
         gcpu_hd=-gcpu + tfcp
         gcpu_1 = gcpu_hd
         GOTO 908
      END IF
      IF(Poisson_Boltz) THEN
         ALLOCATE(fpx(ntap),fpy(ntap),fpz(ntap))
         CALL Fourier_Init(node,nodex,nodey,nodez,ictxt,npy,npz,descQ
     &        ,nprocs,ncube,nbyte,rbyte,nstart_2,nend_2,nlocal_2,ntap
     &        ,xpa,ypa,zpa,xpcma,ypcma,zpcma,chrge,co,oc,volume,alphal
     &        ,pme_order,nfft1,nfft2,nfft3,nfft3_start,nfft3_local
     &        ,nfft2_start,nfft2_local,eer,fpx,fpy,fpz,phi,stressc,atomp
     &        ,grppt,pressure,rkcut)

         CALL Ions_Boltzmann(Boltzmann)
      END IF
      IF(.NOT.clewld) goto 907 
      gcpu=0.d0
      rshell='m'
      IF(pme .AND. rshell .NE. rshk) goto 914
      CALL timer(vfcp,tfcp,elapse)
      gcpu=tfcp

*---  CALLs mts_furier, mts_furipp, mts_furipw 

      CALL mts_furier(node,nodex,nodey,nodez,ictxt,npy,npz,descQ,nprocs
     &     ,ncube,nstart_1,nend_1,nlocal_1,nstart_2,nend_2,nlocal_2,xp0
     &     ,yp0,zp0,xpa,ypa,zpa,xpcma,ypcma,zpcma,urcsp_m,urcs_m,urcp_m
     &     ,virsp_m,virs_m,virp_m,fpx_m,fpy_m,fpz_m,phi,fsin14,fsbend
     &     ,fsbond,fscnstr_slt,fscnstr_slv,coul_bnd_slt,coul_bnd_slv
     &     ,rshell,rshk,eer_m,stressr_m,fudgec,tag_bndg)

      CALL timer(vfcp,tfcp,elapse)
      gcpu=-gcpu + tfcp
914   write(kprint,16003) gcpu,gcpu_md,gcpu+gcpu_md
      gcpu_1=gcpu+gcpu_md

*=======================================================================
*==== L-contribution in reciprocal Space -------------------------------
*=======================================================================

      gcpu=0.d0
      rshell='l'
      IF(pme .AND. rshell .NE. rshk) goto 915
      CALL timer(vfcp,tfcp,elapse)
      gcpu=tfcp

*--   subtract the less accurate k-forces

      CALL mts_furier(node,nodex,nodey,nodez,ictxt,npy,npz,descQ,nprocs
     &     ,ncube,nstart_1,nend_1,nlocal_1,nstart_2,nend_2,nlocal_2,xp0
     &     ,yp0,zp0,xpa,ypa,zpa,xpcma,ypcma,zpcma,urcsp_l,urcs_l,urcp_l
     &     ,virsp_l,virs_l,virp_l,fpx_l,fpy_l,fpz_l,phi,fsin14,fsbend
     &     ,fsbond,fscnstr_slt,fscnstr_slv,coul_bnd_slt,coul_bnd_slv
     &     ,rshell,rshk,eer_l,stressr_l,fudgec,tag_bndg)

      CALL timer(vfcp,tfcp,elapse)
      gcpu=-gcpu + tfcp
915   write(kprint,16002) gcpu,gcpu_ld,gcpu+gcpu_ld
      gcpu_2=gcpu+gcpu_ld


*=======================================================================
*==== H-contribution in reciprocal Space -------------------------------
*=======================================================================

      gcpu=0.d0
      rshell='h'
      IF(pme .AND. (rshell .NE. rshk)) goto 916
      CALL timer(vfcp,tfcp,elapse)
      gcpu=tfcp

*--   subtract the less accurate k-forces

      CALL mts_furier(node,nodex,nodey,nodez,ictxt,npy,npz,descQ,nprocs
     &     ,ncube,nstart_1,nend_1,nlocal_1,nstart_2,nend_2,nlocal_2,xp0
     &     ,yp0,zp0,xpa,ypa,zpa,xpcma,ypcma,zpcma,urcsp_h,urcs_h,urcp_h
     &     ,virsp_h,virs_h,virp_h,fpx_h,fpy_h,fpz_h,phi,fsin14,fsbend
     &     ,fsbond,fscnstr_slt,fscnstr_slv,coul_bnd_slt,coul_bnd_slv
     &     ,rshell,rshk,eer_h,stressr_h,fudgec,tag_bndg)

      CALL timer(vfcp,tfcp,elapse)

      gcpu=-gcpu + tfcp
916   write(kprint,16001) gcpu,gcpu_hd,gcpu+gcpu_hd
      gcpu_3=gcpu+gcpu_hd

#if defined PARALLEL
      WRITE(kprint,16005) GTime_h,GTime_l,GTime_m
#endif
907   theoric_speed_up=(gcpu_1+gcpu_2+gcpu_3)*mrespa*lrespa/(gcpu_1
     &     *mrespa*lrespa+gcpu_2*lrespa+gcpu_3)

908   CONTINUE
c$$$      CALL CompElecPotentialOnGrid(co,xp0,yp0,zp0)

*=======================================================================
*--- Compute force on the co matrix ------------------------------------
*=======================================================================

      IF(cpress) THEN
         CALL comp_stress_conf(stressd_m,stressr_m,prt_m,oc,volume,unitp
     &        ,press_m)
         CALL comp_stress_kinetic(vcax,vcay,vcaz,tmass,co,1,nprot,volume
     &        ,unitp,st_m,press_kin)
         CALL comp_forcep(prt_m,st_m,oc,volume,pext)
         CALL comp_stress_conf(stressd_l,stressr_l,prt_l,oc,volume,unitp
     &        ,press_l)
         CALL comp_stress_conf(stressd_h,stressr_h,prt_h,oc,volume,unitp
     &        ,press_h)
      END IF

*=======================================================================
*--- Print timing for spherical cutoff ---------------------------------
*=======================================================================

      IF(.NOT.clewld) THEN
         gcpu=0.d0
         write(kprint,16001) gcpu,gcpu_hd,gcpu_hd

         write(kprint,16002) gcpu,gcpu_ld,gcpu_ld
         write(kprint,16003) gcpu,gcpu_md,gcpu_md
         gcpu_1=gcpu+gcpu_md
         gcpu_2=gcpu+gcpu_ld
         gcpu_3=gcpu+gcpu_hd
      END IF   
      write(kprint,10067) theoric_speed_up

************************************************************************
***            FAST INTRAMOLECULAR COMPONENTS                        ***
************************************************************************

*=======================================================================
*---- Computes bonded forces of shell n1 -------------------------------
*=======================================================================

      CALL timer(vfcp,tfcp,elapse)
      gcpu=tfcp

      CALL mts_intra_n1(xp0,yp0,zp0,xpcma,ypcma,zpcma,fpx_n1,fpy_n1
     &     ,fpz_n1,fudge,lj_fudge,abmd_dir,puhyd,conf_bnd_slt_n1
     &     ,conf_bnd_slv_n1,coul_bnd_slt_n1,coul_bnd_slv_n1,unb14,cnb14
     &     ,ungrp,cngrp,uptors,uslvtor,stressd_n1,mapdn,nmapdn,uumb,gr
     &     ,nstart_1,nend_1,node,nprocs,ncube)
      IF(abmd) THEN
         IF(abmd_native) CALL comp_abmd_fnative(abmd_unbias,abmd_native
     &        ,abmd_dir,nato_slt,xp0,yp0,zp0,vpx,vpy,vpz,co,atomp,xpcma
     &        ,ypcma,zpcma,uumb,rspset,spring,native_theta,atres_map1
     &        ,atres_map2,nat_listp,nat_list,native_dist,gr,gra,dfree
     &        ,fpx_n1,fpy_n1,fpz_n1,stressd_n1,nstart_1,nend_1,nlocal_1
     &        ,node,nprocs,ncube)
         IF(abmd_cryst) CALL comp_abmd_fcryst(abmd_cryst_dir,rspset
     &        ,abmd_cryst_nvect,abmd_cryst_vect,spring,ntap,xpa,ypa,zpa
     &        ,co,fpx_n1,fpy_n1,fpz_n1,uumb,gr)
         dentr=0.0D0
         DO i=1,neta
            IF(.NOT. near0(qmass(i))) dentr=dentr
     &           -DBLE(ndf_thermos(i))*gascon*t*etap(i)
         END DO
         tn1=tm/DBLE(n1respa)
         a_free=a_free+dfree*tn1
         ts_free=ts_free+dentr*tn1
      END IF

      CALL timer(vfcp,tfcp,elapse)
      gcpu=-gcpu + tfcp
      write(kprint,14004) gcpu
      gcpu_0=gcpu

*=======================================================================
*---- Computes bonded forces of shell n0 -------------------------------
*=======================================================================

      CALL timer(vfcp,tfcp,elapse)
      gcpu=tfcp

      CALL mts_intra_n0(xp0,yp0,zp0,xpcma,ypcma,zpcma,fpx_n0,fpy_n0
     &     ,fpz_n0,ubond,uslvbon,ubend,uslvben,uitors,uslvitor
     &     ,stressd_n0,tag_bndg,nstart_1,nend_1,nlocal_1,ntot_1,node
     &     ,nprocs,ncube)

      CALL timer(vfcp,tfcp,elapse)

*=======================================================================
*---- Compute stress tensors from bonded interaction -------------------
*=======================================================================

      IF(cpress) THEN
         CALL comp_stress_conf(stressd_n1,stressr_n1,prt_n1,oc,volume
     &        ,unitp,press_n1)
         CALL comp_stress_conf(stressd_n0,stressr_n0,prt_n0,oc,volume
     &        ,unitp,press_n0)
      END IF

*=======================================================================
*--- Compute force on the c.of.m of each molecule and reduce the -------
*--- atomic force. Do so for all shell forces --------------------------
*=======================================================================

      ALLOCATE(fcax_n0(nprot),fcay_n0(nprot),fcaz_n0(nprot))
      ALLOCATE(fcax_n1(nprot),fcay_n1(nprot),fcaz_n1(nprot))
      ALLOCATE(fcax_m(nprot),fcay_m(nprot),fcaz_m(nprot))
      ALLOCATE(fcax_l(nprot),fcay_l(nprot),fcaz_l(nprot))
      ALLOCATE(fcax_h(nprot),fcay_h(nprot),fcaz_h(nprot))

      CALL comp_fcm(nstart_cm,nend_cm,protl,fpx_m,fpy_m,fpz_m
     &     ,fcax_m,fcay_m,fcaz_m,mass,tmass,oc)

      CALL comp_fcm(nstart_cm,nend_cm,protl,fpx_l,fpy_l,fpz_l,fcax_l
     &     ,fcay_l,fcaz_l,mass,tmass,oc)

      CALL comp_fcm(nstart_cm,nend_cm,protl,fpx_h,fpy_h,fpz_h,fcax_h
     &     ,fcay_h,fcaz_h,mass,tmass,oc)

      IF(coupl_grp) THEN
         CALL comp_fcm(nstart_cmi,nend_cmi,protl,fpx_n1,fpy_n1,fpz_n1
     &        ,fcax_n1,fcay_n1,fcaz_n1,mass,tmass,oc)
         CALL comp_fcm(nstart_cmi,nend_cmi,protl,fpx_n0,fpy_n0,fpz_n0
     &        ,fcax_n0,fcay_n0,fcaz_n0,mass,tmass,oc)
      ELSE
         CALL zeroa(fcax_n1,fcay_n1,fcaz_n1,nprot,1)
         CALL zeroa(fcax_n0,fcay_n0,fcaz_n0,nprot,1)
      END IF

*=======================================================================
*---- Print timing and simulation CPU time estimates -------------------
*=======================================================================
      
      gcpu=-gcpu + tfcp
      write(kprint,16004) gcpu
      theoric_speed_up=(gcpu+gcpu_0+gcpu_1+gcpu_2+gcpu_3)*n0respa
     &     *n1respa*mrespa*lrespa/(gcpu*n0respa*n1respa*mrespa*lrespa
     &     +gcpu_0*n1respa*mrespa*lrespa+gcpu_1*mrespa*lrespa+gcpu_2
     &     *lrespa+gcpu_3)
      write(kprint,10068) theoric_speed_up
      if(nupdte.gt.0) gcpu_u = gcpu_u /DBLE(nupdte)
      aux = (gcpu*n0respa*n1respa*lrespa*mrespa+gcpu_0*n1respa*lrespa
     &     *mrespa+gcpu_1*mrespa*lrespa+gcpu_2*lrespa + gcpu_3 + gcpu_u)
      gcpu_0 = aux/DBLE(lrespa*mrespa)
      gcpu_1 = aux / time
      aux = aux*(maxstp+(mrject/DBLE(lrespa*mrespa)))
      hours = int(aux/3600.d0)
      min = int((aux-hours*3600)/60.d0)
      write(kprint,80033) hours,min,gcpu_0,gcpu_1
      write(kprint,60030) 

*=======================================================================
*===== MD loop starts here =============================================
*=======================================================================

      CALL timer(vfcp,tfcp,elapse)
      gcpu=tfcp
      elaps=elapse
      lfirst=.true.
      counter=0
100   CONTINUE
      
      CALL timer(vfcp_hh,tfcp_hh,elapse)
      tdelta_hh=tfcp_hh

      nstep=nstep+1
      mstep=mstep+1
      n_q=ntime_q
      IF(MOD(nstep,ntime_q) .NE. 0) n_q=MOD(nstep,ntime_q)
      time=time_q(n_q)
      time2=time*0.5D0

*=======================================================================
*---- Start 4 time steps r-RESPA integrator ----------------------------
*=======================================================================

*=======================================================================
*==== Phony call to forces: Computes only neighbor lists (OLD UPDATE) --
*=======================================================================

      IF(MOD(nstep,nupdte) .EQ. 0 .AND. lupdate) THEN
*=======================================================================
*--  1)  Update shell *H* neighbor list
*=======================================================================
         
#if defined PARALLEL
*----- Now can split and set up correctly
         CALL P_load_balance(P_dyn_update,node,nprocs,ncube,nbyte,rbyte
     &        ,nstart_h,nend_h,nlocal_h,nstart_ah,nend_ah,nlocal_ah,ntap
     &        ,ngrp,grppt,worka,cpu_h,ncpu_h,nmapnl,mapnl,mapnl_save)
#endif

         IF(.not.linked_cell) THEN
            CALL mts_forces('u',xpa,ypa,zpa,xpga,ypga,zpga,xpcma,ypcma
     &           ,zpcma,mapnl,mapdn,nmapdn,ucns_p,ucos_p,virs_p,virsp_p
     &           ,ucnp_p,ucop_p,ucnsp_p,ucosp_p,fpx_p,fpy_p,fpz_p
     &           ,stressd_p,worka,cpu_h,ncpu_h
     &           ,nstart_h,nend_h,nstart_ah,nend_ah,nlocal_ah,node
     &           ,nprocs,ncube,P_dyn_update_shell)
         ELSE   
            aux=rcuth+rtolh+rneih 
            IF(cpress .AND. MOD(nstep,nupdte*nupdte_index) .EQ. 0) THEN
               CALL lc_index(indxyz,ncx,ncy,ncz,nind,indxi,indxj,indxk
     &              ,aux,co)
               WRITE(kprint,2000)
            END IF
            CALL lc_list(ncx,ncy,ncz,nind,indxi,indxj,indxk,aux,co,xpga
     &           ,ypga,zpga,ngrp,nstart_h,nend_h,node,nprocs,ncube
     &           ,worka,kprint,.TRUE.)
         END IF 

         rcuth_save=rcuth
         rtolh_save=rtolh
         rcuth=rcutl
         rtolh=0.0D0


*=======================================================================
*--  2)  Update shell *L* neighbor list
*=======================================================================
         
         CALL mts_forces('h',xpa,ypa,zpa,xpga,ypga,zpga,xpcma,ypcma
     &        ,zpcma,mapnl,mapdn,nmapdn,ucns_p,ucos_p,virs_p,virsp_p
     &        ,ucnp_p,ucop_p,ucnsp_p,ucosp_p,fpx_p,fpy_p,fpz_p,stressd_p
     &        ,worka,cpu_h,ncpu_h,nstart_h
     &        ,nend_h,nstart_ah,nend_ah,nlocal_ah,node,nprocs,ncube
     &        ,P_dyn_update_shell)

         rcutl_save=rcutl
         rtoll_save=rtoll
         rcutl=rcutm
         rtoll=0.0D0
         
*=======================================================================
*--  3)  Update shell *M* neighbor list
*=======================================================================
         
         CALL mts_forces('l',xpa,ypa,zpa,xpga,ypga,zpga,xpcma,ypcma
     &        ,zpcma,mapnl,mapdn,nmapdn,ucns_p,ucos_p,virs_p,virsp_p
     &        ,ucnp_p,ucop_p,ucnsp_p,ucosp_p,fpx_p,fpy_p,fpz_p,stressd_p
     &        ,worka,cpu_h,ncpu_h,nstart_h
     &        ,nend_h,nstart_ah,nend_ah,nlocal_ah,node,nprocs,ncube
     &        ,P_dyn_update_shell)

         rcuth=rcuth_save
         rtolh=rtolh_save
         rcutl=rcutl_save
         rtoll=rtoll_save

      END IF   
      
*=======================================================================
*---  Advances velocities for half time step TIME using H-forces    ----
*=======================================================================
      
      CALL correc(vpx,vpy,vpz,fpx_h,fpy_h,fpz_h,mass,nstart_2,nend_2
     &     ,time)
      
      IF(cnstpp .NE. 0) THEN
         CALL rattle_correc(nstart_2,nend_2,time,xp0,yp0,zp0,vpx
     &        ,vpy,vpz,ntap,cnstp,dssp,coeffp,cnstpp,mass,dnit
     &        ,cnst_protp,cnst_protl,mim_lim,gcpu_rt,iret
     &        ,errmsg)
         IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
      END IF

      CALL correc_stress(cpress,h_skin,nprot,co,oc,vco,vcax,vcay
     &     ,vcaz,fcax_h,fcay_h,fcaz_h,stressd_h,stressr_h,volume
     &     ,press_h,press_kin,pext,tmass,masspp,time,time2,nstart_cm
     &     ,nend_cm,node,nprocs,ncube,rbyte)

      IF(isostress  .OR. FixedAngles_Stress) THEN
         CALL rattle_correc_co(co,dssco,cnstco,vco,masspp,iret,errmsg)
         IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
      END IF
      
*=======================================================================
*---  Reset h-forces to ZERO                                        ----
*=======================================================================
      
      CALL zeroa(fpx_h,fpy_h,fpz_h,ntap,1)
      tl = time/DBLE(lrespa)
      tl2=tl*0.5D0
      CALL timer(vfcp_hh,tfcp_hh,elapse)
      tdelta_hh=tfcp_hh-tdelta_hh
      gcpu_hh=gcpu_hh + tdelta_hh
      
*=======================================================================
*---     Advances velocities for half time step *time* using        ----
*--- |-->         long range L-forces                               ----
*=======================================================================
      
      DO il=1,lrespa
         
         CALL timer(vfcp_ll,tfcp_ll,elapse)
         tdelta_ll=tfcp_ll

         CALL correc(vpx,vpy,vpz,fpx_l,fpy_l,fpz_l,mass,nstart_2,nend_2
     &        ,tl)
         IF(cnstpp .NE. 0) THEN
            CALL rattle_correc(nstart_2,nend_2,tl,xp0,yp0,zp0,vpx
     &           ,vpy,vpz,ntap,cnstp,dssp,coeffp,cnstpp,mass,dnit
     &           ,cnst_protp,cnst_protl,mim_lim,gcpu_rt,iret
     &           ,errmsg)
            IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
         END IF

C Need to change in parallel runs 

         CALL correc_stress(cpress,l_skin,nprot,co,oc,vco,vcax,vcay
     &        ,vcaz,fcax_l,fcay_l,fcaz_l,stressd_l,stressr_l,volume
     &        ,press_l,press_kin,pext,tmass,masspp,tl,tl2,nstart_cm
     &        ,nend_cm,node,nprocs,ncube,rbyte)

         IF(isostress   .OR. FixedAngles_Stress) THEN
            CALL rattle_correc_co(co,dssco,cnstco,vco,masspp
     &           ,iret,errmsg)
            IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
         END IF

*=======================================================================
*---     Reset l-forces to ZERO                                    -----
*=======================================================================
         
         CALL zeroa(fpx_l,fpy_l,fpz_l,ntap,1)
         tm = tl/DBLE(mrespa)
         tm2=tm*0.5D0
         tm4=tm*0.25D0
         
         rshell='m'
         
         CALL timer(vfcp_ll,tfcp_ll,elapse)
         tdelta_ll=tfcp_ll-tdelta_ll
         gcpu_ll=gcpu_ll + tdelta_ll

*=======================================================================
*---        Advances velocities for half time step *time* using     ----
*---    |-->           short range M-forces                         ----
*=======================================================================
         
         DO im=1,mrespa
            CALL timer(vfcp_mm,tfcp_mm,elapse)
            tdelta_mm=tfcp_mm

            ninner=ninner + 1
            nsstt = nsstt + 1
            
            IF(thermos) THEN
c----------  Correct velocities of the NETA thermostats ----------------
c----------  Propagate exp(i L_u tm/4) ---------------------------------
               IF(isostress   .OR. FixedAngles_Stress) THEN
                  CALL rattle_correc_co(co,dssco,cnstco,vco,masspp,iret
     &                 ,errmsg)
                  IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
               END IF
               IF(cnstpp .NE. 0) THEN
                  CALL rattle_correc(nstart_2,nend_2,tm,xp0,yp0,zp0
     &                 ,vpx,vpy,vpz,ntap,cnstp,dssp,coeffp,cnstpp,mass
     &                 ,dnit,cnst_protp,cnst_protl,mim_lim,gcpu_rt
     &                 ,iret,errmsg)
                  IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
               END IF
*--------   Compute forces on the thermostat
               CALL comp_thermos_forces(nprocs,nstart_2,nend_2,nstart_cm
     &              ,nend_cm,ndf_thermos,ss_index,co,vpx,vpy,vpz,vcax
     &              ,vcay,vcaz,mass,tmass,t,fth)
#if defined PARALLEL
               CALL P_merge_vecr8(fth,3)
#endif
*--------   Add contribution from barostat temperature
               IF(cpress) CALL comp_thermos_forces_cm(fth(1),vco,masspr)

               CALL correc_etap(neta,etap,fth,qmass,tm4)
            END IF
      
*----------  Propagate exp(i L_y tm/2) ---------------------------------
            IF(cpress .AND. coupl_mol) THEN
               CALL correc_matr(tm2,co,oc,vco,gmgp,vcax,vcay,vcaz,nprot)
            END IF
            IF(thermos) THEN

c----------  Correct velocities of the NETA thermostats ----------------
c----------  Recompute forces on the thermostat variables --------------
               CALL comp_thermos_forces(nprocs,nstart_2,nend_2,nstart_cm
     &              ,nend_cm,ndf_thermos,ss_index,co,vpx,vpy,vpz,vcax
     &              ,vcay,vcaz,mass,tmass,t,fth)
#if defined PARALLEL
               CALL P_merge_vecr8(fth,3)
#endif
*--------   Add contribution from barostat temperature
               IF(cpress) CALL comp_thermos_forces_cm(fth(1),vco,masspr)

*==========  Propagate exp(i L_u tm/4) =================================
               CALL correc_etap(neta,etap,fth,qmass,tm4)
            END IF
            
            IF(thermos) THEN
*=========  Propagate exp(i L_z tm/4) ==================================
               CALL correc_exp_scale(nstart_2,nend_2,nstart_cm,nend_cm
     &              ,cpress,ss_index,etap,tm4,vcax,vcay,vcaz,vpx,vpy,vpz
     &              ,vco)
            END IF

*========= Propagate exp(i L_x tm/2) -==================================

            CALL correc(vpx,vpy,vpz,fpx_m,fpy_m,fpz_m,mass,nstart_2
     &           ,nend_2,tm)
            IF(cnstpp .NE. 0) THEN
               CALL rattle_correc(nstart_2,nend_2,tm,xp0,yp0,zp0
     &              ,vpx,vpy,vpz,ntap,cnstp,dssp,coeffp,cnstpp,mass,dnit
     &              ,cnst_protp,cnst_protl,mim_lim,gcpu_rt,iret
     &              ,errmsg)
               IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
            END IF
            CALL correc_stress(cpress,m_skin,nprot,co,oc,vco,vcax,vcay
     &           ,vcaz,fcax_m,fcay_m,fcaz_m,stressd_m,stressr_m,volume
     &           ,press_m,press_kin,pext,tmass,masspp,tm,tm2,nstart_cm
     &           ,nend_cm,node,nprocs,ncube,rbyte)
            IF(isostress  .OR. FixedAngles_Stress) THEN
               CALL rattle_correc_co(co,dssco,cnstco,vco,masspp,iret
     &              ,errmsg)
               IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
            END IF
            
*=========  Propagate exp(i L_z tm/4) ==================================
            IF(thermos) THEN
               CALL correc_exp_scale(nstart_2,nend_2,nstart_cm,nend_cm
     &              ,cpress,ss_index,etap,tm4,vcax,vcay,vcaz,vpx,vpy,vpz
     &              ,vco)
            END IF
            
*=======================================================================
*---        Reset m-forces to ZERO
*=======================================================================
            
            CALL zeroa(fpx_m,fpy_m,fpz_m,ntap,1)
            tn1 = tm/DBLE(n1respa)
            tn12=tn1*0.5D0

#if defined PARALLEL           
            nmol_cm=atomp(ntap)

            CALL Get_vp_with_vc(vpx,vpy,vpz,co,nstart_cm,nend_cm,nprot
     &           ,protl,vcax,vcay,vcaz)

            CALL P_Change_Decomposition("N2-N1",ntap,vpx,vpy,vpz
     &           ,nstart_2,nend_2,nstart_1,nend_1,node
     &           ,nprocs)
            CALL P_Change_Decomposition("N2-EX",ntap,vpx,vpy,vpz
     &           ,nstart_2,nend_2,nstart_ex,nend_ex,node
     &           ,nprocs)

            CALL Get_vp_without_vc(vpx,vpy,vpz,oc,nstart_cmi,nend_cmi
     &           ,nprot,protl,mass,tmass,vcax,vcay,vcaz)
            CALL Get_vp_without_vc(vpx,vpy,vpz,oc,nstart_cme,nend_cme
     &           ,nprot,protl,mass,tmass,vcax,vcay,vcaz)

#endif

*=======================================================================
*---                Advance position and velocities for a full      ----
*---       |-->     time step tn1 using intramolecular forces       ----
*=======================================================================

            CALL timer(vfcp_mm,tfcp_mm,elapse)
            tdelta_mm=tfcp_mm-tdelta_mm
            gcpu_mm=gcpu_mm + tdelta_mm

            CALL timer(vfcp_inn,tfcp_inn,elapse)
            tdelta_inn=tfcp_inn

            DO in1=1,n1respa
               ninn1=ninn1+1
               CALL correc(vpx,vpy,vpz,fpx_n1,fpy_n1,fpz_n1,mass
     &              ,nstart_1,nend_1,tn1)
               IF(cnstpp .NE. 0) THEN
                  CALL rattle_correc(nstart_1,nend_1,tn1,lx,ly,lz
     &                 ,vpx,vpy,vpz,ntap,cnstp,dssp,coeffp,cnstpp,mass
     &                 ,dnit,cnst_protp_1,cnst_protl_1,mim_lim,gcpu_rt
     &                 ,iret,errmsg)
                  IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
               END IF
               CALL correc_stress(cpress,n1_skin,nprot,co,oc,vco,vcax
     &              ,vcay,vcaz,fcax_n1,fcay_n1,fcaz_n1,stressd_n1
     &              ,stressr_n1,volume,press_n1,press_kin,pext,tmass
     &              ,masspp,tn1,tn12,nstart_cmi,nend_cmi,node,nprocs
     &              ,ncube,rbyte)

               IF(isostress  .OR. FixedAngles_Stress) THEN
                  CALL rattle_correc_co(co,dssco,cnstco,vco,masspp,iret
     &                 ,errmsg)
                  IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
               END IF

               CALL zeroa(fpx_n1(nstart_1),fpy_n1(nstart_1)
     &              ,fpz_n1(nstart_1),nlocal_1,1)
               tn0 = tn1/DBLE(n0respa)
               tn02=tn0*0.5D0

               DO in0=1,n0respa
                  ninn0=ninn0+1
                  
                  CALL dcopy(nlocal_1,lx(nstart_1),1,xpo(nstart_1),1)
                  CALL dcopy(nlocal_1,ly(nstart_1),1,ypo(nstart_1),1)
                  CALL dcopy(nlocal_1,lz(nstart_1),1,zpo(nstart_1),1)
                  CALL dcopy(nlocal_ex,lx(nstart_ex),1,xpo(nstart_ex),1)
                  CALL dcopy(nlocal_ex,ly(nstart_ex),1,ypo(nstart_ex),1)
                  CALL dcopy(nlocal_ex,lz(nstart_ex),1,zpo(nstart_ex),1)

                  if(.not.start_conf)  then
                     CALL verlet(mass,nstart_1,nend_1,lx,ly,lz,vpx,vpy
     &                    ,vpz,fpx_n0,fpy_n0,fpz_n0,tn0)
                     CALL verlet(mass,nstart_ex,nend_ex,lx,ly,lz,vpx,vpy
     &                    ,vpz,fpx_n0,fpy_n0,fpz_n0,tn0)
                  else
                     CALL starting_verlet(mass,nstart_1,nend_1,lx,ly,lz
     &                    ,vpx,vpy,vpz,fpx_n0,fpy_n0,fpz_n0,tn0,max_dist
     &                    )
                     CALL starting_verlet(mass,nstart_ex,nend_ex,lx,ly
     &                    ,lz,vpx,vpy,vpz,fpx_n0,fpy_n0,fpz_n0,tn0
     &                    ,max_dist)
                  endif
                  IF(cnstpp .NE. 0) THEN
                     CALL rattle_verlet(nstart_1,nend_1,tn0,lx,ly
     &                    ,lz,xpo,ypo,zpo,vpx,vpy,vpz,ntap,cnstp,dssp
     &                    ,coeffp,cnstpp,mass,dnit,cnst_protp_1
     &                    ,cnst_protl_1,mim_lim,gcpu_sh,iret,errmsg)
                     IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
                     CALL rattle_verlet(nstart_ex,nend_ex,tn0,lx,ly
     &                    ,lz,xpo,ypo,zpo,vpx,vpy,vpz,ntap,cnstp,dssp
     &                    ,coeffp,cnstpp,mass,dnit,cnst_protp_ex
     &                    ,cnst_protl_ex,mim_lim,gcpu_sh,iret,errmsg)
                     IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
                  END IF

                  CALL correc_stress_n0(cpress,n0_skin,co,oc,vco,vcax
     &                 ,vcay,vcaz,fcax_n0,fcay_n0,fcaz_n0,stressd_n0
     &                 ,stressr_n0,volume,press_n0,press_kin,pext,tmass
     &                 ,masspp,tn0,tn02,nstart_cmi,nend_cmi,nstart_cme
     &                 ,nend_cme,node,nprocs,ncube,rbyte)

                  IF(isostress  .OR. FixedAngles_Stress) THEN
                     CALL rattle_correc_co(co,dssco,cnstco,vco,masspp
     &                    ,iret,errmsg)
                     IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
                  END IF

*----------  Propagate exp(i L_y tn1/2) ---------------------------------
                  IF(cpress .AND. coupl_grp) THEN
                     CALL correc_matr(tn02,co,oc,vco,gmgp
     &                    ,vcax(nstart_cmi),vcay(nstart_cmi)
     &                    ,vcaz(nstart_cmi),nlocal_cmi)
                     CALL correc_matr(tn02,co,oc,vco,gmgp
     &                    ,vcax(nstart_cme),vcay(nstart_cme)
     &                    ,vcaz(nstart_cme),nlocal_cme)
                  END IF
 
                  CALL verlet_free(nstart_cmi,nend_cmi,xpcma,ypcma,zpcma
     &                 ,vcax,vcay,vcaz,tn0)
                  CALL verlet_free(nstart_cme,nend_cme,xpcma,ypcma,zpcma
     &                 ,vcax,vcay,vcaz,tn0)

*---              Reset n-forces to ZERO
                  
                  CALL zeroa(fpx_n0(nstart_1),fpy_n0(nstart_1)
     &                 ,fpz_n0(nstart_1),nlocal_1,1)
                  
*=======================================================================
*---           Compute N0-forces                                    ----
*=======================================================================
                  
                  IF(cpress) THEN
                     DO i=1,3
                        DO j=1,3
                           coo(i,j)=co(i,j)
                        END DO
                     END DO
                     
*--- Propagate co matrix for step tm -----------------------------------
                     
                     CALL verlet_free(1,3,co(1,1),co(1,2),co(1,3),vco(1
     &                    ,1),vco(1,2),vco(1,3),tn0)
                     IF(isostress  .OR. FixedAngles_Stress) THEN
                        CALL rattle_verlet_co(tm,co,coo,dssco,cnstco,vco
     &                       ,masspp,iret,errmsg)
                        IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
                     END IF
                     
                     CALL matinv(3,3,co,oc,volume)
                     volume=volume*boxl**3
                  END IF

*=======================================================================
*---  Recompute coordinates before recomputing forces ------------------
*=======================================================================

#if defined PARALLEL            

*=======================================================================
*---- Expand coordinates lx ly lz to compute force on n0 and n1 shell --
*=======================================================================


                  CALL P_Comm_Intra(nstart_1,nend_1,node,nprocs,lx,ly,lz
     &                 ,xpcma,ypcma,zpcma)

c --- 
#endif
                  CALL change_coord_inner(ntot_1,ntot_cmi,protl,co,oc
     &                 ,pmass,lx,ly,lz,xp0,yp0,zp0,xpa,ypa,zpa,xpcma
     &                 ,ypcma,zpcma)
                  
                  CALL mts_intra_n0(xp0,yp0,zp0,xpcma,ypcma,zpcma,fpx_n0
     &                 ,fpy_n0,fpz_n0,ubond,uslvbon,ubend,uslvben,uitors
     &                 ,uslvitor,stressd_n0,tag_bndg,nstart_1,nend_1
     &                 ,nlocal_1,ntot_1,node,nprocs,ncube)

                  CALL comp_fcm(nstart_cmi,nend_cmi,protl,fpx_n0,fpy_n0
     &                 ,fpz_n0,fcax_n0,fcay_n0,fcaz_n0,mass,tmass,oc)
                  
                  CALL correc(vpx,vpy,vpz,fpx_n0,fpy_n0,fpz_n0,mass
     &                 ,nstart_1,nend_1,tn0)
                  
                  IF(cnstpp .NE. 0) THEN
                     CALL rattle_correc(nstart_1,nend_1,tn0,lx,ly
     &                    ,lz,vpx,vpy,vpz,ntap,cnstp,dssp,coeffp,cnstpp
     &                    ,mass,dnit,cnst_protp_1,cnst_protl_1,mim_lim
     &                    ,gcpu_rt,iret,errmsg)
                     CALL rattle_correc(nstart_ex,nend_ex,tn0,lx,ly
     &                    ,lz,vpx,vpy,vpz,ntap,cnstp,dssp,coeffp,cnstpp
     &                    ,mass,dnit,cnst_protp_ex,cnst_protl_ex,mim_lim
     &                    ,gcpu_rt,iret,errmsg)
                     IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
                  END IF
                  
*----------  Propagate exp(i L_y tn1/2) ---------------------------------

                  IF(cpress .AND. coupl_grp) THEN
                     CALL correc_matr(tn02,co,oc,vco,gmgp
     &                    ,vcax(nstart_cmi),vcay(nstart_cmi)
     &                    ,vcaz(nstart_cmi),nlocal_cmi)
                     CALL correc_matr(tn02,co,oc,vco,gmgp
     &                    ,vcax(nstart_cme),vcay(nstart_cme)
     &                    ,vcaz(nstart_cme),nlocal_cme)
                  END IF
                  CALL correc_stress_n0(cpress,n0_skin,co,oc,vco,vcax
     &                 ,vcay,vcaz,fcax_n0,fcay_n0,fcaz_n0,stressd_n0
     &                 ,stressr_n0,volume,press_n0,press_kin,pext,tmass
     &                 ,masspp,tn0,tn02,nstart_cmi,nend_cmi,nstart_cme
     &                 ,nend_cme,node,nprocs,ncube,rbyte)
               END DO
               CALL mts_intra_n1(xp0,yp0,zp0,xpcma,ypcma,zpcma,fpx_n1
     &              ,fpy_n1,fpz_n1,fudge,lj_fudge,abmd_dir,puhyd
     &              ,conf_bnd_slt_n1,conf_bnd_slv_n1,coul_bnd_slt_n1
     &              ,coul_bnd_slv_n1,unb14,cnb14,ungrp,cngrp,uptors
     &              ,uslvtor,stressd_n1,mapdn,nmapdn,uumb,gr,nstart_1
     &              ,nend_1,node,nprocs,ncube)
               IF(abmd) THEN
                  IF(abmd_native) CALL comp_abmd_fnative(abmd_unbias
     &                 ,abmd_native,abmd_dir,nato_slt,xp0,yp0,zp0,vpx
     &                 ,vpy,vpz,co,atomp,xpcma,ypcma,zpcma,uumb,rspset
     &                 ,spring,native_theta,atres_map1,atres_map2
     &                 ,nat_listp,nat_list,native_dist,gr,gra,dfree
     &                 ,fpx_n1,fpy_n1,fpz_n1,stressd_n1,nstart_1,nend_1
     &                 ,nlocal_1,node,nprocs,ncube)
                  IF(abmd_cryst) CALL
     &              comp_abmd_fcryst(abmd_cryst_dir,rspset
     &              ,abmd_cryst_nvect,abmd_cryst_vect,spring,ntap,xpa
     &              ,ypa,zpa,co,fpx_n1,fpy_n1,fpz_n1,uumb,gr)
                  dentr=0.0D0
                  DO i=1,neta
                     IF(.NOT. near0(qmass(i))) dentr=dentr
     &                    -DBLE(ndf_thermos(i))*gascon*t*etap(i)
                  END DO
                  a_free=a_free+dfree*tn1
                  ts_free=ts_free+dentr*tn1
               END IF

               CALL comp_fcm(nstart_cmi,nend_cmi,protl,fpx_n1,fpy_n1
     &              ,fpz_n1,fcax_n1,fcay_n1,fcaz_n1,mass,tmass,oc)
               
               CALL correc(vpx,vpy,vpz,fpx_n1,fpy_n1,fpz_n1,mass
     &              ,nstart_1,nend_1,tn1)

               CALL correc_stress(cpress,n1_skin,nprot,co,oc,vco,vcax
     &              ,vcay,vcaz,fcax_n1,fcay_n1,fcaz_n1,stressd_n1
     &              ,stressr_n1,volume,press_n1,press_kin,pext,tmass
     &              ,masspp,tn1,tn12,nstart_cmi,nend_cmi,node,nprocs
     &              ,ncube,rbyte)
            END DO
            
            CALL timer(vfcp_inn,tfcp_inn,elapse)
            tdelta_inn=tfcp_inn-tdelta_inn
            gcpu_inn=gcpu_inn+tdelta_inn

*=======================================================================
*          |-->  END of n-loop                                       ---
*=======================================================================
            
            CALL timer(vfcp_mm,tfcp_mm,elapse)
            tdelta_mm=tfcp_mm

#if defined PARALLEL
*=======================================================================
*---- Expand vpx,vpy,vpz and vcax,vcay,vcaz of the inner loops ---------
*=======================================================================

            CALL Get_vp_with_vc(vpx,vpy,vpz,co,nstart_cme,nend_cme,nprot
     &           ,protl,vcax,vcay,vcaz)
            CALL Get_vp_with_vc(vpx,vpy,vpz,co,nstart_cmi,nend_cmi,nprot
     &           ,protl,vcax,vcay,vcaz)

            CALL P_Change_Decomposition("EX-N2",ntap,vpx,vpy,vpz
     &           ,nstart_ex,nend_ex,nstart_2,nend_2,node
     &           ,nprocs)
            
            CALL P_Change_Decomposition("N1-N2",ntap,vpx,vpy,vpz
     &           ,nstart_1,nend_1,nstart_2,nend_2,node
     &           ,nprocs)
            CALL Get_vp_without_vc(vpx,vpy,vpz,oc,nstart_cm,nend_cm
     &           ,nprot,protl,mass,tmass,vcax,vcay,vcaz)

#endif
*=======================================================================
*---- Change coordinates of the external loop (ex) ---------------------
*=======================================================================

            IF(ptr_ex .LT. ntap) THEN
               CALL change_coord_ex(nstart_ex,nlocal_ex,nstart_cme
     &              ,nlocal_cme,protl,co,oc,lx,ly,lz,xp0,yp0,zp0,xpa,ypa
     &              ,zpa,xpcma,ypcma,zpcma)
            END IF
            IF(thermos) CALL verlet_free_eta(neta,eta,etap,tm)
            
#if defined PARALLEL
            CALL P_expand_r8x3(xp0,yp0,zp0,nstart_1,nend_1,nlocal_1,node
     &           ,nprocs,1)
            IF(ptr_ex .LT. ntap) THEN
               CALL P_expand_r8x3(xp0,yp0,zp0,nstart_ex0
     &              ,nend_ex0,nlocal_ex0,node,nprocs,ptr_ex)
            END IF
#endif

*=======================================================================
*---------- recomputes group position and scaled groups ---------------
*=======================================================================
                  
            CALL appbou(xp0,yp0,zp0,xpg,ypg,zpg,pmass,1,ngrp,grppt)
            CALL inicmp(ss_index,xp0,yp0,zp0,xpcm,ypcm,zpcm,mass,nprot
     &           ,protl)
            CALL change_frame(co,oc,-1,ntap,xp0,yp0,zp0,xpa,ypa,zpa)
            CALL change_frame(co,oc,-1,ngrp,xpg,ypg,zpg,xpga,ypga,zpga)
            CALL change_frame(co,oc,-1,nprot,xpcm,ypcm,zpcm,xpcma,ypcma
     &           ,zpcma)

*=======================================================================
*---        Computes M-forces at new coordinates                     ---
*=======================================================================
            
            IF(polar) THEN
c$$$         xpgg=zp0(1)
c$$$         bin=0.0001D0
c$$$         count=0
c$$$         DO i=-1,1
c$$$            count=count+1
c$$$            zp0(1)=xpgg+bin*i
c$$$            CALL inicmp(ss_index,xp0,yp0,zp0,xpcm,ypcm,zpcm,mass,nprot
c$$$     &           ,protl)
c$$$            CALL appbou(xp0,yp0,zp0,xpg,ypg,zpg,pmass,1,ngrp,grppt)
c$$$
c$$$            CALL change_frame(co,oc,-1,ntap,xp0,yp0,zp0,xpa,ypa,zpa)
c$$$            CALL change_frame(co,oc,-1,ngrp,xpg,ypg,zpg,xpga,ypga,zpga)
c$$$            CALL change_frame(co,oc,-1,nprot,xpcm,ypcm,zpcm,xpcma,ypcma
c$$$     &           ,zpcma)
c$$$
c$$$            CALL Polarization_Forces(fscnstr_slt,fscnstr_slv,ntap
c$$$     &           ,chrge,ss_index,plrzbij,kpol_inp,lbohr,lpnbd,alnbd
c$$$     &           ,betb,nbtype,protl,anprot,annpro,anpoint,polar
c$$$     &           ,pol_type,Ext_ef,xp0,yp0,zp0,xpa,ypa,zpa,xpga,ypga,zpga
c$$$     &           ,xpcma,ypcma,zpcma,co,oc,pmass,mapnl,fpx_m,fpy_m
c$$$     &           ,fpz_m,nnlpp0,npp_m,worka,cpu_h,ncpu_h,nstart_h
c$$$     &           ,nend_h,nlocal_h,nstart_ah,nend_ah,nlocal_ah,node
c$$$     &           ,nodex,nodey,nodez,ictxt,npy,npz,descQ,nprocs,ncube
c$$$     &           ,fstep,tag_bndg,ngrp,grppt,ingrpp,ingrp,ingrp_x
c$$$     &           ,errmsg,iret,nstep,unitc,efact,pi,alphal,U_conf
c$$$     &           ,U_ele,Udirect,Urecip,Uind,uself_dip,uself,Ugrp
c$$$     &           ,U_Thole,ucop_m,ucos_m,ucosp_m,ucnp_m,ucns_m,ucnsp_m
c$$$     &           ,U_solv,Old_dipoles,Fixt_Dipoles,mesos,Mesos_Rho
c$$$     &           ,volume,Polarization_Model)
c$$$
c$$$            energ(count)=Udirect+Ugrp+uself+uself_dip+Uind
c$$$     &           +Urecip+U_Thole+U_conf
c$$$            deriv(count)=fpz_m(1)
c$$$            WRITE(*,*) energ(count)
c$$$         END DO
c$$$         WRITE(*,*) 'Analitical:         ',deriv(2)
c$$$         WRITE(*,*) 'Finite Differences: ',-0.5D0*(energ(3)-energ(1))
c$$$     &        /bin
c$$$         WRITE(*,*) DABS(deriv(2)+0.5D0*(energ(3)-energ(1))/bin)
c$$$         WRITE(*,*) 'Derivative end'
c$$$         STOP
               CALL   Polarization_Forces(fscnstr_slt,fscnstr_slv,ntap
     &              ,chrge,ss_index,plrzbij,kpol_inp,lbohr,lpnbd,alnbd
     &              ,betb,nbtype,protl,anprot,annpro,anpoint,polar
     &              ,pol_type,Ext_ef,xp0,yp0,zp0,xpa,ypa,zpa,xpga,ypga
     &              ,zpga,xpcma,ypcma,zpcma,co,oc,pmass,mapnl,fpx_m
     &              ,fpy_m,fpz_m,worka,cpu_h,ncpu_h
     &              ,nstart_h,nend_h,nlocal_h,nstart_ah,nend_ah
     &              ,nlocal_ah,node,nodex,nodey,nodez,ictxt,npy,npz
     &              ,descQ,nprocs,ncube,fstep,tag_bndg,ngrp,grppt,ingrpp
     &              ,ingrp,ingrp_x,errmsg,iret,nstep,unitc,efact,pi
     &              ,alphal,U_conf,U_ele,Udirect,Urecip,Uind,uself_dip
     &              ,uself,Ugrp,U_Thole,ucop_m,ucos_m,ucosp_m,ucnp_m
     &              ,ucns_m,ucnsp_m,U_solv,Old_Dipoles,Fixt_Dipoles
     &              ,mesos,Mesos_Rho,volume,Polarization_Model)
c$$$               fpx(1)=DSQRT(SUM(fpx_m(1:ntap)**2+fpy_m(1:ntap)**2
c$$$     &              +fpz_m(1:ntap)**2))
c$$$               WRITE(95,'(f12.4,6e18.8)')  time*DBLE(ninner)
c$$$     &              /DBLE(mrespa*lrespa),fpx(1)

               eer_m=Urecip
            ELSE
               CALL mts_forces(rshell,xpa,ypa,zpa,xpga,ypga,zpga,xpcma
     &              ,ypcma,zpcma,mapnl,mapdn,nmapdn,ucns_m,ucos_m,virs_m
     &              ,virsp_m,ucnp_m,ucop_m,ucnsp_m,ucosp_m,fpx_m,fpy_m
     &              ,fpz_m,stressd_m,worka,cpu_h
     &              ,ncpu_h,nstart_h,nend_h,nstart_ah,nend_ah,nlocal_ah
     &              ,node,nprocs,ncube,P_dyn_update_shell)
#if defined PARALLEL
               IF(pme) THEN
                  IF(rshk .NE. rshell) THEN
                     CALL P_Reduce_Forces(fpx_m,fpy_m,fpz_m)
                  END IF
               ELSE
                  CALL P_Reduce_Forces(fpx_m,fpy_m,fpz_m)
               END IF
#endif
               IF(abmd .AND. (associate .OR. dissociate)) CALL
     &              comp_abmd_fdiss(abmd_dir,rspset,diss_list,spring,xpa
     &              ,ypa,zpa,co,fpx_m,fpy_m,fpz_m,uumb,gr)
               
               IF(clewld) THEN
                  CALL mts_furier(node,nodex,nodey,nodez,ictxt,npy,npz
     &                 ,descQ,nprocs,ncube,nstart_1,nend_1,nlocal_1
     &                 ,nstart_2,nend_2,nlocal_2,xp0,yp0,zp0,xpa,ypa,zpa
     &                 ,xpcma,ypcma,zpcma,urcsp_m,urcs_m,urcp_m,virsp_m
     &                 ,virs_m,virp_m,fpx_m,fpy_m,fpz_m,phi,fsin14
     &                 ,fsbend,fsbond,fscnstr_slt,fscnstr_slv
     &                 ,coul_bnd_slt,coul_bnd_slv,rshell,rshk,eer_m
     &                 ,stressr_m,fudgec,tag_bndg)
               END IF
            END IF

*=======================================================================
*--- Compute force on the c.of.m of each molecule and reduce the -------
*--- atomic force. Do so for   m-shell forces --------------------------
*=======================================================================
            
            CALL comp_fcm(nstart_cm,nend_cm,protl,fpx_m,fpy_m,fpz_m
     &           ,fcax_m,fcay_m,fcaz_m,mass,tmass,oc)
            
*=======================================================================
*---        Scale velocities according to the thermostat velocities ----
*=======================================================================
            
*=========  Propagate exp(i L_z tm/4) ==================================
            IF(thermos) THEN
               CALL correc_exp_scale(nstart_2,nend_2,nstart_cm,nend_cm
     &              ,cpress,ss_index,etap,tm4,vcax,vcay,vcaz,vpx,vpy,vpz
     &              ,vco)
            END IF
            
*========= Propagate exp(i L_x tm/2) -==================================

            CALL correc(vpx,vpy,vpz,fpx_m,fpy_m,fpz_m,mass,nstart_2
     &           ,nend_2,tm)


            CALL correc_stress(cpress,m_skin,nprot,co,oc,vco,vcax,vcay
     &           ,vcaz,fcax_m,fcay_m,fcaz_m,stressd_m,stressr_m,volume
     &           ,press_m,press_kin,pext,tmass,masspp,tm,tm2,nstart_cm
     &           ,nend_cm,node,nprocs,ncube,rbyte)

*=======================================================================
*---        Scale velocities a second time -----------------------------
*=======================================================================
            
*=========  Propagate exp(i L_z tm/4) ==================================
            IF(thermos) THEN
               CALL correc_exp_scale(nstart_2,nend_2,nstart_cm,nend_cm
     &              ,cpress,ss_index,etap,tm4,vcax,vcay,vcaz,vpx,vpy,vpz
     &              ,vco)
            END IF
            
*=======================================================================
*---------  Correct velocities of the NETA thermostats -----------------
*=======================================================================
            
            IF(thermos) THEN
               IF(isostress  .OR. FixedAngles_Stress) THEN
                  CALL rattle_correc_co(co,dssco,cnstco,vco,masspp,iret
     &                 ,errmsg)
                  IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
               END IF
               IF(cnstpp .NE. 0) THEN
                  CALL rattle_correc(nstart_2,nend_2,tm,xp0,yp0,zp0
     &                 ,vpx,vpy,vpz,ntap,cnstp,dssp,coeffp,cnstpp,mass
     &                 ,dnit,cnst_protp,cnst_protl,mim_lim,gcpu_rt
     &                 ,iret,errmsg)
                  IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
               END IF
*--------   Compute forces on the thermostat
               CALL comp_thermos_forces(nprocs,nstart_2,nend_2,nstart_cm
     &              ,nend_cm,ndf_thermos,ss_index,co,vpx,vpy,vpz,vcax
     &              ,vcay,vcaz,mass,tmass,t,fth)
#if defined PARALLEL
               CALL P_merge_vecr8(fth,3)
#endif
*--------   Add contribution from barostat temperature
               IF(cpress) CALL comp_thermos_forces_cm(fth(1),vco,masspr)

               CALL correc_etap(neta,etap,fth,qmass,tm4)
            END IF
            
*=======================================================================
*---------  Correct velocities of c.of.m for the barostat --------------
*=======================================================================
            
*==========  Propagate exp(i L_y tm/2) =================================
            IF(cpress .AND. coupl_mol) THEN
               CALL correc_matr(tm2,co,oc,vco,gmgp,vcax,vcay,vcaz,nprot)
            END IF
            
*=======================================================================
*---------  Correct velocities of the NETA thermostats ----------------
*=======================================================================
            
            IF(thermos) THEN
*--------   Compute forces on the thermostat
               CALL comp_thermos_forces(nprocs,nstart_2,nend_2,nstart_cm
     &              ,nend_cm,ndf_thermos,ss_index,co,vpx,vpy,vpz,vcax
     &              ,vcay,vcaz,mass,tmass,t,fth)
#if defined PARALLEL
               CALL P_merge_vecr8(fth,3)
#endif
*--------   Add contribution from barostat temperature
               IF(cpress) CALL comp_thermos_forces_cm(fth(1),vco,masspr)

               CALL correc_etap(neta,etap,fth,qmass,tm4)
               CALL comp_thermos_energy(neta,ndf_thermos,t,qmass,eta
     &              ,etap,uceh,hpot,temph)
            END IF

*=======================================================================
*----------  Compute averages and do some analysis at time step M ------
*=======================================================================


            INCLUDE 'mtsmd_avg_inc.CPP.f'
            IF(abmd_native) THEN
               IF(MOD(ninner,nprint) .EQ. 0) THEN
                  fstep=time*DBLE(ninner)/DBLE(mrespa*lrespa)
                  WRITE(kprint,'('' FREE '',f10.2,2e16.8)') fstep
     &                 ,a_free*efact/1000.0D0,ts_free/1000.0D0
               END IF
            END IF

            CALL timer(vfcp_mm,tfcp_mm,elapse)
            tdelta_mm=tfcp_mm-tdelta_mm
            gcpu_mm=gcpu_mm+tdelta_mm

         END DO
         
*=======================================================================
*---     <--|  END of M-loop                                        ----
*=======================================================================
         
*=======================================================================
*---     Computes L-forces at new coordinates                       ----
*=======================================================================
         
         CALL timer(vfcp_ll,tfcp_ll,elapse)
         tdelta_ll=tfcp_ll

         rshell='l'
         IF(.NOT. polar) THEN
            CALL mts_forces(rshell,xpa,ypa,zpa,xpga,ypga,zpga,xpcma
     &           ,ypcma,zpcma,mapnl,mapdn,nmapdn,ucns_l,ucos_l,virs_l
     &           ,virsp_l,ucnp_l,ucop_l,ucnsp_l,ucosp_l,fpx_l,fpy_l
     &           ,fpz_l,stressd_l,worka,cpu_h
     &           ,ncpu_h,nstart_h,nend_h,nstart_ah,nend_ah,nlocal_ah
     &           ,node,nprocs,ncube,P_dyn_update_shell)
#if defined PARALLEL
            IF(pme) THEN
               IF(rshk .NE. rshell) THEN
                  CALL P_Reduce_Forces(fpx_l,fpy_l,fpz_l)
               END IF
            ELSE
               CALL P_Reduce_Forces(fpx_l,fpy_l,fpz_l)
            END IF
#endif         
            IF(clewld) THEN 
               CALL mts_furier(node,nodex,nodey,nodez,ictxt,npy,npz
     &              ,descQ,nprocs,ncube,nstart_1,nend_1,nlocal_1
     &              ,nstart_2,nend_2,nlocal_2,xp0,yp0,zp0,xpa,ypa,zpa
     &              ,xpcma,ypcma,zpcma,urcsp_l,urcs_l,urcp_l,virsp_l
     &              ,virs_l,virp_l,fpx_l,fpy_l,fpz_l,phi,fsin14,fsbend
     &              ,fsbond,fscnstr_slt,fscnstr_slv,coul_bnd_slt
     &              ,coul_bnd_slv,rshell,rshk,eer_l,stressr_l,fudgec
     &              ,tag_bndg)
            END IF
         END IF

*=======================================================================
*--- Compute force on the c.of.m of each molecule and reduce the -------
*--- atomic force. Do so for   l-shell forces --------------------------
*=======================================================================
         
         CALL comp_fcm(nstart_cm,nend_cm,protl,fpx_l,fpy_l,fpz_l,fcax_l
     &        ,fcay_l,fcaz_l,mass,tmass,oc)
         
*---     corrects velocities
         
         CALL correc(vpx,vpy,vpz,fpx_l,fpy_l,fpz_l,mass,nstart_2,nend_2
     &        ,tl)
         
         CALL correc_stress(cpress,l_skin,nprot,co,oc,vco,vcax,vcay
     &        ,vcaz,fcax_l,fcay_l,fcaz_l,stressd_l,stressr_l,volume
     &        ,press_l,press_kin,pext,tmass,masspp,tl,tl2,nstart_cm
     &        ,nend_cm,node,nprocs,ncube,rbyte)

         CALL timer(vfcp_ll,tfcp_ll,elapse)
         tdelta_ll=tfcp_ll-tdelta_ll
         gcpu_ll=gcpu_ll + tdelta_ll
      END DO
      
*=======================================================================
*---  <--| END of L-loop                                            ----
*=======================================================================
      
      
*=======================================================================
*---     Computes H-forces (reciprocal +direct lattice)             ----
*---               at new coordinates                               ----
*=======================================================================
      
      CALL timer(vfcp_hh,tfcp_hh,elapse)
      tdelta_hh=tfcp_hh
      rshell='h'
      IF(.NOT. polar) THEN
         CALL mts_forces(rshell,xpa,ypa,zpa,xpga,ypga,zpga,xpcma,ypcma
     &        ,zpcma,mapnl,mapdn,nmapdn,ucns_h,ucos_h,virs_h,virsp_h
     &        ,ucnp_h,ucop_h,ucnsp_h,ucosp_h,fpx_h,fpy_h,fpz_h,stressd_h
     &        ,worka,cpu_h,ncpu_h,nstart_h,nend_h
     &        ,nstart_ah,nend_ah,nlocal_ah,node,nprocs,ncube
     &        ,P_dyn_update_shell)
#if defined PARALLEL
         IF(pme) THEN
            IF(rshk .NE. rshell) THEN
               CALL P_Reduce_Forces(fpx_h,fpy_h,fpz_h)
            END IF
         ELSE
            CALL P_Reduce_Forces(fpx_h,fpy_h,fpz_h)
         END IF
#endif      
         IF(clewld) THEN 
            CALL mts_furier(node,nodex,nodey,nodez,ictxt,npy,npz,descQ
     &           ,nprocs,ncube,nstart_1,nend_1,nlocal_1,nstart_2,nend_2
     &           ,nlocal_2,xp0,yp0,zp0,xpa,ypa,zpa,xpcma,ypcma,zpcma
     &           ,urcsp_h,urcs_h,urcp_h,virsp_h,virs_h,virp_h,fpx_h
     &           ,fpy_h,fpz_h,phi,fsin14,fsbend,fsbond,fscnstr_slt
     &           ,fscnstr_slv,coul_bnd_slt,coul_bnd_slv,rshell,rshk
     &           ,eer_h,stressr_h,fudgec,tag_bndg)
         END IF
      END IF
      
*=======================================================================
*--- Compute force on the c.of.m of each molecule and reduce the -------
*--- atomic force. Do so for   h-shell forces --------------------------
*=======================================================================
      
      CALL comp_fcm(nstart_cm,nend_cm,protl,fpx_h,fpy_h,fpz_h,fcax_h
     &     ,fcay_h,fcaz_h,mass,tmass,oc)
      
*=======================================================================
*---  Corrects velocities                                            ---
*=======================================================================
      
      CALL correc(vpx,vpy,vpz,fpx_h,fpy_h,fpz_h,mass,nstart_2,nend_2
     &     ,time)

      CALL correc_stress(cpress,h_skin,nprot,co,oc,vco,vcax,vcay
     &     ,vcaz,fcax_h,fcay_h,fcaz_h,stressd_h,stressr_h,volume
     &     ,press_h,press_kin,pext,tmass,masspp,time,time2,nstart_cm
     &     ,nend_cm,node,nprocs,ncube,rbyte)
      
*=======================================================================
*---  Interaction with a stochastic bath                            ---
*=======================================================================
      
      IF(landersen) THEN
         CALL comp_vel_labframe(vpx,vpy,vpz,vpx,vpy,vpz,co,nprot,protl
     &        ,vcax,vcay,vcaz)
         CALL collision(ntap,vpx,vpy,vpz,mass,nutime,t,time)
         CALL comp_vcm(vpx,vpy,vpz,oc,nprot,protl,mass,tmass,vcax,vcay
     &        ,vcaz)
      END IF


*=======================================================================
*----------  Dump restart file and do tests at timestep H --------------
*=======================================================================

      INCLUDE 'mtsmd_dump_inc.CPP.f'

*=======================================================================
*--------         MD LOOP ENDS   !!                    -----------------
*=======================================================================

      CALL timer(vfcp_hh,tfcp_hh,elapse)
      tdelta_hh=tfcp_hh-tdelta_hh
      gcpu_hh=gcpu_hh + tdelta_hh
      IF(nstep .LT. nrject) GOTO 100
      IF(nstep .LT. maxstp) GOTO 100

      IF(annealing) THEN
         grad_max=-1.0D0
         IF(.NOT. ASSOCIATED(fpx)) THEN
            ALLOCATE(fpx(ntap),fpy(ntap),fpz(ntap))
         END IF
         DO i=1,ntap
            fpx(i)=fpx_n0(i)+fpx_n1(i)+fpx_m(i)+fpx_l(i)+fpx_h(i)
            fpy(i)=fpy_n0(i)+fpy_n1(i)+fpy_m(i)+fpy_l(i)+fpy_h(i)
            fpz(i)=fpz_n0(i)+fpz_n1(i)+fpz_m(i)+fpz_l(i)+fpz_h(i)
            IF(DABS(fpx(i)) .GT. grad_max) grad_max=DABS(fpx(i))
            IF(DABS(fpy(i)) .GT. grad_max) grad_max=DABS(fpy(i))
            IF(DABS(fpz(i)) .GT. grad_max) grad_max=DABS(fpz(i))
         END DO
         WRITE(kprint,20000) DABS(grad_max)
         CALL prtfrc(kprint,ngrp,grppt,nres,M1,prsymb,beta,xp0,yp0,zp0
     &        ,fpx,fpy,fpz)
      END IF

*=======================================================================
*     Write timing
*=======================================================================

      CALL timer(vfcp,tfcp,elapse)
      gcpu=-gcpu + tfcp
      elaps= -elaps + elapse
      
      write(kprint,*)
      write(kprint,60030)
      write(kprint,17000) gcpu
      write(kprint,18000) elaps
      write(kprint,60200) gcpu/time_fs
      write(kprint,60210) gcpu_inn/time_fs
      write(kprint,60220) gcpu_mm/time_fs
      write(kprint,60230) gcpu_ll/time_fs
      write(kprint,60240) gcpu_hh/time_fs
      write(kprint,60260) gcpu_rt/time_fs
      write(kprint,60270) gcpu_sh/time_fs
      write(kprint,60300) elaps/time_fs
      write(kprint,60030)

      IF(wrtgyr) THEN
         CLOSE(kgyr)
      END IF

*================= END OF EXECUTABLE STATEMENTS ========================

1000  FORMAT(/ / / /20('*'),'  M o l e c u l a r   T o p o l o g y  ',
     &     20('*')/ / / / /)
1100  FORMAT('*',13(' '),
     &     ' I n p u t   O p e r a t i o n s   C o m p l e t e d ',
     &     12(' '),'*'/'*',78(' '),'*')
1110  FORMAT('*',13(' '),
     &     '   M . D.   I n t e r m e d i a t e   R e s u l t s  ',
     &     12(' '),'*'/'*',78(' '),'*')
1200  FORMAT(80('*'))
1300  FORMAT('*',78(' '),'*')
2000  FORMAT(/'<------ Linked cell list indexing updated -------->'/)
13000 FORMAT(/22x,' Temperature has been rescaled ',i5,' times '/)
14000 FORMAT(/22x,'Adjusting bond length to Force Field.'/
     &     22x,'     This will take a while...'/ /) 
14010 FORMAT(/22x,'Setting up Parallel Decompositions.'/
     &       22x,'     This will take a while...'/ /) 
20000 FORMAT(
     &     21x,'***************************************'/
     &     21x,'*        GradMax    ',e12.5,'      *'/
     &     21x,'***************************************')
70200 FORMAT('<------ Dumping Restart File ------->'/)
70300 FORMAT(/ /'<------ Reading Restart File ------->')
70400 FORMAT('<------ Restart File Read in ------->'/ /)
70100 FORMAT('Velocities have been rescaled ---->'/)
      
17000 FORMAT(  15x,' Total cpu time for the run       = ',f10.3)
18000 FORMAT(  15x,' Total elapsed time for the run   = ',f10.3)
60200 FORMAT(  15x,' Averaged time per femtosecond    = ',3x,f7.3)
60210 FORMAT(  15x,' Averaged time in inner loop      = ',3x,f7.3)
60220 FORMAT(  15x,' Averaged time in m loop          = ',3x,f7.3)
60230 FORMAT(  15x,' Averaged time in l loop          = ',3x,f7.3)
60240 FORMAT(  15x,' Averaged time in h loop          = ',3x,f7.3)
60250 FORMAT(  15x,' Averaged time inside the m loop  = ',3x,f7.3)
60260 FORMAT(  15x,' Averaged time to rattle vel.     = ',3x,f7.3)
60270 FORMAT(  15x,' Averaged time to rattle coord.   = ',3x,f7.3)
60300 FORMAT(  15x,' Averaged elapsed per femtosecond = ',3x,f7.3/ /)
15000 FORMAT(/ / / 10x,'* * * * r - R E S P A  i s  O N  * * * *'  / /)
15500 FORMAT(/   10x,'======= timing ========================='  /)
16011 FORMAT(/5x,'CPUtime for update          =',f10.3)
15011 FORMAT(/5x,'CPUtime for linked cell indexing =',f10.3)
16001 FORMAT(/5x,'CPUtime for h-contribution: RECP =',f8.2,
     &     ' DIR =',f7.3,' TOT =',f7.3)
16002 FORMAT(/5x,'CPUtime for l-contribution: RECP =',f8.2,
     &     ' DIR =',f7.3,' TOT =',f7.3)
16003 FORMAT(/5x,'CPUtime for m-contribution: RECP =',f8.2,
     &     ' DIR =',f7.3,' TOT =',f7.3)
16005 FORMAT(/5x,'CPUtime for Communications: H-Sh =',f9.4,
     &     ' L-Sh = ',f9.4,' M-Sh =',f9.4)
14004 FORMAT(/5x,'CPUtime for n1-contribution  =',f9.4) 
16004 FORMAT(/5x,'CPUtime for n0-contribution  =',f9.4) 
10067 FORMAT(/5x,'THEORIC SPEED UP FOR NON BONDED PART =',f8.2)
10068 FORMAT(/5x,'OVERALL THEORIC SPEED UP =',f8.2/)
80033 format(/5x,'Expected CPU time for the RUN:',I4,
     &     ' hours and ',I2, ' min',/,  
     &     /5x,' Expected average time per M step:',f8.2,' sec.'/  
     &     /5x,' Expected average time per femto :',f8.2,' sec.'/)  
80035 FORMAT(/5x,'***** Communication time direct forces ******'/
     &     /5x,'     Shell h : ',f9.4/
     &     /5x,'     Shell l : ',f9.4/
     &     /5x,'     Shell m : ',f9.4/)
60030 FORMAT(/10x,'==========================================='/)
78410 FORMAT
     &     (/ /' *******ERROR: MAXT for PME too small. INCREASE MAXT.'/
     &     /)
70500 FORMAT(/ /15x,' <-------- Time Limit Reached  ------->'/ /)
70700 FORMAT(/ /15x,' Program Stops Smoothly. Restart Dumped'/ /)
70120 FORMAT('Velocities of the barostat have been rescaled   ---->'/)
70130 FORMAT('Velocities of the thermostats have been rescaled   ---->'/
     &     )
20900 FORMAT(/ /' ****** WARNING ! WARNING ! WARNING ***************' /
     &     /' ******   drift_remove ON           ***************' /
     &     /' ****** WARNING ! WARNING ! WARNING ***************' / /)
10072 FORMAT(/ /' --- Total energy remove per particle = ', f12.4,
     &     ' (KJ/mole)'/
     &     ' --- Number of dirty scaling          = ', I10/
     &     ' --- Frequency of scaling             = ',f10.3,
     &     ' 1/ps '/)
10977 FORMAT(/ /'*******WARNING: NO COFACTOR ATOMS SELECTED '/ 
     &     ' NCOFACTOR IS SET TO ZERO AND NO FIELD I COMPUTED'/ /)
80000 FORMAT('REMARK   Rigid body fit on CA atoms')
80100 FORMAT('REMARK   Rigid body fit on heavy atoms')
80101 FORMAT('REMARK   Rigid body fit on backbone atoms')
93410 FORMAT(/ /
     &     21x,'***************************************'/
     &     21x,'*                                     *'/
     &     21x,'*       MD run on ',i3,' CPUs            *'/
     &     21x,'*                                     *'/
     &     21x,'*                                     *'/
     &     21x,'***************************************'/ /)
89000 FORMAT(/ /
     &     21x,'***************************************'/
     &     21x,'*                                     *'/
     &     21x,'*     Force Field Parameters for      *'/
     &     21x,'*       the current molecules         *'/
     &     21x,'*                                     *'/
     &     21x,'*                                     *'/
     &     21x,'***************************************'/ /)
89500 FORMAT(/ /21x,'***************************************'/ /)
90000 FORMAT('BOND')
91000 FORMAT('BENDING')
92000 FORMAT('TORSION PROPER')
93000 FORMAT('TORSION IMPROPER')
94000 FORMAT('END')
95000 FORMAT('#   Constrained stretchings follow')
96000 FORMAT('NONBONDED MIXRULE')
97000 FORMAT('#   Warning! Mixrules always assumed')
98000 FORMAT(/38('>'),38('<')/
     &     ' Warning number of groups in the restart file is ',
     &     'different from current'/
     &     5x,' Restart group number is ',i7,'  current is  ',i7/38('>')
     &     ,38('<')/)
99000 FORMAT(/38('>'),38('<')/
     &     ' Warning number of molecules in the restart file is ',
     &     'different from current.'/
     &     5x,' Restart molecule number is ',i7,'  current is
     &     ',i7/38('>'),38('<')/)
87000 FORMAT(/ /
     &     21x,'***************************************'/
     &     21x,'*                                     *'/
     &     21x,'*        Polarization Model is        *'/
     &     21x,'*                                     *'/
     &     21x,'*           G a u s s i a n           *'/
     &     21x,'*                                     *'/
     &     21x,'***************************************'/ /)
88000 FORMAT(/ /
     &     21x,'***************************************'/
     &     21x,'*                                     *'/
     &     21x,'*        Polarization Model is        *'/
     &     21x,'*                                     *'/
     &     21x,'*              F u l l                *'/
     &     21x,'*                                     *'/
     &     21x,'***************************************'/ /)
81000 FORMAT(/ /10x,
     &     '<------ Initializing thermostat coordinates ------->')
82000 FORMAT(/ /10x,
     &     '<------- Initializing barostat coordinates -------->')
      RETURN
      END
