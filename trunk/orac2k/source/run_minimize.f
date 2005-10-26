#define _FP_        ip_fp
#define _RP_        ip_rp
#define _WA_        ip_wa
#define _IWA_       ip_iwa
#define _NBD_       ip_nbd
      SUBROUTINE run_minimize(mapnl,xp0,yp0,zp0,xpg,ypg,zpg,eta,xpcm
     &     ,ypcm,zpcm)

************************************************************************
*                                                                      *
*     MTSMD is the driver of the MD run when multiple time scales      *
*     methods are used. The integration algorithm is r-RESPA           *
*     of order o(dt^3) (Tuckermann et al. JCP 97 1990 (1992))          *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi                                        *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*   RUN_MINIMIZE externals:    	       	       	       	       	       *  
*       add_energies appbou change_frame chkewl comp_dynamic_mat       *
*       comp_molmass comp_vel_labframe copy_protl dumprs	       *
*       erfc_spline_init erf_corr_cutoff fft_pme_init fndgrp	       *
*       gen_stress_tag_bnd get_total_energy inicmp int_corr_erf_spline *
*       lc_index lc_list linmin_total matinv mts_forces		       *
*       plotc prtacc prtat prtba prtbnd prtcn			       *
*       prtfrc prtit prtite_min prtpt prtsq readrs		       *
*       timer tr_inbox write_pot_bond write_pot_nbond xerror	       *
*       zero zero3x3 zeroa					       *
*                                                                      *
************************************************************************


*======================= DECLARATIONS ==================================

      IMPLICIT none
      
      INCLUDE 'parst.h'

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
      INCLUDE 'unit.h'
      INCLUDE 'pme.h'
      INCLUDE 'lc_list.h'
      INCLUDE 'parallel.h'

*-------------------- VARIABLES FOR L-BFGS-B  ---------------------------

      INTEGER  nmax,mmax,lenwa
      PARAMETER( nmax=3*mb,mmax = 17)
      PARAMETER(lenwa = 2*mmax*nmax +  4*nmax + 11*mmax*mmax + 8*mmax)
      REAL*8  dsave(29),factr,low,up
      REAL(8), DIMENSION(:), POINTER :: rp,fp,wa
      INTEGER isave(44),iprint,kbfgs,len1,len2
      INTEGER, DIMENSION(:), POINTER :: nbd,iwa

      CHARACTER*60 task, csave
      LOGICAL lsave(4),bfgs_restart
      INTEGER bfgs_nf,nf_press,bfgs_nf_old,bfgs_m_old
      CHARACTER*80 bfgs_file

*-------------------- DEFINITION OF AN EXTERNAL FUNCTION ---------------

      REAL*8   ddot,dnrm2
      LOGICAL near0,ok,exist
      EXTERNAL  near0,ddot,dnrm2,f1dim_der

*-------------------- LOCAL VARIABLES ----------------------------------

      CHARACTER*80 errmsg
      character*1  rshk

      INTEGER    kp,kt,nato_slt,iter_avg

      INTEGER i,j,nsstt,tstep,nstep,iret
      INTEGER flag
      LOGICAL lupdate
      REAL*8  ucns,ucos,ucnsp,ucosp,ucnp,ucop,fudgec,fstep
      REAl*8 urcs,urcp,urcsp,eer,stressr(3,3),stressd(3,3),stressd_p(3,3
     &     ),stress_conf(3,3),stress_tot(3,3),prt(3,3)
     &     ,press_conf,press_kin,errca,errhe,errbc
     &     ,erral,drpca,drpbc,drphe,drpal

      REAL*8  temp,elapse,puconf,pucoul,puhyd,pubnd,ubend,uptors,uitors
     &     ,uconf,ucoul,ureal,urecp,ucek,pucek,tempt,fsrtal,tempr,temppr
     &     ,tcm,rcm,elaps,rms,unb14,cnb14,fsbond,fsbend,fsin14,purecp
     &     ,vfcp,tfcp,uumb,upconf,upcoul,gcpu,temph,uceh,hpot
     &     ,ubond,aux,hstep,uslvbon,uslvben,ucepr,gcpu_u
      REAL*8  gr,pueng,gra
      INTEGER navg
      PARAMETER (navg = 60)

      REAL*8  sumarray(navg),ssmarray(navg)
      REAL*8 sum_econf,sum_ecoul,sum_enbnd,sum_etotpot,sum_tote,sum_ucek
     &     ,sum_temp,sum_tempt,sum_tempr,sum_temppr,sum_temph,sum_rms
     &     ,sum_pecek,sum_pehyd,sum_peconf,sum_pecoul,sum_percip
     &     ,sum_enb14,sum_ebend,sum_ebond,sum_eitor,sum_eptor,sum_pnbd
     &     ,sum_pebnd,sum_pepot,sum_ptote,sum_gr,sum_epcoul,sum_epconf
     &     ,sum_presst,sum_press,sum_pressc,sum_pressk,sum_st(3,3)
     &     ,sum_co(3,3),sum_volume,sum_pv,sum_temppra,sum_eslvint
     &     ,sum_eebond,sum_eebend,sum_eeptors,sum_eeitors

      REAL*8 ssm_econf,ssm_ecoul,ssm_enbnd,ssm_etotpot,ssm_tote,ssm_ucek
     &     ,ssm_temp,ssm_tempt,ssm_tempr,ssm_temppr,ssm_temph,ssm_rms
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
     &     (sum_temph, sumarray(11)), (sum_rms, sumarray(12)),
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
     &     (ssm_temph, ssmarray(11)), (ssm_rms, ssmarray(12)),
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

      REAL*8  temppra,energy

*----------- ARRAYS USED BY COMP_DYNAMIC_MAT ---------------------------
      
      REAL*8 fpx2(m1,cheb_order),fpy2(m1,cheb_order),fpz2(m1,cheb_order)
     &     ,d_mat(n_mat*(n_mat+1)/2),wk(n_mat*(n_mat+1)/2+n_mat)
     &     ,eigvl(n_mat),eigvc(n_mat,n_mat)
      INTEGER mad,mbd
      
*----------- LOCAL WORK ARRAYS FOR THE RUN -----------------------------

      INTEGER ngrp_old,nprot_old,nind(2),indxyz,ind_a

      INTEGER, DIMENSION(:), POINTER ::  indxi,indxj,indxk

      REAL*8  gpx(mb),gpy(mb),gpz(mb),hpx(mb),hpy(mb),hpz(mb),fpx(mb)
     &     ,fpy(mb),fpz(mb),fpx1(mb),fpy1(mb),fpz1(mb),vpx(mb),vpy(mb)
     &     ,vpz(mb),vpx1(mb),vpy1(mb),vpz1(mb)
      REAL*8  etap(hoov),vh1(hoov),vco(3,3)
      REAL*8  zz1(3,3),lzz,lzz1
      REAL*8  xpo(mb),ypo(mb),zpo(mb),co2(3,3),oc2(3,3)
      REAL*8  xp1(mb),yp1(mb),zp1(mb),xpa(mb),ypa(mb),zpa(mb),xpga(m11)
     &     ,ypga(m11),zpga(m11),xpcma(npm),ypcma(npm),zpcma(npm)
     &     ,tmass(numpr),tmassb(numpr)
      INTEGER dthree
      INTEGER numatoms,igrn,krdf(maxint*g1)
      INTEGER mapdn(2,mf),nmapdn(mf),worka(m1,m10)
      INTEGER npp,npp_u,mpp_a,M_get_length
      INTEGER cnstp(2,mbs),cnstpp,cnstpp_slv,cnst_protl(2*m1)
     &     ,cnst_protl_1(2*m1),cnst_protl_ex(2*m1),mapnl_save(m8),nmapnl
     &     ,cnst_protp,cnst_protp_1,cnst_protp_ex

*--- DYNAM is a scratch common block: here to save storage; not passed
*    to any of the external

      COMMON /dynam/ rp,fp,wa,dsave,xp1,yp1,zp1,xpa,ypa
     &     ,zpa,xpcma,ypcma,zpcma,xpga,ypga,zpga,gpx,gpy,gpz,hpx,hpy,hpz
     &     ,fpx,fpy,fpz,fpx1,fpy1,fpz1,vpx,vpy,vpz,vpx1,vpy1,vpz1,fpx2
     &     ,fpy2,fpz2,d_mat,wk,eigvl,eigvc,mapdn,nmapdn,worka,cnstp
     &     ,nbd,iwa,isave,cnstpp,cnstpp_slv,cnst_protl
     &     ,cnst_protl_1,cnst_protl_ex,mapnl_save,nmapnl,cnst_protp
     &     ,cnst_protp_1,cnst_protp_ex,krdf

*     Phony forces and energies for neighbor list

      REAL*8  ucns_p,ucos_p,virs_p,virsp_p,ucnp_p,ucop_p,ucnsp_p,ucosp_p
     &     ,fpx_p(1),fpy_p(1),fpz_p(1),conf_bnd_slt,coul_bnd_slt
     &     ,conf_bnd_slv,coul_bnd_slv,self_slt,self_slv,uslvtor,uslvitor
     &     ,conf_bnd_slt_n1,coul_bnd_slt_n1,conf_bnd_slv_n1
     &     ,coul_bnd_slv_n1
      
      LOGICAL dpress,dhoover,grflag
      INTEGER tag_bndg(m2)
      REAL*8  grad_max
      REAL*8  work(mspline),utotal,gg,dgg,gamma,eps2,eps,fret,ffwork(2)
     &     ,utotal_old

      INTEGER count

*==================== EXECUTABLE STATEMENTS ============================



*===  Check if the dimension of the work array are sufficient 


      IF(mb.LT.ntap) THEN
         errmsg=' While in MTSMD: PARAMETER MB dimensions the work'
     &        / /' arrays is sufficient. Abort. '
         CALL xerror(errmsg,80,1,2)
         STOP
      END IF

*=======================================================================
*----- Initialize some stuff -------------------------------------------
*=======================================================================

*===  set few variable to zero

      flag=nflag(1)
      hstep=0.3D0
      urcs=0.0D0
      urcp=0.0D0
      urcsp=0.0D0
      eer=0.0D0
      rcutl=0.0D0
      rtoll=0.0D0
      ucek=0.0D0
      pucek=0.0D0
      uceh=0.0D0
      ucepr=0.0D0
      fsrtal=0.0D0
      t=0.0D0
      CALL zero3x3(stressr)
      CALL zero3x3(stressd)
      uumb=0.0D0
      conf_bnd_slv=0.0D0
      conf_bnd_slt=0.0D0
      tstep=0
      fudgec=1.d0-fudge
      nsstt=0
      iret=0
      errmsg=' '
      nato_slt=ntap-nmol*nato_slv
      iter_avg=0
      eps=eps_energy
      eps2=eps**2
      cnstpp=0
      cnst_protl(1)=0
      cnst_protp=0
      iprint=0
      task='START'
      press_kin=0.0D0
      temppra=0.0D0
      indxyz=indmax

*=======================================================================
*---- Allocate memory for neighbor list                             ----
*=======================================================================

      IF(linked_cell) THEN
         indxyz=(2*(ncx-1)+1)*(2*(ncy-1)+1)*(2*(ncz-1)+1)
         ALLOCATE(indxi(indxyz),indxj(indxyz),indxk(indxyz))
      END IF

      CALL comp_molmass(nprot,protl,mass,tmass)

*=======================================================================
*----- Copy old CO and OC matrix to temporary arrays -------------------
*=======================================================================

      IF(change_cell) THEN
         CALL dcopy(9,co,1,co2,1)
         CALL dcopy(9,oc,1,oc2,1)
      END IF

*=======================================================================
*----- Read input from restart file if needed --------------------------
*=======================================================================
      
      kp=0
      kt=-1
      IF(flag.GT.0) THEN
         WRITE(kprint,70300)
         IF(restart_old) THEN
            CALL readrs_old(kdump_in,restart_in,nstep,temp,ntap,ngrp_old
     &           ,nprot_old,xp0,yp0,zp0,vpx,vpy,vpz,vpx1,vpy1,vpz1,vpx1
     &           ,vpy1,vpz1,vpx1,vpy1,vpz1,vpx1,vpy1,vpz1,eta,etap,vpx1
     &           ,dthree,dhoover,sumarray,ssmarray,navg,navg,co,oc
     &           ,dpress,vco,vpx1,vpy1,vpz1,grflag,krdf,igrn,maxint*g1)
         ELSE
            CALL readrs(kdump_in,restart_in,nstep,temp,ntap,ngrp_old
     &           ,nprot_old,xp0,yp0,zp0,vpx,vpy,vpz,eta,etap,dthree
     &           ,dhoover,sumarray,ssmarray,navg,navg,co,oc,dpress,vco
     &           ,grflag,krdf,igrn,maxint*g1)
         END IF
         WRITE(kprint,70400)
         temp=0.0D0
      END IF
      CALL zeroa(vpx,vpy,vpz,ntap,1)
      CALL zero3x3(vco)
      CALL zero0(eta,hoov)
      CALL zero0(etap,hoov)


      IF(flag .EQ. 1 .AND. nstep .GE. maxstp ) THEN
         WRITE(errmsg,'(''In MTSMD: Step number in restart = '',
     &        i6,'' while in input = '',i6,'' Increase TIME.'')')
     &        nstep,maxstp
         CALL xerror(errmsg,80,1,2)
      END IF

*=======================================================================
*---  Write a Banner Page                                            ---
*=======================================================================

      WRITE(kprint,'(/ /)')
      WRITE(kprint,1200)
      WRITE(kprint,1300)
      WRITE(kprint,1100)
      WRITE(kprint,1100)
      WRITE(kprint,1300)
      WRITE(kprint,1200)
      WRITE(kprint,'(/ /)')

      CALL zero(sumarray,ssmarray,navg,navg)
      nstep=0

      IF(coupl_atm .OR. coupl_grp) THEN
         CALL copy_protl(protl,protlb,nprot,nprotb)
         CALL comp_molmass(nprot,protl,mass,tmassb)
         nprot=ntap
         count=0
         DO i=1,ntap
            protl(1+count)=1
            protl(2+count)=i
            count=count+2
         END DO
         IF(nprot .GT. npm) THEN
            errmsg=
     &           'Number of molecules exceeds physical dimensions. '
     &           / /'Increase _N_PROT_MAX_ to '
            WRITE(errmsg(75:80),'(i6)') nprot
            CALL xerror(errmsg,80,1,2)            
         END IF
         CALL comp_molmass(nprot,protl,mass,tmass)
      END IF

*=======================================================================
*----- Calculate solute center of mass coordinates and velocities ------
*=======================================================================

      CALL inicmp(ss_index,xp0,yp0,zp0,xpcm,ypcm,zpcm,mass,nprot,protl)

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
      END IF
      IF(nflag(1) .NE. 0 .AND. nprot_old .NE. nprot ) THEN
         WRITE(kprint,99000) nprot_old,nprot
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
*--- Change frame to get xpa, ypa, zpa etc in box fractions ------------
*=======================================================================

      CALL change_frame(co,oc,-1,ntap,xp0,yp0,zp0,xpa,ypa,zpa)
      CALL change_frame(co,oc,-1,ngrp,xpg,ypg,zpg,xpga,ypga,zpga)
      CALL change_frame(co,oc,-1,nprot,xpcm,ypcm,zpcm,xpcma,ypcma,zpcma)

*=======================================================================
*----- Get for each node its atomic decomposition ----------------------
*=======================================================================

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
     &     ,ntap,ngrp,grppt,worka,xpa,ypa,zpa
     &     ,xpga,ypga,zpga)


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
*----- Change the cell according to input commands ---------------------
*=======================================================================

      IF(change_cell) THEN
         CALL change_frame(co,oc,-1,nprot,xpcm,ypcm,zpcm,xpcma,ypcma
     &        ,zpcma)
         CALL change_origin(1,nprot,protl,xp0,yp0,zp0,xpo,ypo,zpo,xpcma
     &        ,ypcma,zpcma,co)

         CALL dcopy(9,co2,1,co,1)
         CALL dcopy(9,oc2,1,oc,1)

         CALL change_origin(-1,nprot,protl,xp0,yp0,zp0,xpo,ypo,zpo,xpcma
     &        ,ypcma,zpcma,co)
         CALL change_origin(1,nprot,protl,xp0,yp0,zp0,xpo,ypo,zpo,xpcma
     &        ,ypcma,zpcma,co)
         CALL change_frame(co,oc,1,nprot,xpcma,ypcma,zpcma,xpcm,ypcm
     &        ,zpcm)
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
     &           ,erf_arr_corr,work)
            END IF
            CALL Pme_init(node,nodex,nodey,nodez,npy,npz,ictxt,descQ
     &           ,fftable,nfft1,nfft2,nfft3,nfft3_local,iret,errmsg)
            IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
         ELSE
            aux=(rtolh+rcuth)**2
            CALL chkewl(oc,aux,alphal,rkcut,volume)
            rshk=shell_pme
            IF(erf_corr) THEN  
              CALL int_corr_erf_spline(rlew,ruew,nbinew,alphal,rkcut
     &           ,erf_arr_corr,work)
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

*=======================================================================
*---  set fake neighbor lists
*=======================================================================

*=======================================================================
*----- Print titles for the run ----------------------------------------
*=======================================================================


      CALL prtite_min
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

*==== Phony call to forces: Computes only neighbor lists (OLD UPDATE)
c 
      
      CALL timer(vfcp,tfcp,elapse)
      gcpu=tfcp
      IF(lupdate.and.(.not.linked_cell)) THEN
*--      update shell h neighbor list
         CALL mts_forces('u',xpa,ypa,zpa,xpga,ypga,zpga,xpcma,ypcma
     &        ,zpcma,mapnl,mapdn,nmapdn,ucns_p,ucos_p,virs_p,virsp_p
     &        ,ucnp_p,ucop_p,ucnsp_p,ucosp_p,fpx_p,fpy_p,fpz_p
     &        ,stressd_p,worka,cpu_h,ncpu_h
     &        ,nstart_h,nend_h,nstart_ah,nend_ah,nlocal_ah,node
     &        ,nprocs,ncube,P_dyn_update_shell)
      ELSE
         aux= rcuth+rtolh+rneih
         CALL lc_index(indxyz,ncx,ncy,ncz,nind,indxi,indxj,indxk,aux
     &        ,co)
         CALL timer(vfcp,tfcp,elapse)
         gcpu=-gcpu + tfcp
         write(kprint,15011) gcpu
         CALL timer(vfcp,tfcp,elapse)
         gcpu=tfcp
         CALL lc_list(ncx,ncy,ncz,nind,indxi,indxj,indxk,aux,co,xpga
     &        ,ypga,zpga,ngrp,nstart_h,nend_h,node,nprocs,ncube
     &        ,worka,kprint,.TRUE.)
      END IF   
      
      CALL timer(vfcp,tfcp,elapse)
      gcpu=-gcpu + tfcp
      gcpu_u=gcpu
      write(kprint,16011) gcpu
      
*=======================================================================
*---- Zeroes all forces ------------------------------------------------
*=======================================================================
      
      CALL zeroa(fpx,fpy,fpz,ntap,1)
      
      CALL get_total_energy(.TRUE.,mapnl,mapdn,nmapdn,tag_bndg,1,fudgec
     &     ,xp0,yp0,zp0,fpx,fpy,fpz,stressd,stressr,utotal,ucns,ucos
     &     ,urcs,coul_bnd_slv,conf_bnd_slv_n1,coul_bnd_slv_n1,self_slv
     &     ,fsbond,fsbend,fsin14,unb14,cnb14,uslvbon,uslvben,uslvtor
     &     ,uslvitor,uumb,uptors,uitors,ubond,ubend,ucnp,ucop,urcp
     &     ,conf_bnd_slt_n1,coul_bnd_slt,coul_bnd_slt_n1,self_slt,ucnsp
     &     ,ucosp,urcsp,eer)

c$$$
c$$$    TEST
c$$$
      CALL add_energies(pme,pressure,slv_exist,slt_exist,ucns
     &     ,ucos,urcs,coul_bnd_slv,conf_bnd_slv_n1
     &     ,coul_bnd_slv_n1,self_slv,uslvbon,uslvben,uslvtor
     &     ,uslvitor,uumb,uptors,uitors,ubond,ubend,ucnp,ucop
     &     ,urcp,conf_bnd_slt_n1,coul_bnd_slt,coul_bnd_slt_n1
     &     ,self_slt,ucnsp,ucosp,urcsp,eer,uconf,ucoul,ureal
     &     ,urecp,pubnd,purecp,puconf,pucoul,upconf,upcoul
     &     ,prt,stressd,stressr,stress_conf,stress_tot,co,oc,volume
     &     ,unitp,press_conf)
      CALL prtacc(node,pucek,puhyd,puconf,pueng,pucoul,self_slt
     &     ,fsbond,purecp,fsbend,fsin14,unb14,cnb14,ubend,ubond
     &     ,uitors,uptors,pubnd,uceh,hpot,ucoul,uconf,urecp,ureal
     &     ,fsrtal,ucek,upconf,upcoul,uslvbon,uslvben,uslvtor
     &     ,uslvitor,eer,uumb,temp,temph,tcm,rcm,tempt,tempr
     &     ,temppr,gr,gra,ucepr,stress_tot,press_conf,pressc
     &     ,press_kin,temppra,errca,errhe,errbc,erral,drpca,drpbc
     &     ,drphe,drpal,sum_econf,sum_ecoul,sum_enbnd,sum_etotpot
     &     ,sum_eslvint,sum_eebond,sum_eebend,sum_eeptors
     &     ,sum_eeitors,sum_tote,sum_ucek,sum_temp,sum_tempt
     &     ,sum_tempr,sum_temppr,sum_temph,sum_pecek,sum_pehyd
     &     ,sum_peconf,sum_pecoul,sum_percip,sum_enb14,sum_ebend
     &     ,sum_ebond,sum_eitor,sum_eptor,sum_pnbd,sum_pebnd
     &     ,sum_pepot,sum_ptote,sum_gr,sum_epcoul,sum_epconf
     &     ,sum_co,sum_st,sum_presst,sum_press,sum_pressc
     &     ,sum_pressk,sum_volume,sum_pv,sum_temppra,ssm_econf
     &     ,ssm_ecoul,ssm_enbnd,ssm_etotpot,ssm_eslvint,ssm_eebond
     &     ,ssm_eebend,ssm_eeptors,ssm_eeitors,ssm_tote,ssm_ucek
     &     ,ssm_temp,ssm_tempt,ssm_tempr,ssm_temppr,ssm_temph
     &     ,ssm_pecek,ssm_pehyd,ssm_peconf,ssm_pecoul,ssm_percip
     &     ,ssm_enb14,ssm_ebend,ssm_ebond,ssm_eitor,ssm_eptor
     &     ,ssm_pnbd,ssm_pebnd,ssm_pepot,ssm_ptote,ssm_gr
     &     ,ssm_epcoul,ssm_epconf,ssm_co,ssm_st,ssm_presst
     &     ,ssm_press,ssm_pressc,ssm_pressk,ssm_volume,ssm_pv
     &     ,ssm_temppra,energy,nstep,nstep)


c$$$
c$$$   TEST
c$$$
      
      utotal_old=0.0D0
      IF(l_bfgs_b) THEN
         bfgs_nf=3*ntap

*=======================================================================
*---- Allocate memory for BFGS minimizer                            ----
*=======================================================================

         len1=bfgs_nf
         ALLOCATE(rp(len1),fp(len1),nbd(len1))
         len1=3*bfgs_nf
         ALLOCATE(iwa(len1))
         len1=2*bfgs_nf*bfgs_m+4*bfgs_nf+11*bfgs_m*bfgs_m+8*bfgs_m
         ALLOCATE(wa(len1))

         IF(task(1:5) .EQ. 'START') THEN
            CALL zero0(fp,bfgs_nf)
            CALL zeroi(nbd,bfgs_nf)
         END IF
         CALL dcopy(ntap,xp0,1,rp(1),1)
         CALL dcopy(ntap,yp0,1,rp(1+ntap),1)
         CALL dcopy(ntap,zp0,1,rp(1+2*ntap),1)
      END IF
      IF(minimize) THEN
         ok=.FALSE.
         grad_max=-1.0D0
         DO i=1,ntap
            IF(DABS(fpx(i)) .GT. grad_max) grad_max=DABS(fpx(i))
            IF(DABS(fpy(i)) .GT. grad_max) grad_max=DABS(fpy(i))
            IF(DABS(fpz(i)) .GT. grad_max) grad_max=DABS(fpz(i))
         END DO
         
         IF(DABS(grad_max) .LE. eps) THEN
            ok=.TRUE.
            WRITE(kprint,21000) DABS(grad_max)
         ELSE
            WRITE(kprint,21500) DABS(grad_max)
         END IF
         
         DO i=1,ntap
            gpx(i)=-fpx(i)
            gpy(i)=-fpy(i)
            gpz(i)=-fpz(i)
            hpx(i)=gpx(i)
            hpy(i)=gpy(i)
            hpz(i)=gpz(i)
            fpx(i)=hpx(i)
            fpy(i)=hpy(i)
            fpz(i)=hpz(i)
         END DO
         
*=======================================================================
*===== MINIMIZATION LOOP ===============================================
*=======================================================================
         
         CALL timer(vfcp,tfcp,elapse)
         gcpu=tfcp
         elaps=elapse
         nstep=0
         utotal_old=utotal
         DO WHILE(.NOT. ok .AND. nstep .LT. maxstp)
            nstep=nstep+1
            IF(conj_grad .OR. steepest) THEN
               CALL linmin_total(mapnl,mapdn,nmapdn
     &              ,tag_bndg,fudgec,xp0,yp0,zp0,xp1,yp1,zp1,fpx,fpy,fpz
     &              ,fpx1,fpy1,fpz1,ntap,fret)
               task(1:2) = 'FG'
            END IF

200         CONTINUE
            IF(l_bfgs_b) THEN
               CALL setulb(bfgs_nf,bfgs_m,rp,low,up,nbd,utotal
     &           ,fp,BFGS_factr,eps,wa,iwa,task,iprint,csave,lsave,isave
     &           ,dsave)
            END IF
            IF(task(1:2) .EQ. 'FG') THEN 
               IF(l_bfgs_b) THEN
                  CALL dcopy(ntap,rp(1),1,xp0,1)
                  CALL dcopy(ntap,rp(1+ntap),1,yp0,1)
                  CALL dcopy(ntap,rp(1+2*ntap),1,zp0,1)
                  CALL zeroa(fpx,fpy,fpz,ntap,1)
               END IF
               utotal=0.0D0
               CALL get_total_energy(.TRUE.,mapnl,mapdn,nmapdn,tag_bndg
     &              ,1,fudgec,xp0,yp0,zp0,fpx,fpy,fpz,stressd,stressr
     &              ,utotal,ucns,ucos,urcs,coul_bnd_slv,conf_bnd_slv_n1
     &              ,coul_bnd_slv_n1,self_slv,fsbond,fsbend,fsin14,unb14
     &              ,cnb14,uslvbon,uslvben,uslvtor,uslvitor,uumb,uptors
     &              ,uitors,ubond,ubend,ucnp,ucop,urcp,conf_bnd_slt_n1
     &              ,coul_bnd_slt,coul_bnd_slt_n1,self_slt,ucnsp,ucosp
     &              ,urcsp,eer)
            
               CALL add_energies(pme,pressure,slv_exist,slt_exist,ucns
     &              ,ucos,urcs,coul_bnd_slv,conf_bnd_slv_n1
     &              ,coul_bnd_slv_n1,self_slv,uslvbon,uslvben,uslvtor
     &              ,uslvitor,uumb,uptors,uitors,ubond,ubend,ucnp,ucop
     &              ,urcp,conf_bnd_slt_n1,coul_bnd_slt,coul_bnd_slt_n1
     &              ,self_slt,ucnsp,ucosp,urcsp,eer,uconf,ucoul,ureal
     &              ,urecp,pubnd,purecp,puconf,pucoul,upconf,upcoul
     &              ,prt,stressd,stressr,stress_conf,stress_tot,co,oc
     &              ,volume,unitp,press_conf)
               grad_max=-1.0D0
               DO i=1,ntap
                  IF(DABS(fpx(i)) .GT. grad_max) grad_max=DABS(fpx(i))
                  IF(DABS(fpy(i)) .GT. grad_max) grad_max=DABS(fpy(i))
                  IF(DABS(fpz(i)) .GT. grad_max) grad_max=DABS(fpz(i))
               END DO
            
               IF(.NOT. l_bfgs_b .AND. DABS(grad_max) .LE. eps)
     &              ok=.TRUE.
               CALL prtacc(node,pucek,puhyd,puconf,pueng,pucoul,self_slt
     &              ,fsbond,purecp,fsbend,fsin14,unb14,cnb14,ubend,ubond
     &           ,uitors,uptors,pubnd,uceh,hpot,ucoul,uconf,urecp,ureal
     &           ,fsrtal,ucek,upconf,upcoul,uslvbon,uslvben,uslvtor
     &           ,uslvitor,eer,uumb,temp,temph,tcm,rcm,tempt,tempr
     &           ,temppr,gr,gra,ucepr,stress_tot,press_conf,pressc
     &           ,press_kin,temppra,errca,errhe,errbc,erral,drpca,drpbc
     &           ,drphe,drpal,sum_econf,sum_ecoul,sum_enbnd,sum_etotpot
     &           ,sum_eslvint,sum_eebond,sum_eebend,sum_eeptors
     &           ,sum_eeitors,sum_tote,sum_ucek,sum_temp,sum_tempt
     &           ,sum_tempr,sum_temppr,sum_temph,sum_pecek,sum_pehyd
     &           ,sum_peconf,sum_pecoul,sum_percip,sum_enb14,sum_ebend
     &           ,sum_ebond,sum_eitor,sum_eptor,sum_pnbd,sum_pebnd
     &           ,sum_pepot,sum_ptote,sum_gr,sum_epcoul,sum_epconf
     &           ,sum_co,sum_st,sum_presst,sum_press,sum_pressc
     &           ,sum_pressk,sum_volume,sum_pv,sum_temppra,ssm_econf
     &           ,ssm_ecoul,ssm_enbnd,ssm_etotpot,ssm_eslvint,ssm_eebond
     &           ,ssm_eebend,ssm_eeptors,ssm_eeitors,ssm_tote,ssm_ucek
     &           ,ssm_temp,ssm_tempt,ssm_tempr,ssm_temppr,ssm_temph
     &           ,ssm_pecek,ssm_pehyd,ssm_peconf,ssm_pecoul,ssm_percip
     &           ,ssm_enb14,ssm_ebend,ssm_ebond,ssm_eitor,ssm_eptor
     &           ,ssm_pnbd,ssm_pebnd,ssm_pepot,ssm_ptote,ssm_gr
     &           ,ssm_epcoul,ssm_epconf,ssm_co,ssm_st,ssm_presst
     &           ,ssm_press,ssm_pressc,ssm_pressk,ssm_volume,ssm_pv
     &           ,ssm_temppra,energy,nstep,nstep)
               IF(MOD(nstep,nprint) .EQ. 0) THEN
                  WRITE(kprint,20000) DABS(grad_max),utotal-utotal_old
               END IF
               IF(l_bfgs_b) THEN
                  DO i=1,ntap
                     fpx(i)=-fpx(i)
                     fpy(i)=-fpy(i)
                     fpz(i)=-fpz(i)
                  END DO
                  CALL dcopy(ntap,fpx,1,fp(1),1)
                  CALL dcopy(ntap,fpy,1,fp(1+ntap),1)
                  CALL dcopy(ntap,fpz,1,fp(1+2*ntap),1)
               END IF
            ELSE IF(task(1:5) .EQ. 'NEW_X') THEN
               GOTO 200
            ELSE
               ok=.TRUE.
            END IF
            IF(conj_grad .AND. (.NOT. ok)) THEN
               gg=0.0D0
               dgg=0.0D0
               DO i=1,ntap
                  gg=gg+gpx(i)**2+gpy(i)**2+gpz(i)**2
                  dgg=dgg+(fpx(i)+gpx(i))*fpx(i)+(fpy(i)+gpy(i))
     &                 *fpy(i)+(fpz(i)+gpz(i))*fpz(i)
c               dgg=dgg+fpx(i)**2+fpy(i)**2+fpz(i)**2
               END DO
               IF(DABS(gg) .EQ. 0.0D0) GOTO 100
               
               gamma=dgg/gg
               
               DO i=1,ntap
                  gpx(i)=-fpx(i)
                  gpy(i)=-fpy(i)
                  gpz(i)=-fpz(i)
                  hpx(i)=gpx(i)+gamma*hpx(i)
                  hpy(i)=gpy(i)+gamma*hpy(i)
                  hpz(i)=gpz(i)+gamma*hpz(i)
                  fpx(i)=hpx(i)
                  fpy(i)=hpy(i)
                  fpz(i)=hpz(i)
               END DO
            ELSE IF(steepest .AND. (.NOT. ok)) THEN
               DO i=1,ntap
                  fpx(i)=-fpx(i)
                  fpy(i)=-fpy(i)
                  fpz(i)=-fpz(i)
               END DO
            END IF
            
*=======================================================================
*---- Write coordinates to a pdb files when required -------------------
*=======================================================================
            
            IF(nascii .NE. 0) THEN
               IF(MOD(nstep,nascii) .EQ. 0) THEN
                  CALL change_frame(co,oc,-1,ntap,xp0,yp0,zp0,xpa,ypa
     &                 ,zpa)
                  IF(ascii_wsc) THEN
                     CALL tr_wsc(co,xpa,ypa,zpa,xpo,ypo,zpo,mass,nprot
     &                    ,protl)
                  ELSE
                     CALL tr_inbox(xpa,ypa,zpa,xpo,ypo,zpo,mass,nprot
     &                    ,protl)
                  END IF
                  CALL change_frame(co,oc,1,ntap,xpo,ypo,zpo,xpo
     &                    ,ypo,zpo)
                  fstep=time*DBLE(nstep)
                  CALL plotc(co,abmd,gr,gra,fstep,beta,xpo,ypo,zpo,ntap
     &                 ,nres,m1,prsymb,chrge)
               END IF
            END IF
            
*=======================================================================
*---- Dump the restart file when required ------------------------------
*=======================================================================
            
            IF(nsave .NE. 0) THEN
               IF(nstep .NE. 0 .AND. MOD(nstep,nsave).EQ.0)THEN
                  grflag= igr .OR. gofr
                  WRITE(kprint,70200)
                  IF(restart_write) THEN
                     CALL dumprs(kdump_out,restart_out,nstep,temp,ntap
     &                    ,ngrp,nprot,xp0,yp0,zp0,vpx,vpy,vpz,eta,etap,3
     &                    ,thermos,sumarray,ssmarray,navg,navg,co,oc
     &                    ,cpress,vco,grflag,krdf,igrn)
                  END IF
               END IF
            END IF
            IF(lupdate .AND. nupdte .NE. 0) THEN
               IF(.NOT. linked_cell) THEN
*--      update shell h neighbor list
                  IF(MOD(nstep,nupdte) .EQ. 0) CALL mts_forces('u',xpa
     &                 ,ypa,zpa,xpga,ypga,zpga,xpcma,ypcma,zpcma,mapnl
     &                 ,mapdn,nmapdn,ucns_p,ucos_p,virs_p,virsp_p,ucnp_p
     &                 ,ucop_p,ucnsp_p,ucosp_p,fpx_p,fpy_p,fpz_p
     &                 ,stressd_p,worka,cpu_h
     &                 ,ncpu_h,nstart_h,nend_h,nstart_ah,nend_ah
     &                 ,nlocal_ah,node,nprocs,ncube,P_dyn_update_shell)
               ELSE
                  IF(MOD(nstep,nupdte) .EQ. 0) THEN                  
                     aux= rcuth+rtolh+rneih
                     CALL lc_index(indxyz,ncx,ncy,ncz,nind,indxi,indxj
     &                    ,indxk,aux,co)
                     CALL timer(vfcp,tfcp,elapse)
                     gcpu=-gcpu + tfcp
                     write(kprint,15011) gcpu
                     CALL timer(vfcp,tfcp,elapse)
                     gcpu=tfcp
                     CALL lc_list(ncx,ncy,ncz,nind,indxi,indxj,indxk,aux
     &                    ,co,xpga,ypga,zpga,ngrp,nstart_h,nend_h,node
     &                    ,nprocs,ncube,worka,kprint,
     &                    .TRUE.)
                  END IF
               END IF
            END IF
            utotal_old=utotal
         
*=======================================================================
*--- Stop the run smoothly if a file called STOP is found --------------
*=======================================================================

            INQUIRE(FILE='./STOP',EXIST=exist)
            IF(exist) THEN
               IF(node .EQ. 0) THEN
                  grflag= igr .OR. gofr   
                  WRITE(kprint,70200)
                  CALL matinv(3,3,co,oc,volume)
                  volume=volume*boxl**3
                  IF(restart_write) THEN
                     CALL dumprs(kdump_out,restart_out,nstep,temp,ntap
     &                    ,ngrp,nprot,xp0,yp0,zp0,vpx,vpy,vpz,eta,etap,3
     &                    ,thermos,sumarray,ssmarray,navg,navg,co,oc
     &                    ,cpress,vco,grflag,krdf,igrn)
                  END IF

                  IF(wrtgyr) THEN
                     CLOSE(kgyr)
                  END IF
                  WRITE(kprint,70700)
               END IF
               STOP
            END IF
         END DO

100      CONTINUE

      END IF

*=======================================================================
*---  Print Final Results ----------------------------------------------
*=======================================================================
         
      IF(minimize) WRITE(kprint,22000)
      CALL add_energies(pme,pressure,slv_exist,slt_exist,ucns,ucos
     &     ,urcs,coul_bnd_slv,conf_bnd_slv_n1,coul_bnd_slv_n1
     &     ,self_slv,uslvbon,uslvben,uslvtor,uslvitor,uumb,uptors
     &     ,uitors,ubond,ubend,ucnp,ucop,urcp,conf_bnd_slt_n1
     &     ,coul_bnd_slt,coul_bnd_slt_n1,self_slt,ucnsp,ucosp,urcsp
     &     ,eer,uconf,ucoul,ureal,urecp,pubnd,purecp,puconf,pucoul
     &     ,upconf,upcoul,prt,stressd,stressr,stress_conf,stress_tot,co
     &     ,oc,volume,unitp,press_conf)
      
      nprint=1
      CALL prtacc(node,pucek,puhyd,puconf,pueng,pucoul,self_slt,fsbond
     &     ,purecp,fsbend,fsin14,unb14,cnb14,ubend,ubond,uitors,uptors
     &     ,pubnd,uceh,hpot,ucoul,uconf,urecp,ureal,fsrtal,ucek,upconf
     &     ,upcoul,uslvbon,uslvben,uslvtor,uslvitor,eer,uumb,temp,temph
     &     ,tcm,rcm,tempt,tempr,temppr,gr,gra,ucepr,stress_tot
     &     ,press_conf,pressc,press_kin,temppra,errca,errhe,errbc,erral
     &     ,drpca,drpbc,drphe,drpal,sum_econf,sum_ecoul,sum_enbnd
     &     ,sum_etotpot,sum_eslvint,sum_eebond,sum_eebend,sum_eeptors
     &     ,sum_eeitors,sum_tote,sum_ucek,sum_temp,sum_tempt,sum_tempr
     &     ,sum_temppr,sum_temph,sum_pecek,sum_pehyd,sum_peconf
     &     ,sum_pecoul,sum_percip,sum_enb14,sum_ebend,sum_ebond
     &     ,sum_eitor,sum_eptor,sum_pnbd,sum_pebnd,sum_pepot,sum_ptote
     &     ,sum_gr,sum_epcoul,sum_epconf,sum_co,sum_st,sum_presst
     &     ,sum_press,sum_pressc,sum_pressk,sum_volume,sum_pv
     &     ,sum_temppra,ssm_econf,ssm_ecoul,ssm_enbnd,ssm_etotpot
     &     ,ssm_eslvint,ssm_eebond,ssm_eebend,ssm_eeptors,ssm_eeitors
     &     ,ssm_tote,ssm_ucek,ssm_temp,ssm_tempt,ssm_tempr,ssm_temppr
     &     ,ssm_temph,ssm_pecek,ssm_pehyd,ssm_peconf,ssm_pecoul
     &     ,ssm_percip,ssm_enb14,ssm_ebend,ssm_ebond,ssm_eitor,ssm_eptor
     &     ,ssm_pnbd,ssm_pebnd,ssm_pepot,ssm_ptote,ssm_gr,ssm_epcoul
     &     ,ssm_epconf,ssm_co,ssm_st,ssm_presst,ssm_press,ssm_pressc
     &     ,ssm_pressk,ssm_volume,ssm_pv,ssm_temppra,energy,nstep,nstep)
      
*=======================================================================
*---- Write coordinates to a pdb files and dump restart ----------------
*=======================================================================

      IF(minimize) THEN
         IF(nsave .NE. 0) THEN
            grflag= igr .OR. gofr
            WRITE(kprint,70200)
            IF(restart_write) THEN
               CALL dumprs(kdump_out,restart_out,nstep,temp,ntap,ngrp
     &              ,nprot,xp0,yp0,zp0,vpx,vpy,vpz,eta,etap,3,thermos
     &              ,sumarray,ssmarray,navg,navg,co,oc,cpress,vco,grflag
     &              ,krdf,igrn)
            END IF
         END IF
         
         IF(nascii .NE. 0) THEN 
            CALL change_frame(co,oc,-1,ntap,xp0,yp0,zp0,xpa,ypa,zpa) 
            IF(ascii_wsc) THEN 
               CALL tr_wsc(co,xpa,ypa,zpa,xpo,ypo,zpo,mass,nprot,protl) 
            ELSE 
               CALL tr_inbox(xpa,ypa,zpa,xpo,ypo,zpo,mass,nprot 
     &              ,protl) 
            END IF 
            CALL change_frame(co,oc,1,ntap,xpo,ypo,zpo,xpo 
     &           ,ypo,zpo) 
            fstep=time*DBLE(nstep) 
            CALL plotc(co,abmd,gr,gra,fstep,beta,xpo,ypo,zpo,ntap,nres
     &           ,m1,prsymb,chrge) 
         END IF 
      END IF 

      grad_max=-1.0D0
      DO i=1,ntap
         IF(DABS(fpx(i)) .GT. grad_max) grad_max=DABS(fpx(i))
         IF(DABS(fpy(i)) .GT. grad_max) grad_max=DABS(fpy(i))
         IF(DABS(fpz(i)) .GT. grad_max) grad_max=DABS(fpz(i))
      END DO
      WRITE(kprint,20001) DABS(grad_max)
      IF(write_grad) CALL prtfrc(kprint,ngrp,grppt,nres,M1,prsymb
     &     ,beta,xp0,yp0,zp0,fpx,fpy,fpz)
      
*=======================================================================
*---- Compute dynamical matrix -----------------------------------------
*=======================================================================
      
      IF(frequencies) THEN
         WRITE(kprint,23000) 
         mad=m1
         mbd=n_mat
         CALL comp_dynamic_mat(mapnl,mapdn,nmapdn,tag_bndg
     &        ,fudgec,xp0,yp0,zp0,mad,mbd,fpx2,fpy2,fpz2,d_mat,wk,eigvl
     &        ,eigvc)
      END IF
      
*=======================================================================
*===== print potential and equilibirum value  (by matteo) ==============
*=======================================================================

      if(equi.and.slt_exist) THEN
         gcpu = unite*avogad/1000.d0/4.184d0
         write(kequi,*) lbond, ' ===== BOND TABLE  '
         CALL equi_bond(kequi,gcpu,betb,lbnd,lbond,xp0,yp0,zp0,
     &        potbo(1,1),potbo(1,2))
         write(kequi,*) lbend, ' ===== BEND TABLE  '
         CALL equi_bend(kequi,gcpu,betb,lbndg,lbend,xp0,yp0,zp0,
     &        potbe(1,1),potbe(1,2),potbe(1,3),potbe(1,4))
         write(kequi,*) ltors,
     &        ' ===== TORSION TABLE '
         CALL equi_tors(kequi,gcpu,betb,ltor,ltors,xp0,yp0,zp0,
     &        potto(1,1),potto(1,2))
         write(kequi,*) litor,
     &        ' ===== IMPROPER TABLE '
         CALL equi_impro(kequi,gcpu,betb,litr,litor,xp0,yp0,zp0,
     &        potit(1,1),potit(1,2),potit(1,3))
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
      IF(nstep .NE. 0) write(kprint,60200) gcpu/DBLE(nstep)
      IF(nstep .NE. 0) write(kprint,60300) elaps/DBLE(nstep)
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
     &     '      I n t e r m e d i a t e   R e s u l t s        ',
     &     12(' '),'*'/'*',78(' '),'*')
1200  FORMAT(80('*'))
1300  FORMAT('*',78(' '),'*')
13000 FORMAT(/22x,' Temperature has been rescaled ',i5,' times '/)
14000 FORMAT(/22x,'Adjusting bond length to Force Field.'/
     &     22x,'     This will take a while...'/ /) 
70200 FORMAT('<------ Dumping Restart File ------->'/)
70300 FORMAT(/ /'<------ Reading Restart File ------->')
70400 FORMAT('<------ Restart File Read in ------->'/ /)
70100 FORMAT('Velocities have been rescaled ---->'/)

17000 FORMAT(  15x,' Total cpu time for the run     = ',f10.3)
18000 FORMAT(  15x,' Total elapsed time for the run = ',f10.3)
60200 FORMAT(  15x,' Averaged time per step         = ',3x,f7.3)
60300 FORMAT(  15x,' Averaged elapsed per step      = ',3x,f7.3/ /)
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
14004 FORMAT(/5x,'CPUtime for n1-contribution  =',f9.4) 
16004 FORMAT(/5x,'CPUtime for n0-contribution  =',f9.4) 
10067 FORMAT(/5x,'THEORIC SPEED UP FOR NON BONDED PART =',f8.2)
10068 FORMAT(/5x,'OVERALL THEORIC SPEED UP =',f8.2/)
20000 FORMAT(
     &11x,'***********************************************'
     &     ,'**************'/
     &11x,'*      GradMax    ',e12.5,'   Delta_E   ',e12.5,'     *'/
     &11x,'***********************************************'
     &     ,'**************'/)
20001 FORMAT(
     &     21x,'***************************************'/
     &     21x,'*        GradMax    ',e12.5,'      *'/
     &     21x,'***************************************'/)
21000 FORMAT(/ /
     &     21x,'***************************************'/
     &     21x,'*        Already at minimum !!        *'/
     &     21x,'*        GradMax    ',e12.5,'      *'/
     &     21x,'***************************************'/ /)
21500 FORMAT(/ /
     &     21x,'***************************************'/
     &     21x,'*        Not yet at minimum !!        *'/
     &     21x,'*        GradMax    ',e12.5,'      *'/
     &     21x,'***************************************'/ /)
22000 FORMAT(/ /
     &     21x,'***************************************'/
     &     21x,'*                                     *'/
     &     21x,'*        M I N I M I Z A T I O N      *'/
     &     21x,'*                                     *'/
     &     21x,'*       F I N A L   R E S U L T S     *'/
     &     21x,'*                                     *'/
     &     21x,'***************************************'/ /)
23000 FORMAT(/ /
     &     21x,'***************************************'/
     &     21x,'*     Computing Dynamical Matrix      *'/
     &     21x,'*      and Harmonic Frequencies       *'/
     &     21x,'***************************************'/ /)
80033 format(/5x,'Expected CPU time for the RUN:',I4,
     &     ' hours and ',I2, ' min',/,  
     &       /5x,' Expected average time per M step:',f8.2,' sec.'/  
     &       /5x,' Expected average time per femto :',f8.2,' sec.'/)  
60030 FORMAT(/10x,'==========================================='/)
78410 FORMAT
     &     (/ /' *******ERROR: MAXT for PME too small. INCREASE MAXT.'/
     &     /)
70700 FORMAT(/ /15x,' Program Stops Smoothly. Restart Dumped'/ /)
70120 FORMAT('Velocities of the barostat have been rescaled   ---->'/)
70130 FORMAT('Velocities of the thermostats have been rescaled   ---->'/
     &     )
20900 FORMAT(/ /' ****** WARNING ! WARNING ! WARNING ***************' /
     &        /' ******   drift_remove ON           ***************' /
     &     /' ****** WARNING ! WARNING ! WARNING ***************' / /)
10072 FORMAT(/ /' --- Total energy remove per particle = ', f12.4,
     &         ' (KJ/mole)'/
     &         ' --- Number of dirty scaling          = ', I10/
     &         ' --- Frequency of scaling             = ',f10.3,
     &         ' 1/ps '/)
10977 FORMAT(/ /'*******WARNING: NO COFACTOR ATOMS SELECTED '/ 
     &     ' NCOFACTOR IS SET TO ZERO AND NO FIELD I COMPUTED'/ /)
80000 FORMAT('REMARK   Rigid body fit on CA atoms')
80100 FORMAT('REMARK   Rigid body fit on heavy atoms')
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
      RETURN
      END
