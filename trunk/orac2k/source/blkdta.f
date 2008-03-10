      SUBROUTINE blkdta
************************************************************************
*                                                                      *
*     Subroutine BLKDTA set up defaults for the simulation             *
*     parameters.                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************

      IMPLICIT none

      INTEGER i,j

      include 'parst.h'
      include 'cpropar.h'
      include 'parallel.h'

      INCLUDE 'unit.h'
      COMMON /card / jrec,nrec,jerr,nerr,nline,istrt
      INTEGER jrec,nrec,jerr,nerr,nline,istrt(80)

*=======================================================================
*------- Set defaults parameters ---------------------------------------
*=======================================================================

      knlist=5
      kprint=6
      kdump_in=0
      kdump_out=0
      kgdata=0
      kcnfig_i=0
      kcnfig_o=0
      kgofr=0
      kplot=0
      kconf=0
      kpot=0
      kout=0
      kfield=0
      kvaf=0
      krms=0
      kxrms=0
      kxrms_atm=0
      ksol1=0
      ksol2=0
      kgroup=0
      ktop=0
      ktemplate=0
      ktest=0
      iurest=0
      kspec=0
      kgr=0
      kplot_fragm=0
      kplot_center=0
      kgofr_sk=0
      kcoord_slt=0
      kcoord_slv=0
      ktpgprm_read=0
      ktpgprm_write=0
      kfreq=0
      kout_ef=0
      kout_dp=0
c-----
      klda_eend=0
      klda_inst=0
      klda_rmin=0
      kcalc_cofm=0
      klda_flu=0
      klda_hyd=0
      kprot_hyd=0
      kprot_rest=0
      kpol_inp=0
      kout_dp=0
      kout_dp_sta=0
      kout_dp_tot=0
c-----      
      annpro=0
      jerr=1
      P_nei_fact=2.0D0
      P_dyn_update=.FALSE.
      P_dyn_update_shell='l'
      frequencies=.FALSE.
      debug=.FALSE.
      debug_rs=.FALSE.
      debug_ct=.FALSE.
      debug_st=.FALSE.
      debug_bt=.FALSE.
      debug_pt=.FALSE.
      debug_it=.FALSE.
      steepest=.FALSE.
      conj_grad=.FALSE.
      l_bfgs_b=.FALSE.
      energy_then_die=.FALSE.
      cutoff=.TRUE.
      skip_step=.FALSE.
      scan_traj=.FALSE.
      grpcut=.FALSE.
      protei=.FALSE.
      solven=.FALSE.
      electr=.FALSE.
      clewld=.FALSE.
      linser=.FALSE.
      AddTime=.FALSE.
      recstrc=.FALSE.
      stinit = .FALSE.
      tpgfil = .FALSE.
      rescm  = .FALSE.
      tpgwbn = .FALSE.
      prttpg = .FALSE.
      sction=.FALSE.
      prtatl=.FALSE.
      prthbl=.FALSE.
      prtnbl=.FALSE.
      prtn14=.FALSE.
      prtcnl=.FALSE.
      prtbal=.FALSE.
      prtptl=.FALSE.
      prtitl=.FALSE.
      prtipw=.FALSE.
      hmass=.FALSE.
      prttor=.FALSE.
      prtben=.FALSE.
      prtstr=.FALSE.
      igr=.FALSE.
      gofr=.FALSE.
      restart_old=.FALSE.
      restart_read=.FALSE.
      restart_write=.FALSE.
      old_tpg=.FALSE.
      read_co=.FALSE.
      gofr_nprint=0
      gofr_navg=0
      gofr_ncomp=0
      gofr_neighbor=.FALSE.
      gofr_intra=.FALSE.
      gofr_avg=.FALSE.
      gofr_cut=0.0D0
      delrg=0.05D0
      lenerg=.FALSE.
      uenerg=.FALSE.
      fold=.FALSE.
      wrtgyr=.FALSE.
      noneq=.FALSE.
      lfield=.FALSE.
      scharge=.FALSE.
      pfix=.FALSE.
      nonbnd=.TRUE.
      hoover=.FALSE.
      thermos=.FALSE.
      vacf=.FALSE.
      mdsim=.TRUE.
      minimize=.FALSE.
      write_grad=.FALSE.
      annealing=.FALSE.
      change_cell=.FALSE.
      resize_cell=.FALSE.
      analys=.FALSE.
      anxrms=.FALSE.
      anxrms_cell=.FALSE.
      ascii_nocell=.FALSE.
      anxca=.FALSE.
      anxbc=.FALSE.
      anxhe=.FALSE.
      anrms=.FALSE.
      anslv=.FALSE.
      stretch=.FALSE.
      stretch_heavy=.FALSE.
      stoprun=.FALSE.
      md_respa=.FALSE.
      adihed=.TRUE.
      hydbnd=.FALSE.
      spectra=.FALSE.
      lgr=.FALSE.
      linked_cell=.FALSE.
      adjust_cnstr=.TRUE.
      readjust_cnstr=.FALSE.
      slv_read=.FALSE.
      slv_exist=.FALSE.
      slt_exist=.FALSE.
      slv_generate=.FALSE.
      slv_create=.FALSE.
      pdb_read=.FALSE.
      slt_create=.FALSE.
      abmd=.FALSE.
      dissociate=.FALSE.
      associate=.FALSE.
      diss_mol=.FALSE.
      diss_atoms=.FALSE.
      abmd_tors=.FALSE.
      occupy_space=.FALSE.
      diffusion=.FALSE.
      time_corr=.FALSE.
      not_time_corr=.FALSE.
      erf_corr=.FALSE.
      erfc_spline=.FALSE.
      remove_momentum=.FALSE.
      anprot=.FALSE.
      Virtual_residue=.FALSE.
      virtual_atoms(1)=0
      divide_records=1
      divide_spline=1
      ns1=256
      nv1=256
      ns2=512
      nv2=512
      ns3=1024
      nv3=1024
      ns4=2048
      nv4=2048
      niop=1
      nios=1
      rpass =0.005
      ngprint = 10000
      nprot_charges=0
      corr_atoms(1)=0
      do i=1,5
         ntmtsp(i)=-1
      end do
      do i=1,3
         ntmtss(i)=-1
      end do
      bfgs_m=4
      bfgs_factr=1.0D7
      totvaf=300
      novaf=3
      pext=0.0D0
      taut=0.0D0
      taup=0.0D0
      volumepr=0.0D0
      compressibility=5.3D-04
      isostress=.FALSE.
      cpress=.FALSE.
      pressure=.FALSE.
      lberendsen=.FALSE.
      landersen=.FALSE.
      IF(nprocs .EQ. 1) THEN
         coupl_mol=.TRUE.
         coupl_grp=.FALSE.
         coupl_atm=.FALSE.
      ELSE
         coupl_mol=.FALSE.
         coupl_grp=.TRUE.
         coupl_atm=.FALSE.
      END IF

      native=.FALSE.
      check_native=.FALSE.
      wpr=40.0D0
      ecut=1.0D0
      epcut=1.0D0
      echrg=-11.0D0
      pec=8
      nec=100
      rspcut=1.0D0
      rspoff=0.0D0
      rspon=1.0D0
      hrcut=1.0D0
      hacut=10.0D0
      hrson=3.50D0
      hrsoff=4.50D0
      hanon=50.0D0
      hanoff=70.0D0
      rcut_hb=4.5D0
      acut_hb=0.0D0
      a2cut_hb=0.0D0
      TimeLimit=1.0D+30
      nhskip=3
      nupdte=20
      nupdte_index=1
      alphal=0.4D0
      omega=1.0D0
      spring=0.0D0
      lmbdt=0.5D0
      rgnato=0
      fudge=1.0D0
      lj_fudge=1.0D0
      ancorr=100
      antstar=0.0D0
      anfact=1.0D0
      dtemph=5000.0D0
      dtemppr=5000.0D0
      hstep_freq=0.03D0
      nstep_freq=6
      native_dist=3.2D0
      native_theta=10.0D0
      neta=3
      DO i=1,neta
         qmass(i)=0.0D0
      END DO
      scale=1E8
      iseed=12345667
      compressibility=4.9D-04
      eps_energy=0.5D0
      itor_ptype=0
      nospt=3
      nplot_fragm=0
      nplot_center=0
      ncofactor=0
      nfragm=0
      nrespa=1
      n0respa=1
      n1respa=1
      mrespa=1
      lrespa=1
      ltest_times=.false.
      clean=.true.
      start_conf= .false.
      distmax = 0.25
      rcuth=9.7
      rtolh=0.3
      rcutl=7.3
      rtoll=0.3
      rcutm=4.1
      rtolm=0.3
      rneih=1.5 
      rneil=0.45
      rneim=0.35
      kl=0.0
      klt=0.0
      km=0.0
      kmt=0.0
      erfc_bin=0.007D0
      hhisto_bin=0.007D0
      pme=.false.
      shell_pme='0'
      nfft1=0
      nfft2=0
      nfft3=0
      mim_lim=4
      pme_order=0
      DO 10 i=1,8
          DO 20 j=1,3
              rmol(j,i)=0.0D0
20        CONTINUE
10    CONTINUE
      mixrww=.true.
      mixrwp=.false.
      amphi=.FALSE.
      bending=.TRUE.
      replicate=.FALSE.
      write_ff_pars=.FALSE.
      nplot=0
      nascii=0
      nfield=0
      mfield=1
      nvi = 0
      ndipole = 0
      nhoov=1
      whoov=500.0D0
      time=1.0D0
      radius=1.6d+0
      rkcut=1.0D9
      nflag(1)=0
      nflag(2)=0
      maxstp=1
      maxrun=0
      nrject=1
      nprop=200
      nsave=0
      nconf=0
      nprint=1
      menerg=0
      t=300.0D0
      baxis=0.0D0
      caxis=0.0D0
      nform=1
      icl=1
      icm=1
      icn=1
      nmol=0
      nato_slv=0
      start_anl=0
      stop_anl=2**(24)
      update_anl=0
      nhbonds=0
      nhhisto=0
      nx_cav=0
      ny_cav=0
      nz_cav=0
      size_atom_cav=2.5D0
      bin_size_cav=0.01D0
      rmax_size_cav=3.0D0
      ncavities=0
      ChargeIon=0.0D0
      SigmaIon=0.0D0
      mass_pfix=0.0D0
      SmoothFactor=0.0D0
      DoFreeEnergy=.FALSE.
      cavities=.FALSE.
      hbonds_res=.FALSE.
      hbonds_tot=.FALSE.
      hbonds_vor=.FALSE.
      prttopl=.FALSE.
      fragm_dist=.FALSE.
c---
      efield=.FALSE.
      lmol_ef=.FALSE.
      print_ef=.FALSE.
      nmol_ef=0
      nfreq_ef=0
      print_dp=.FALSE.
      nfreq_dp=0
     
      latms=.FALSE.
      natms=0
      sel_lda=.FALSE.
      nlda_mol=0
      nlda_atm=0
      nlda_zero=0
      avg_lda=.FALSE.
      ninst_lda=0
      calc_cofm=.FALSE.
      ninst_fit=0
      inst_fit=.FALSE.
c---
      lda_flu=.FALSE.
      lda_hyd=.FALSE.
      nlda_flu=0
      nlda_hyd=0
      prot_hyd=.FALSE.
      res_time=.FALSE.
      SecStructure=.FALSE.
      nprot_hyd=0
      coeff_hyd = 1.3
      min_hyd=1
      max_hyd=0
      polar=.FALSE.
      NFreqPotential=1
      polar_scale = 1.0d0
      nfreq_polar = 1
      Ext_ef(1) = 0.0d0
      Ext_ef(2) = 0.0d0
      Ext_ef(3) = 0.0d0

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END

