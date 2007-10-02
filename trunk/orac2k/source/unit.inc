      INTEGER nfile_dump
      PARAMETER(nfile_dump=40)
      INTEGER kread,kprint,kdump_in,kdump_out,kgdata,kcnfig_i,kcnfig_o
     &     ,knlist,kgofr,kplot,kconf,ksofk,kpot,kout,kfield,kvaf,krms
     &     ,kxrms,kxrms_atm,ksol1,ksol2,kgroup,ktop,ktemplate,ktest
     &     ,iurest,kspec,kgr,kplot_fragm,kfragm_dist,kplot_center
     &     ,kgofr_sk,kvi,kcoord_slv,kcoord_slt,ktpgprm_write
     &     ,ktpgprm_read,kavg,kavg_xrms,kwrite_dump_i(nfile_dump)
     &     ,kwrite_dump_o(nfile_dump),no_records_i(nfile_dump)
     &     ,no_records_o(nfile_dump),nwrite_dump_i,nwrite_dump_o
     &     ,recl1_dump_o,recl1_dump_i,recl2_dump_o,recl2_dump_i
     &     ,kvoronoi,kcavities,ktopol,khbonds,khhisto,kabmd,kdipole
     &     ,knative,knative_tpl,kdiff,kfreq,kequi,kout_ef,kout_dp
     &     ,kout_dp_sta,kout_dp_tot,kpol_inp,kxslt,kfit,kcalc_cofm
     &     ,klda_inst,klda_eend,klda_rmin,klda_flu,klda_hyd,kprot_hyd
     &     ,kprot_rest,kprot_lda,klda_rest,KElecPot,khydration
      CHARACTER*80 file_names_i(nfile_dump),file_names_o(nfile_dump)
     &     ,FileElecPot
      REAL*8 pi,efact,avogad,boltz,gascon,planck,elechg,epso,twopi,boxl
      REAL*8 unitm,unitl,unitt,uniti,unite,unitc,unitp,ttstep
      REAL*8 unitefield,unitepot,lbohr
      COMMON /cmcons/ pi,efact,avogad,boltz,gascon,planck
     &  ,elechg,epso,twopi,boxl
      COMMON /cmconv/ unitm,unitl,unitt,uniti,unite,unitc,unitp,ttstep
     &               ,unitefield,unitepot,lbohr
      COMMON /cmunit/ kread,kprint,kdump_in,kdump_out,kgdata,kcnfig_i
     &     ,kcnfig_o,knlist,kgofr,kplot,kconf,kpot,kout,kfield,kvaf,krms
     &     ,kxrms,kxrms_atm,ksol1,ksol2,kgroup,ktop,ktemplate,ktest
     &     ,iurest,kspec,kgr,kplot_fragm,kfragm_dist,kplot_center
     &     ,kgofr_sk,ksofk,kvi,kcoord_slv,kcoord_slt,ktpgprm_write
     &     ,ktpgprm_read,kavg,kavg_xrms,kwrite_dump_i,kwrite_dump_o
     &     ,kvoronoi,kcavities,ktopol,khbonds,khhisto,kabmd,kdipole
     &     ,knative,knative_tpl,kdiff,nwrite_dump_i,nwrite_dump_o
     &     ,no_records_i,no_records_o,recl1_dump_o,recl1_dump_i
     &     ,recl2_dump_o,recl2_dump_i,kfreq,kequi,file_names_i
     &     ,file_names_o,FileElecPot,kout_ef,kout_dp,kout_dp_sta
     &     ,kout_dp_tot,kpol_inp,kxslt,kfit,kcalc_cofm,klda_inst
     &     ,klda_eend,klda_rmin,klda_flu,klda_hyd,kprot_hyd,kprot_rest
     &     ,kprot_lda,klda_rest,KElecPot,khydration


