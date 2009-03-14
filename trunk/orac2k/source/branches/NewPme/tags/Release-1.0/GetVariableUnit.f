      SUBROUTINE GetVariableUnit(Pkread,Pkprint,Pkdump_in,Pkdump_out
     &     ,Pkgdata,Pkcnfig_i,Pkcnfig_o,Pknlist,Pkgofr,Pkplot,Pkconf
     &     ,Pksofk,Pkpot,Pkout,Pkfield,Pkvaf,Pkrms,Pkxrms,Pkxrms_atm
     &     ,Pksol1,Pksol2,Pkgroup,Pktop,Pktemplate,Pktest,Piurest,Pkspec
     &     ,Pkgr,Pkplot_fragm,Pkfragm_dist,Pkplot_center,Pkgofr_sk,Pkvi
     &     ,Pkcoord_slv,Pkcoord_slt,Pktpgprm_write,Pktpgprm_read,Pkavg
     &     ,Pkavg_xrms,Pnwrite_dump_i,Pnwrite_dump_o,Pkvoronoi
     &     ,Pkcavities,Pktopol,Pkhbonds,Pkhhisto,Pkabmd,Pkdipole
     &     ,Pknative,Pknative_tpl,Pkdiff,Pkfreq,Pkequi,Pkout_ef,Pkout_dp
     &     ,Pkout_dp_sta,Pkout_dp_tot,Pkpol_inp,Pkxslt,Pkfit,Pkcalc_cofm
     &     ,Pklda_inst,Pklda_eend,Pklda_rmin,Pklda_flu,Pklda_hyd
     &     ,Pkprot_hyd,Pkprot_rest,Pkprot_lda,Pklda_rest,Ppi,Pefact
     &     ,Pavogad,Pboltz,Pgascon,Pplanck,Pelechg,Pepso,Ptwopi,Pboxl
     &     ,Punitm,Punitl,Punitt,Puniti,Punite,Punitc,Punitp,Pttstep
     &     ,Punitefield,Punitepot,Plbohr)

!!$***********************************************************************
!!$   Time-stamp: <01/04/02 15:24:09 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Mon Apr  2 2001 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*

      IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*


!!$----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'unit,h'

!!$------------------------- LOCAL VARIABLES ----------------------------*


!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

      IF(PRESENT(Pkread)) THEN
         Pkread=kread
      END IF
      IF(PRESENT(Pkprint)) THEN
         Pkprint=kprint
      END IF
      IF(PRESENT(Pkdump_in)) THEN
         Pkdump_in=kdump_in
      END IF
      IF(PRESENT(Pkdump_out)) THEN
         Pkdump_out=kdump_out
      END IF
      IF(PRESENT(Pkgdata)) THEN
         Pkgdata=kgdata
      END IF
      IF(PRESENT(Pkcnfig_i)) THEN
         Pkcnfig_i=kcnfig_i
      END IF
      IF(PRESENT(Pkcnfig_o)) THEN
         Pkcnfig_o=kcnfig_o
      END IF
      IF(PRESENT(Pknlist)) THEN
         Pknlist=knlist
      END IF
      IF(PRESENT(Pkgofr)) THEN
         Pkgofr=kgofr
      END IF
      IF(PRESENT(Pkplot)) THEN
         Pkplot=kplot
      END IF
      IF(PRESENT(Pkconf)) THEN
         Pkconf=kconf
      END IF
      IF(PRESENT(Pksofk)) THEN
         Pksofk=ksofk
      END IF
      IF(PRESENT(Pkpot)) THEN
         Pkpot=kpot
      END IF
      IF(PRESENT(Pkout)) THEN
         Pkout=kout
      END IF
      IF(PRESENT(Pkfield)) THEN
         Pkfield=kfield
      END IF
      IF(PRESENT(Pkvaf)) THEN
         Pkvaf=kvaf
      END IF
      IF(PRESENT(Pkrms)) THEN
         Pkrms=krms
      END IF
      IF(PRESENT(Pkxrms)) THEN
         Pkxrms=kxrms
      END IF
      IF(PRESENT(Pkxrms_atm)) THEN
         Pkxrms_atm=kxrms_atm
      END IF
      IF(PRESENT(Pksol1)) THEN
         Pksol1=ksol1
      END IF
      IF(PRESENT(Pksol2)) THEN
         Pksol2=ksol2
      END IF
      IF(PRESENT(Pkgroup)) THEN
         Pkgroup=kgroup
      END IF
      IF(PRESENT(Pktop)) THEN
         Pktop=ktop
      END IF
      IF(PRESENT(Pktemplate)) THEN
         Pktemplate=ktemplate
      END IF
      IF(PRESENT(Pktest)) THEN
         Pktest=ktest
      END IF
      IF(PRESENT(Piurest)) THEN
         Piurest=iurest
      END IF
      IF(PRESENT(Pkspec)) THEN
         Pkspec=kspec
      END IF
      IF(PRESENT(Pkgr)) THEN
         Pkgr=kgr
      END IF
      IF(PRESENT(Pkplot_fragm)) THEN
         Pkplot_fragm=kplot_fragm
      END IF
      IF(PRESENT(Pkfragm_dist)) THEN
         Pkfragm_dist=kfragm_dist
      END IF
      IF(PRESENT(Pkplot_center)) THEN
         Pkplot_center=kplot_center
      END IF
      IF(PRESENT(Pkgofr_sk)) THEN
         Pkgofr_sk=kgofr_sk
      END IF
      IF(PRESENT(Pkvi)) THEN
         Pkvi=kvi
      END IF
      IF(PRESENT(Pkcoord_slv)) THEN
         Pkcoord_slv=kcoord_slv
      END IF
      IF(PRESENT(Pkcoord_slt)) THEN
         Pkcoord_slt=kcoord_slt
      END IF
      IF(PRESENT(Pktpgprm_write)) THEN
         Pktpgprm_write=ktpgprm_write
      END IF
      IF(PRESENT(Pktpgprm_read)) THEN
         Pktpgprm_read=ktpgprm_read
      END IF
      IF(PRESENT(Pkavg)) THEN
         Pkavg=kavg
      END IF
      IF(PRESENT(Pkavg_xrms)) THEN
         Pkavg_xrms=kavg_xrms
      END IF
      IF(PRESENT(Pkwrite_dump_i)) THEN
         Pkwrite_dump_i=kwrite_dump_i
      END IF
      IF(PRESENT(Pkwrite_dump_o)) THEN
         Pkwrite_dump_o=kwrite_dump_o
      END IF
      IF(PRESENT(Pno_records_i)) THEN
         Pno_records_i=no_records_i
      END IF
      IF(PRESENT(Pno_records_o)) THEN
         Pno_records_o=no_records_o
      END IF
      IF(PRESENT(Pnwrite_dump_i)) THEN
         Pnwrite_dump_i=nwrite_dump_i
      END IF
      IF(PRESENT(Pnwrite_dump_o)) THEN
         Pnwrite_dump_o=nwrite_dump_o
      END IF
      IF(PRESENT(Pkvoronoi)) THEN
         Pkvoronoi=kvoronoi
      END IF
      IF(PRESENT(Pkcavities)) THEN
         Pkcavities=kcavities
      END IF
      IF(PRESENT(Pktopol)) THEN
         Pktopol=ktopol
      END IF
      IF(PRESENT(Pkhbonds)) THEN
         Pkhbonds=khbonds
      END IF
      IF(PRESENT(Pkhhisto)) THEN
         Pkhhisto=khhisto
      END IF
      IF(PRESENT(Pkabmd)) THEN
         Pkabmd=kabmd
      END IF
      IF(PRESENT(Pkdipole)) THEN
         Pkdipole=kdipole
      END IF
      IF(PRESENT(Pknative)) THEN
         Pknative=knative
      END IF
      IF(PRESENT(Pknative_tpl)) THEN
         Pknative_tpl=knative_tpl
      END IF
      IF(PRESENT(Pkdiff)) THEN
         Pkdiff=kdiff
      END IF
      IF(PRESENT(Pkfreq)) THEN
         Pkfreq=kfreq
      END IF
      IF(PRESENT(Pkequi)) THEN
         Pkequi=kequi
      END IF
      IF(PRESENT(Pkout_ef)) THEN
         Pkout_ef=kout_ef
      END IF
      IF(PRESENT(Pkout_dp)) THEN
         Pkout_dp=kout_dp
      END IF
      IF(PRESENT(Pkout_dp_sta)) THEN
         Pkout_dp_sta=kout_dp_sta
      END IF
      IF(PRESENT(Pkout_dp_tot)) THEN
         Pkout_dp_tot=kout_dp_tot
      END IF
      IF(PRESENT(Pkpol_inp)) THEN
         Pkpol_inp=kpol_inp
      END IF
      IF(PRESENT(Pkxslt)) THEN
         Pkxslt=kxslt
      END IF
      IF(PRESENT(Pkfit)) THEN
         Pkfit=kfit
      END IF
      IF(PRESENT(Pkcalc_cofm)) THEN
         Pkcalc_cofm=kcalc_cofm
      END IF
      IF(PRESENT(Pklda_inst)) THEN
         Pklda_inst=klda_inst
      END IF
      IF(PRESENT(Pklda_eend)) THEN
         Pklda_eend=klda_eend
      END IF
      IF(PRESENT(Pklda_rmin)) THEN
         Pklda_rmin=klda_rmin
      END IF
      IF(PRESENT(Pklda_flu)) THEN
         Pklda_flu=klda_flu
      END IF
      IF(PRESENT(Pklda_hyd)) THEN
         Pklda_hyd=klda_hyd
      END IF
      IF(PRESENT(Pkprot_hyd)) THEN
         Pkprot_hyd=kprot_hyd
      END IF
      IF(PRESENT(Pkprot_rest)) THEN
         Pkprot_rest=kprot_rest
      END IF
      IF(PRESENT(Pkprot_lda)) THEN
         Pkprot_lda=kprot_lda
      END IF
      IF(PRESENT(Pklda_rest)) THEN
         Pklda_rest=klda_rest
      END IF
      IF(PRESENT(Pfile_names_i)) THEN
         Pfile_names_i=file_names_i
      END IF
      IF(PRESENT(Ppi)) THEN
         Ppi=pi
      END IF
      IF(PRESENT(Pefact)) THEN
         Pefact=efact
      END IF
      IF(PRESENT(Pavogad)) THEN
         Pavogad=avogad
      END IF
      IF(PRESENT(Pboltz)) THEN
         Pboltz=boltz
      END IF
      IF(PRESENT(Pgascon)) THEN
         Pgascon=gascon
      END IF
      IF(PRESENT(Pplanck)) THEN
         Pplanck=planck
      END IF
      IF(PRESENT(Pelechg)) THEN
         Pelechg=elechg
      END IF
      IF(PRESENT(Pepso)) THEN
         Pepso=epso
      END IF
      IF(PRESENT(Ptwopi)) THEN
         Ptwopi=twopi
      END IF
      IF(PRESENT(Pboxl)) THEN
         Pboxl=boxl
      END IF
      IF(PRESENT(Punitm)) THEN
         Punitm=unitm
      END IF
      IF(PRESENT(Punitl)) THEN
         Punitl=unitl
      END IF
      IF(PRESENT(Punitt)) THEN
         Punitt=unitt
      END IF
      IF(PRESENT(Puniti)) THEN
         Puniti=uniti
      END IF
      IF(PRESENT(Punite)) THEN
         Punite=unite
      END IF
      IF(PRESENT(Punitc)) THEN
         Punitc=unitc
      END IF
      IF(PRESENT(Punitp)) THEN
         Punitp=unitp
      END IF
      IF(PRESENT(Pttstep)) THEN
         Pttstep=ttstep
      END IF
      IF(PRESENT(Punitefield)) THEN
         Punitefield=unitefield
      END IF
      IF(PRESENT(Punitepot)) THEN
         Punitepot=unitepot
      END IF
      IF(PRESENT(Plbohr)) THEN
         Plbohr=lbohr
      END IF



!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      END SUBROUTINE GetVariableUnit
