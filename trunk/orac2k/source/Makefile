# Makefile for ORAC. Works with ULTRIX, OSF, HPUX, AIX, 
# IRIX and IRIX64 machines
# 
# 

SHELL = /bin/sh
DESCRIPTION = Merge version MAX-PIE. Unit of time is fs.
STATUS = beta version 
RELEASE = 1.0
CPP = /lib/cpp -traditional
CC = cc
LIB = nolib
CURDIR = `pwd`
DIST_DIR = tmp_dir/orac/source

#==========================================================================
# Define default rules
#==========================================================================

.f.o :
	$(FC) $(FFLAGS) -D$(TARGET0) -c $<
.SUFFIXES: .o .f90
.f90.o :
	$(F90) $(FFLAGS) -D$(TARGET0) -c $<

#==========================================================================
# Define include files
#==========================================================================


INCLUDE = cpropar.inc parst.inc orac.inc fourier.inc unit.inc \
	resid.inc vacf.inc pbc.inc analysis.inc mts_spectra.inc\
	mts_vacf.inc pme.inc parallel.inc parameters.inc\
	lc_list.inc voronoi.inc\
	iobuffer.inc pressure.h DXMLDEF.FOR mtsmd_avg_inc.f \
	mtsmd_dump_inc.f

INCPRE = $(patsubst %.inc, %.h, $(INCLUDE)) 

#==========================================================================
# Object modules
#==========================================================================

OBJ_NP = xerror_mod.o Module_Neighbors.o RMS_Matrix_Mod.o\
	EUL_Mod.o ElecPotential_Mod.o Gauss_Overlap_Mod.o \
	Neighbors_Mod.o DENSITY_Mod.o Class_Connect.o Class_Bending.o \
	Class_Gauss_Param.o Module_Fourier.o Class_Gcharges.o\
	Module_Extra_Forces.o reset_nbond.o Poisson_Boltzmann.o\
	Module_Stress.o Module_Thole.o CENTER_SOL_Mod.o\
	INPUT_Mod.o PDBs_Mod.o RMS_Subtract_Mod.o\
	add_energies.o add_int1.o add_int2.o add_real1.o add_real2.o\
	add_char1.o\
	add_concta.o add_mapnl.o add_solvent_tpg.o add_solvent_coord.o\
	add1h.o add2h.o add3h.o add_str.o addh.o addh2o.o\
	addlnk.o addtor.o addtpg.o adjust_bonds.o\
	anneal.o appatm.o appbkn.o appbou.o appct.o appgrp.o\
	asng_xrmsw.o asstpe.o atoms.o atthd.o\
	bcnstp.o bend.o blkdta.o bond.o brent.o brot.o CalcDistConf.o\
	calc_avg_str.o calc_avg2_str.o calc_avg_xrms.o calc_gofr.o\
	calc_sk.o calc_sk_eet.o calc_xslt.o calc_inst_fit.o \
	calc_RotMat.o calc_cofm.o\
	calc_lda.o calc_lda_avg.o calc_lda_flu.o calc_lda_hyd.o calc_hyd.o\
	calc_xrms.o CalcXrmsSecStruct.o caldr.o cc.o change_coord_ex.o\
	change_coord_inner.o change_virtual_potential.o virtual_extract.o\
	calc_lda_rest.o change_frame.o change_origin.o\
	change_tors.o change_step_sim.o\
	check_length_sim.o check_read_columns.o check_topology.o\
	check_zero_coord.o chk14.o chkewl.o\
	chlab.o chmass.o chnpr.o close_files.o clpcm.o clt14.o\
	collision.o comd0.o \
	comp_abmd_fdiss.o comp_abmd_fcryst.o\
	comp_ang_vel.o comp_cavities_vor.o\
	comp_abmd_tors.o comp_cell_center.o comp_concth.o\
	comp_contacts.o comp_dip.o comp_dynamic_mat.o\
	comp_fcm.o comp_forcep.o comp_gmgp.o comp_molmass.o\
	comp_native.o comp_neigh_vor.o comp_noneq_free.o comp_nmol.o\
	comp_stress_conf.o\
	comp_stress_kinetic.o comp_stress_bnd.o comp_thermos_energy.o\
	comp_rmsq.o comp_thermos_forces.o\
	comp_thermos_forces_cm.o comp_vcm.o comp_vel_labframe.o\
	correc.o correc_etap.o Correct_Ewald.o\
	coordinate_spline.o coordinate_spline_init.o copy_protl.o\
	correc3x3.o correc_gmgp.o correc_matr.o correc_exp_scale.o\
	correc_scale.o correc_scale_2.o correc_scale_co.o cov_thermos.o\
	covbd.o covchg.o covfod.o covnb.o covpr.o covuen.o crdrtf.o cself.o\
	cself_dipole.o \
	csort.o dblas.o dcopy3.o delrfg.o dist_ijk.o dlineq.o\
	DCopy_V3.o \
	dspat.o dumptp.o duni.o dupl_char1.o dupl_concta.o\
	dupl_int1.o dupl_int2.o dupl_mapnl.o dupl_real1.o dupl_real2.o\
	dupl_reali2.o erf_corr_cutoff.o erfc_spline.o erfc_spline_init.o\
	equi_bond.o equi_bend.o equi_tors.o equi_impro.o\
	f1dim_adjust_bonds.o \
	ferrf.o fit.o ferrf_tag.o fftpack.o fft3-public.o fft_default.o\
	fft_pme_init.o fill_charge_grid.o fill_dipole_grid.o\
	find_diss_mol.o find_igint_vor.o find_length_run.o\
	find_thermos.o fitors.o fix_molecules.o \
	fnb14.o fnb14_2.o fnbgrp.o fnbgrp2.o fndfmt.o fndgrp.o fpbend.o fpbond.o\
	fpbond_adjust_bonds.o fptors.o gauss.o gen_abmd_kvect.o\
	gen_stress_tag_bnd.o generate_slv.o generate_rnd_rot.o\
	get_atres.o get_bspline_coeffs.o get_bspline_coeffs2.o GetKineticTest.o\
	GetHardSphere.o get_virtual_energy.o\
	Get_vp_with_vc.o Get_vp_without_vc.o \
	get_pointer_protl.o get_native.o get_prot_cnstr.o get_sgr.o\
	get_spectra_vacf.o get_type_slv.o iceil.o\
	get_velocities.o getchr.o GetSimpleCO.o grad_sum.o grad_sum2.o\
	grad_sum_dipole.o hapacp.o hapdon.o igmap.o\
	indexx.o iindexx.o inicmp.o initialize_cavities.o\
	Init_Fixed_Molecule.o\
	int_corr_erf_spline.o interpol_dyna.o chebft_dyn.o\
	chder_dyn.o chebev_dyn.o join.o kinetic.o\
	L-BFGS-B.o latt_func.o lc_index.o link.o LinkedCellIndex.o\
	linmin_adjust_bonds.o LoadBspModuli.o low_up.o lsqqtn.o mappa.o\
	mapnl_divide.o matinv.o min_pack.o\
	minv.o mnbrakf.o Make_Pol_Models.o\
	mts_plot_fragm.o mts_plotp.o mts_test.o\
	mts_umbrl.o near0.o normal.o\
	parse.o parse_numbers.o\
	pbend.o pbond.o pckben.o pckfrg.o pckitor.o\
	pcktor.o phapp.o pitors.o plgndr.o plot_center.o plotc.o plot_fragm.o\
	plotd.o pntres.o pptors.o prtacc.o prtat.o\
	prtba.o prtbnd.o prtcn.o prtfrc.o prtit.o prtite.o prtite_min.o\
	print_title_analysis.o prtpt.o prtsq.o\
	prttors.o ranf.o ranvel.o rattle_verlet.o rattle_verlet_co.o\
	rattle_correc.o rattle_correc_co.o rdaprm.o rdatpg.o rdcmac.o\
	rdrgd.o read_analysis.o \
	read_inout.o read_integrator.o read_parameters.o\
	read_potential.o\
	read_properties.o read_run.o read_setup.o read_simulation.o\
	read_solute.o read_solvent.o read_space_group.o readco.o readtp.o\
	redefine_slv.o redmss.o remove_mv.o repl_coord.o repl_tpg.o\
	ResetCofM.o restim.o rigfit.o\
	rmbend.o rmpbond.o rotb.o RotCoordColumn.o RotCoordRow.o \
	sbend.o sbond.o ScalarSum.o scalar_sum_normal.o scalar_sum_transp.o \
	scalar_sum.o \
	scale_charges.o sctab.o search_clsth.o serchd.o set_const_co.o\
	set_hbonds_masks.o set_hhisto.o set_ss_array.o set_tempp.o\
	set_tempt.o setuns.o setup_fake_map.o Setup_fake_decompa.o\
	setup_skin_shell.o\
	spec14.o starting_verlet.o strblk.o\
	time_correlation.o timcor.o timcor0.o tors.o\
	tr_inbox.o tr_wsc.o trans_center.o twofft.o up_low.o update.o\
	uscrpl.o usort.o\
	verify_variables.o verlet.o verlet_free.o verlet_free_eta.o windows.o\
	wrenc.o write_contacts.o write_fragm_dist.o\
	write_pot_bond.o write_pot_nbond.o write_gyr.o write_gofrp.o\
	write_gofrw.o write_hbonds.o write_hbonds_vor.o write_hhisto.o\
	write_bends.o\
	write_bonds.o write_tors.o write_vacf_phi.o write_diffusion.o\
	write_voronoi_header.o write_rms.o write_sk.o\
	write_xrms.o write_xrms_atm.o xfit.o xerror.o\
	zero0.o zero.o zero3.o zero3x3.o zero_gofr.o zero_voronoi.o\
	zeroa.o zeroi.o write_mol_dipole.o write_electric_field.o\
	P_cm_split.o P_get_prot.o P_get_intra.o Pme_init.o\
	Polarization_Forces.o Do_Electric_Fields.o

OBJ_PP = P_omit_intra.o P_broadcast_i.o P_check_decomp.o\
	P_expand_r8.o P_expand_r8x3.o P_expand_i.o\
	P_fold_r8.o P_fold_r8x3.o P_fold_i.o P_get_errmsg.o \
	P_Reduce_Forces.o P_get_gofr.o \
	P_load_balance.o P_merge_i.o P_merge_r8.o\
	P_merge_vecr8.o P_merge_veci.o P_split_scalar0.o P_split_scalar.o\
	P_get_iret.o P_adjust_decomp.o P_Buffer_attach.o \
	P_Buffer_detach.o P_get_neighbors.o P_Write_neighbors.o


ifeq ($(MODE), PARALLEL)
OBJ = $(OBJ_NP) $(OBJ_PP)
else
OBJ = $(OBJ_NP)
endif

ifeq ($(TARGET0), Linux)
OBJTIM  = vstart.o timer_linux.o
else
OBJTIM = vstart.o
endif

PRECOMP_NP = analyse_voronoi.f bldbox.f start.f furippa.f insert_solute.f\
	mts_forpp.f comp_voronoi.f comp_abmd_fnative.f\
	write_cavities_vor.f \
	pme_code.f timer.f verify_input.f open_file_dump.f orac.f\
	read_input.f free_mem.f M_free.f M_free_s.f\
	run_minimize.f get_electric_field.f\
	get_memory_iobuffer.f get_total_energy.f\
	f1dim_der.f comp_total_energy.f\
	f1dim_total.f linmin_total.f mnbrakf_total.f brent_total.f\
	mts_intra_n0.f mts_intra_n1.f  mts_furier.f\
	mts_furier2.f\
	do_pmesh_kspace.f do_pmesh_kspace_dipole.f \
	write_confc.f read_confc_rows.f\
	read_confc_columns.f setup_box.f M_get_length.f\
	M_newlen.f M_memory.f memory.f M_memory_s.f fft_pme.f fft_pme_dipole.f\
	openf.f fft_setup.f fft_forward.f fft_back.f correc_stress.f\
	correc_stress_n0.f lc_list.f  lc_list2.f lc_list3.f mts_forces.f\
	mts_forces2.f mts_forpp2.f\
	read_parallel.f P_setup_decomp.f P_setup_decompa.f\
	P_atoms_split_intra.f P_split_intra.f P_split_constr.f\
	P_atoms_split_inter.f 

PRECOMP_NKL = mtsmd.f drive_analysis.f 

PRECOMP_PP = P_whoami.f P_open_input.f P_nwrite.f P_nread.f\
	P_npread.f \
	P_barrier.f \
	P_pme_assign_node.f 

PRE_F90 = P_Change_Decomposition.f90 P_Comm_Intra.f90 do_rfft3d.f90\
	HYDRATION_Mod.f90


PRECOMP_F90 = $(PRE_F90) 

ifeq ($(MODE), PARALLEL)
PRECOMP = $(PRECOMP_NP) $(PRECOMP_PP)
else
PRECOMP = $(PRECOMP_NP)
endif


PRECOMPT = $(patsubst %.f, %.CPP.f, $(PRECOMP))

PRECOMPT_NKL = $(patsubst %.f, %.CPP.f, $(PRECOMP_NKL))

PRECOMPT_F90 = $(patsubst %.f90, %.CPP.f90, $(PRECOMP_F90))

OBJPRE = $(patsubst %.f, %.o, $(PRECOMPT)) 

OBJPRE_NKL = $(patsubst %.f, %.o, $(PRECOMPT_NKL)) 

OBJPRE_F90 = $(patsubst %.f90, %.o, $(PRECOMPT_F90))

LIBCOMP = ffct4.f fft4rc.f

LIBCOMPT = $(patsubst %.f, %.LIB.f, $(LIBCOMP))

OBJLIB = $(patsubst %.f, %.o, $(LIBCOMPT))

ifeq ($(TARGET0), Linux)
OBJREC = recpr_linux.o perc_bondh_linux.o
OBJREC1 = recpr.o perc_bondh.o
else
OBJREC = recpr.o perc_bondh.o
OBJREC1 = recpr_linux.o perc_bondh_linux.o
endif

EIGRS = eigrs.o

#==========================================================================
# Uncomment the folling if an old restart has to be used
# 
#DUMP = dumprs_old.o readrs_old.o
#==========================================================================

DUMP = dumprs.o readrs.o readrs_old.o

#==========================================================================
# List of sources 
#==========================================================================

ifeq ($(TARGET0), Linux)
LIST = Makefile $(INCLUDE) $(patsubst %.o, %.f, $(OBJ_NP))\
	$(patsubst %.o, %.f, $(OBJ_PP)) $(PRECOMP_NP)\
	$(PRECOMP_PP) $(LIBCOMP) $(patsubst %.o, %.f, $(EIGRS))\
	$(patsubst %.o, %.c, $(OBJREC)) $(patsubst %.o, %.f, $(OBJREC1))\
	$(patsubst %.o, %.f, $(DUMP)) $(patsubst %.o, %.c, $(OBJTIM))\
	config.h_dist mtsmd_avg_inc.f mtsmd_dump_inc.f
else
LIST = Makefile $(INCLUDE) $(patsubst %.o, %.f, $(OBJ_NP))\
	$(patsubst %.o, %.f, $(OBJ_PP)) $(PRECOMP_NP)\
	$(PRECOMP_PP) $(LIBCOMP) $(patsubst %.o, %.f, $(EIGRS))\
	$(patsubst %.o, %.f, $(OBJREC)) $(patsubst %.o, %.c, $(OBJREC1))\
	$(patsubst %.o, %.f, $(DUMP)) $(patsubst %.o, %.c, $(OBJTIM))\
	timer_linux.c config.h_dist mtsmd_avg_inc.f mtsmd_dump_inc.f
endif


include Make.what

#==========================================================================
# Make Default 
#==========================================================================

all: 
	${MAKE} $(TARGET0)	


#==========================================================================
# Make include
#==========================================================================

parst.h: config.h parst.inc
	@$(CPP) parst.inc | awk '!/^ *$$/ && !/^#/' > temp
	./python/include_f90.py temp
	file=parst;\
	test -f $$file.h && ( cmp -s $$file.h temp || \
	cp temp $$file.h) || cp temp $$file.h ;\
	rm -f temp temp2 >/dev/null

vacf.h: config.h vacf.inc
	@$(CPP) vacf.inc | awk '!/^ *$$/ && !/^#/' > temp
	./python/include_f90.py temp
	file=vacf;\
	test -f $$file.h && ( cmp -s $$file.h temp || \
	cp temp $$file.h) || cp temp $$file.h ;\
	rm -f temp temp2 >/dev/null

mts_vacf.h: config.h mts_vacf.inc
	@$(CPP) mts_vacf.inc | awk '!/^ *$$/ && !/^#/' > temp
	./python/include_f90.py temp
	file=mts_vacf;\
	test -f $$file.h && ( cmp -s $$file.h temp || \
	cp temp $$file.h) || cp temp $$file.h ;\
	rm -f temp temp2 >/dev/null

lc_list.h: config.h lc_list.inc
	@$(CPP) lc_list.inc | awk '!/^ *$$/ && !/^#/' > temp
	./python/include_f90.py temp
	file=lc_list;\
	test -f $$file.h && ( cmp -s $$file.h temp || \
	cp temp $$file.h) || cp temp $$file.h ;\
	rm -f temp temp2 >/dev/null

mts_spectra.h: config.h mts_spectra.inc
	@$(CPP) mts_spectra.inc | awk '!/^ *$$/ && !/^#/' > temp
	./python/include_f90.py temp
	file=mts_spectra;\
	test -f $$file.h && ( cmp -s $$file.h temp || \
	cp temp $$file.h) || cp temp $$file.h ;\
	rm -f temp temp2 >/dev/null

analysis.h: config.h analysis.inc
	@$(CPP) -D$(MEMORY) analysis.inc | awk '!/^ *$$/ && !/^#/' > temp
	./python/include_f90.py temp
	file=analysis;\
	test -f $$file.h && ( cmp -s $$file.h temp || \
	cp temp $$file.h) || cp temp $$file.h ;\
	rm -f temp temp2 >/dev/null

resid.h: config.h resid.inc
	@$(CPP) resid.inc | awk '!/^ *$$/ && !/^#/' > temp
	./python/include_f90.py temp
	file=resid;\
	test -f $$file.h && ( cmp -s $$file.h temp || \
	cp temp $$file.h) || cp temp $$file.h ;\
	rm -f temp temp2 >/dev/null

fourier.h: config.h fourier.inc
	@$(CPP) fourier.inc | awk '!/^ *$$/ && !/^#/' > temp
	./python/include_f90.py temp
	file=fourier;\
	test -f $$file.h && ( cmp -s $$file.h temp || \
	cp temp $$file.h) || cp temp $$file.h ;\
	rm -f temp temp2 >/dev/null

pbc.h: config.h pbc.inc
	@$(CPP) pbc.inc | awk '!/^ *$$/ && !/^#/' > temp
	./python/include_f90.py temp
	file=pbc;\
	test -f $$file.h && ( cmp -s $$file.h temp || \
	cp temp $$file.h) || cp temp $$file.h ;\
	rm -f temp temp2 >/dev/null

pme.h: config.h pme.inc
	@$(CPP) -D$(TARGET1) -D$(MODE) -D$(FFTS) pme.inc | awk '!/^ *$$/ && !/^#/' > temp
	./python/include_f90.py temp
	file=pme;\
	test -f $$file.h && ( cmp -s $$file.h temp || \
	cp temp $$file.h) || cp temp $$file.h ;\
	rm -f temp temp2 >/dev/null

parallel.h: parallel.inc
	@$(CPP) -D$(TARGET1) -D$(MODE) parallel.inc | awk '!/^ *$$/ && !/^#/' > temp
	./python/include_f90.py temp
	file=parallel;\
	test -f $$file.h && ( cmp -s $$file.h temp || \
	cp temp $$file.h) || cp temp $$file.h ;\
	rm -f temp temp2 >/dev/null

orac.h: orac.inc
	@$(CPP) -D$(TARGET0) orac.inc | awk '!/^ *$$/ && !/^#/' > temp
	./python/include_f90.py temp
	file=orac;\
	test -f $$file.h && ( cmp -s $$file.h temp || \
	cp temp $$file.h) || cp temp $$file.h ;\
	rm -f temp temp2 >/dev/null


cpropar.h: cpropar.inc
	@$(CPP) -D$(TARGET0) cpropar.inc | awk '!/^ *$$/ && !/^#/' > temp
	./python/include_f90.py temp
	file=cpropar;\
	test -f $$file.h && ( cmp -s $$file.h temp || \
	cp temp $$file.h) || cp temp $$file.h ;\
	rm -f temp temp2 >/dev/null

iobuffer.h: iobuffer.inc
	@$(CPP) iobuffer.inc | awk '!/^ *$$/ && !/^#/' > temp
	./python/include_f90.py temp
	file=iobuffer;\
	test -f $$file.h && ( cmp -s $$file.h temp || \
	cp temp $$file.h) || cp temp $$file.h ;\
	rm -f temp temp2 >/dev/null

voronoi.h: config.h voronoi.inc
	@$(CPP) -D$(MEMORY) voronoi.inc | awk '!/^ *$$/ && !/^#/' > temp
	./python/include_f90.py temp
	file=voronoi;\
	test -f $$file.h && ( cmp -s $$file.h temp || \
	cp temp $$file.h) || cp temp $$file.h ;\
	rm -f temp temp2 >/dev/null

unit.h: unit.inc
	@$(CPP) -D$(TARGET0) unit.inc | awk '!/^ *$$/ && !/^#/' > temp
	./python/include_f90.py temp
	file=unit;\
	test -f $$file.h && ( cmp -s $$file.h temp || \
	cp temp $$file.h) || cp temp $$file.h ;\
	rm -f temp temp2 >/dev/null

parameters.h: parameters.inc
	@$(CPP) -D$(TARGET0) parameters.inc | awk '!/^ *$$/ && !/^#/' > temp
	./python/include_f90.py temp
	file=parameters;\
	test -f $$file.h && ( cmp -s $$file.h temp || \
	cp temp $$file.h) || cp temp $$file.h ;\
	rm -f temp temp2 >/dev/null

mpi_size.h: mpi_size.inc
	@$(CPP) -D$(MPI_TYPE) mpi_size.inc | awk '!/^ *$$/ && !/^#/' > temp
	./python/include_f90.py temp
	file=mpi_size;\
	test -f $$file.h && ( cmp -s $$file.h temp || \
	cp temp $$file.h) || cp temp $$file.h ;\
	rm -f temp temp2 >/dev/null


#==========================================================================
# Make precompiled files
#==========================================================================

%.CPP.f: %.f
	@${CPP} -P -D$(TARGET0) -D$(TARGET1) -D$(MODE) -D$(PARA_MODE) -D$(MEMORY) -D$(FFTS) $< | awk '!/^ *$$/ && !/^#/' > temp;\
	file=$@; \
	test -f $$file && ( cmp -s $$file temp || \
	cp temp $$file) || cp temp $$file ;\
	rm -f temp temp2 >/dev/null

%.CPP.f90: %.f90
	@${CPP} -P -D$(TARGET0) -D$(TARGET1) -D$(MODE) -D$(PARA_MODE) -D$(MEMORY) -D$(FFTS) $< | awk '!/^ *$$/ && !/^#/' > temp;\
	file=$@; \
	test -f $$file && ( cmp -s $$file temp || \
	cp temp $$file) || cp temp $$file ;\
	rm -f temp temp2 >/dev/null

$(OBJPRE): %.CPP.o: %.f


$(OBJPRE_NKL): %.CPP.o: %.f

$(OBJPRE_F90): %.CPP.o: %.f90

%.LIB.f: %.f
	@${CPP} -P -D$(LIB) $< | awk '!/^ *$$/ && !/^#/' > temp;\
	file=$@; \
	test -f $$file && ( cmp -s $$file temp || \
	cp temp $$file) || cp temp $$file ;\
	rm -f temp temp2 >/dev/null

$(OBJLIB): %.LIB.o: %.f

#==========================================================================
# Decide what to make
#==========================================================================


PGI :
	${MAKE} obj_PGI "FFLAGS = -O2 "
	${MAKE} prepro_PGI "FFLAGS = -O2" 
	${MAKE} prepronkl_PGI "FFLAGS = -O0 "
	${MAKE} eigr_PGI "FFLAGS = -O2"
	${MAKE} recur_PGI "FFLAGS =  "
	${MAKE} orac_PGI

INTEL :
	${MAKE} obj_INTEL "FFLAGS = -w95 -O2 " 
	${MAKE} prepro_INTEL "FFLAGS = -w95 -O3" 
	${MAKE} prepronkl_INTEL "FFLAGS = -warn noalignments -O0 "
	${MAKE} eigr_INTEL "FFLAGS = -w95 -O1 "
	${MAKE} recur_INTEL "FFLAGS = -w95  "
	${MAKE} orac_INTEL 

OSF1 :
	${MAKE} obj_OSF "FFLAGS = -fast " 
	${MAKE} prepro_OSF "FFLAGS = -O5 -fast" 
	${MAKE} prepronkl_OSF "FFLAGS = -warn noalignments -O0 "
	${MAKE} eigr_OSF "FFLAGS = -O1 "
	${MAKE} recur_OSF "FFLAGS = -recursive "
	${MAKE} orac_OSF 

XLF :
	${MAKE} obj_XLF "FFLAGS = -O5 -qnosave -qhot -qipa" 
	${MAKE} prepro_XLF "FFLAGS = -O5 -qnosave -qhot -qipa" 
	${MAKE} prepronkl_XLF "FFLAGS = -O0 -qsave -qstrict "
	${MAKE} eigr_XLF "FFLAGS = -O5 -qsave -qipa "
	${MAKE} recur_XLF "FFLAGS = -O4 -qnosave -qipa"
	${MAKE} orac_XLF 

G95 :
	${MAKE} obj_G95 
	${MAKE} prepro_G95
	${MAKE} prepronkl_G95
	${MAKE} eigr_G95
	${MAKE} recur_G95
	${MAKE} orac_G95 


#========================================================================

DISTRIBUTION : 
	@rm -r -f $(CURDIR)/tmp_dir
	mkdir $(CURDIR)/tmp_dir ;\
	mkdir $(CURDIR)/tmp_dir/orac ;\
	mkdir $(CURDIR)/tmp_dir/orac/source ;\
	ln -s $(CURDIR)/RCS $(CURDIR)/tmp_dir/orac/source/RCS ;\
	cd $(CURDIR)/tmp_dir/orac/source ; co $(LIST) ; rm RCS ;\
	mv config.h_dist config.h;\
	chmod 644 * ; cd ../

#	for j in $(LIST);\
#	do cp $(CURDIR)/$$j $(CURDIR)/$(DIST_DIR)/$$j;done;\
#	cp $(CURDIR)/config.h $(CURDIR)/$(DIST_DIR)/config.h;\
#	chmod 600 $(CURDIR)/$(DIST_DIR)/* ; cd ../



#==========================================================================
# Compile and load ORAC on OSF   
#==========================================================================


obj_OSF : $(OBJ) $(DUMP) $(OBJLIB) $(OBJTIM)

recur_OSF: $(OBJREC)

eigr_OSF: $(EIGRS)

prepro_OSF: $(OBJPRE) $(OBJPRE_F90)

prepronkl_OSF: $(OBJPRE_NKL)

orac_OSF : $(OBJ) $(OBJREC) $(DUMP) $(OBJPRE) \
	$(OBJPRE_NKL) $(OBJPRE_F90) $(OBJLIB) $(EIGRS) $(OBJTIM)
	$(LD) -o orac_OSF $(OBJ) $(OBJREC) $(DUMP) $(OBJPRE) \
	$(OBJPRE_NKL) $(OBJPRE_F90) $(OBJLIB) $(EIGRS) $(OBJTIM)\
	-L/usr/local/assist/lib -L$(FFTW_LIB) $(LLFLAGS) -lmpi

#==========================================================================
# Compile and load ORAC on OSF   
#==========================================================================


obj_PGI : $(OBJ) $(DUMP) $(OBJLIB) $(OBJTIM)

recur_PGI: $(OBJREC)

eigr_PGI: $(EIGRS)

prepro_PGI: $(OBJPRE) $(OBJPRE_F90) 

prepronkl_PGI: $(OBJPRE_NKL)

orac_PGI : $(OBJ) $(OBJREC) $(DUMP) $(OBJPRE) $(OBJPRE_NKL)\
	$(OBJPRE_F90) $(OBJLIB) $(EIGRS) $(OBJTIM)
	$(LD) -Mvect -Mipa=fast -fastsse  -o orac_PGI $(OBJ) $(OBJREC) $(DUMP) $(OBJPRE) $(OBJPRE_NKL)\
	$(OBJPRE_F90) $(OBJLIB) $(EIGRS) $(OBJTIM)\
	-L$(FFTW_LIB) $(LLFLAGS) -lm


#==========================================================================
# Compile and load ORAC on INTEL 
#==========================================================================

obj_INTEL : $(OBJ) $(DUMP) $(OBJLIB) $(OBJTIM)

recur_INTEL: $(OBJREC)

eigr_INTEL: $(EIGRS)

prepro_INTEL: $(OBJPRE) $(OBJPRE_F90) 

prepronkl_INTEL: $(OBJPRE_NKL)

orac_INTEL: $(OBJ) $(OBJREC) $(DUMP) $(OBJPRE) $(OBJPRE_NKL)\
	$(OBJPRE_F90) $(OBJLIB) $(EIGRS) $(OBJTIM)
	$(LD) -o orac_INTEL $(OBJ) $(OBJREC) $(DUMP) $(OBJPRE) $(OBJPRE_NKL)\
	$(OBJPRE_F90) $(OBJLIB) $(EIGRS) $(OBJTIM)\
	-L$(FFTW_LIB) $(LLFLAGS)

#==========================================================================
# Compile and load ORAC on XLF compiler
#==========================================================================

obj_XLF : $(OBJ) $(DUMP) $(OBJLIB) $(OBJTIM)

recur_XLF: $(OBJREC)

eigr_XLF: $(EIGRS)

prepro_XLF: $(OBJPRE) $(OBJPRE_F90) 

prepronkl_XLF: $(OBJPRE_NKL)

orac_XLF: $(OBJ) $(OBJREC) $(DUMP) $(OBJPRE) $(OBJPRE_NKL)\
	$(OBJPRE_F90) $(OBJLIB) $(EIGRS) $(OBJTIM)
	$(LD) -O5 -qnosave -qhot -qipa -o orac_XLF $(OBJ) $(OBJREC)\
	$(DUMP) $(OBJPRE) $(OBJPRE_NKL)\
	$(OBJPRE_F90) $(OBJLIB) $(EIGRS) $(OBJTIM)\
	-L$(FFTW_LIB) $(LLFLAGS)

#==========================================================================
# Compile and load ORAC on G95 compiler
#==========================================================================

obj_G95 : $(OBJ) $(DUMP) $(OBJLIB) $(OBJTIM)

recur_G95: $(OBJREC)

eigr_G95: $(EIGRS)

prepro_G95: $(OBJPRE) $(OBJPRE_F90) 

prepronkl_G95: $(OBJPRE_NKL)

orac_G95: $(OBJ) $(OBJREC) $(DUMP) $(OBJPRE) $(OBJPRE_NKL)\
	$(OBJPRE_F90) $(OBJLIB) $(EIGRS) $(OBJTIM)
	$(LD) -O3 -o orac_G95 $(OBJ) $(OBJREC)\
	$(DUMP) $(OBJPRE) $(OBJPRE_NKL)\
	$(OBJPRE_F90) $(OBJLIB) $(EIGRS) $(OBJTIM)\
	-L$(FFTW_LIB) $(LLFLAGS)
#==========================================================================
# Clean Objects
#==========================================================================

clean :
	rm -f *.o *.CPP.f *.CPP.f90 *.LIB.f *.mod

#==========================================================================
# Clean Objects and executable
#==========================================================================

extra_clean :
	rm -f *.o *.CPP.f *.LIB.f orac_INTEL orac_OSF orac_PGI
	@cp config.h temp ; rm -f *.h ; mv temp config.h > /dev/null

#==========================================================================
# Get all *.f files,  checkin all files, tar 
#==========================================================================

get_orac :
	co -l $(LIST)

checkout :
	co -l $(LIST)

checkin : 
	echo $(DESCRIPTION) > release.tmp
	ci -r$(RELEASE) -f -l -m"$(STATUS)"  -trelease.tmp $(LIST)
	rm -f release.tmp

tar :   
	tar zcvf orac.tar.gz $(LIST)

#==========================================================================
# Dependencies 
#==========================================================================

accdif.o: parst.h pbc.h
add_solvent_tpg.o: parst.h cpropar.h
adjust_solvent.o: parst.h
adjust_bonds.o: parst.h
analyse_voronoi.CPP.o:analyse_voronoi.CPP.f
analyse_voronoi.CPP.o:parst.h voronoi.h
bcnstp.o: pbc.h
bldbox.CPP.o:bldbox.CPP.f cpropar.h parst.h
setup_box.CPP.o: setup_box.CPP.f cpropar.h parst.h unit.h
blkdta.o: unit.h cpropar.h parst.h
brot.o: unit.h
CalcXrmsSecStruct.o: parst.h
calc_xrms.o: parst.h
calc_xslt.o: parst.h
calc_lda.o: parst.h
calc_hyd.o: pbc.h
calc_lda_rest.o: pbc.h
calc_lda_avg.o: parst.h
calc_lda_hyd.o: pbc.h
calc_inst_fit.o: parst.h
calc_RotMat.o: parst.h
calc_avg_str.o: parst.h
calc_avg2_str.o: parst.h
calc_cofm.o: parst.h
calc_gofr.o: parst.h pbc.h
calc_sk.o: unit.h parst.h
calc_sk_eet.o: unit.h parst.h
cewct.o: unit.h
change_tors.o: parameters.h parst.h
chkewl.o: fourier.h unit.h
claslv.o: parst.h resid.h pbc.h
close_files.o: unit.h
collision.o: unit.h parst.h
comp_abmd_fcryst.o: parst.h
comp_abmd_fnative.CPP.o: comp_abmd_fnative.CPP.f parst.h
comp_contacts.o: pbc.h
comp_native.o: pbc.h
comp_neigh_vor.o: parst.h pbc.h voronoi.h
comp_cavities_vor.o: parst.h voronoi.h pbc.h analysis.h
comp_dip.o: pbc.h
comp_dynamic_mat.o: parst.h cpropar.h unit.h
comp_total_energy.CPP.o: parst.h cpropar.h unit.h
comp_thermos_energy.o: unit.h
comp_thermos_forces.o: unit.h
comp_voronoi.CPP.o: pbc.h parst.h voronoi.h
correc_stress.CPP.o: unit.h
correc_stress_n0.CPP.o: unit.h
cov_thermos.o: unit.h
covbd.o: unit.h
covchg.o: unit.h
covfod.o: unit.h
covhov.o: unit.h
covnb.o: unit.h
covpr.o: unit.h
covpw.o: unit.h
covsol.o: unit.h
covsol_int.o: unit.h
covuen.o: unit.h
cpfld.o: unit.h pbc.h
crdrtf.o: parst.h cpropar.h unit.h
cself.o: unit.h
cself_dipole.o: unit.h
delrfg.o: parst.h
drive_analysis.CPP.o: drive_analysis.CPP.f cpropar.h unit.h fourier.h\
	parst.h iobuffer.h lc_list.h voronoi.h analysis.h HYDRATION_Mod.CPP.o\
	DENSITY_Mod.o CENTER_SOL_Mod.o PDBs_Mod.o RMS_Subtract_Mod.o
dspat.o: parst.h pbc.h
dumptp.o: cpropar.h parst.h
erf_corr_cutoff.o: unit.h
errmrk.o: unit.h
errprt.o: unit.h
ewself_slv.o: unit.h
f1dim_adjust_bonds.o: cpropar.h parst.h
f1dim_der.CPP.o: cpropar.h parst.h
f1dim_total.CPP.o: cpropar.h parst.h
ferrf.o: unit.h
ferrf_tag.o: unit.h
ffct4.o: DXMLDEF.FOR
fft4rc.o: fft4rc.LIB.f DXMLDEF.FOR
fft_pme.CPP.o: fft_pme.CPP.f pme.h
fft_pme2.CPP.o: fft_pme2.CPP.f pme.h
fhbondb.o: parst.h pbc.h
find_igint_vor.o: parst.h voronoi.h
find_length_run.o: unit.h
fitors.o: unit.h 
fnb14.o: unit.h
fnb14_2.o: unit.h
fnbgrp.o: unit.h pbc.h
fnbgrp2.o: unit.h pbc.h polar.h
furippa.CPP.o: furippa.CPP.f fourier.h parst.h unit.h
gauss.o: unit.h
gaussc.o: unit.h
gausse.o: unit.h
getchr.o: unit.h
get_total_energy.CPP.o: get_total_energy.CPP.f parst.h cpropar.h unit.h\
	parallel.h
grphpw.o: unit.h
graphw.o: unit.h
insert_solute.CPP.o: insert_solute.CPP.f unit.h pbc.h
join.o: parst.h parameters.h unit.h
kinetic.o: unit.h parst.h
lc_list.CPP.o: pbc.h parst.h lc_list.h
lc_list2.CPP.o: pbc.h parst.h lc_list.h
lc_list3.CPP.o: pbc.h parst.h lc_list.h
mappc.o: parst.h
mts_fnb14.o: unit.h pbc.h
mts_forces.CPP.o: mts_forces.CPP.f unit.h cpropar.h parst.h mpi_size.h
mts_forces2.CPP.o: mts_forces2.CPP.f unit.h cpropar.h parst.h
mts_forpp.CPP.o: mts_forpp.CPP.f mts_forpp.f force_ele_nospline.f \
	force_ele_spline.f force_noele.f force_ele_nospline_sw.f\
	force_ele_spline_sw.f force_noele_sw.f pressure.h\
	mts_forpp.CPP.f parst.h unit.h pbc.h
mts_forpp2.CPP.f: pressure.h polar.h
mts_forpp2.CPP.o: mts_forpp2.CPP.f parst.h unit.h pbc.h polar.h
mts_furier.CPP.o: mts_furier.CPP.f parst.h cpropar.h pme.h
mts_furier2.CPP.o: mts_furier2.CPP.f parst.h cpropar.h pme.h
mts_intra_n0.CPP.o: cpropar.h parst.h
mts_intra_n1.CPP.o: cpropar.h parst.h
mts_plot_fragm.o: unit.h
plot_fragm.o: unit.h
mts_plotp.o: unit.h
mts_plots.o: unit.h
mts_spectra.o: mts_spectra.h
mtsmd.CPP.o: mtsmd.CPP.f cpropar.h unit.h fourier.h \
	parst.h pme.h lc_list.h parallel.h mtsmd_avg_inc.CPP.f\
	mtsmd_dump_inc.CPP.f Module_Extra_Forces.o
orac.CPP.o: orac.CPP.f parst.h cpropar.h orac.h unit.h parallel.h
open_file_dump.CPP.o: open_file_dump.CPP.f unit.h

P_get_neighbors.o: pbc.h
P_split_intra.CPP.o: parst.h cpropar.h
P_split_constr.CPP.o: parst.h cpropar.h
P_open_input.CPP.o: P_open_input.CPP.f
P_pme_assign_node.CPP.o: P_pme_assign_node.CPP.f

pckben.o: parst.h
pcktor.o: parst.h
pckitor.o: parst.h
perc_bondh.o: parst.h cpropar.h
pitors.o: parst.h
plot.o: unit.h
plot_center.o: unit.h
plotb.o: unit.h
plotc.o: unit.h
pltbck.o: unit.h
pme_code.CPP.o: pme_code.CPP.f
ifeq ($(shell uname), HP-UX)
	f77 +OP4 +U77 -c pme_code.CPP.f -o pme_code.CPP.o
endif
pothoov.o: unit.h
pptors.o: parst.h
prtacc.o: parst.h cpropar.h unit.h
prtat.o: unit.h 
prtbnd.o: unit.h
prtba.o: unit.h 
prtcn.o: unit.h 
prtecf.o: unit.h
prtit.o: unit.h 
prtite.o: parst.h cpropar.h unit.h
prtite_min.o: parst.h cpropar.h unit.h
print_title_analysis.o: parst.h cpropar.h unit.h
prtpt.o: unit.h pbc.h
prtrms.o: unit.h
prtsq.o: unit.h
prttors.o: unit.h parst.h pbc.h
pstopm.o: pbc.h
ranvel.o: unit.h parst.h
rdcmac.o: parst.h
rdcnc.o: unit.h
rdcsol.o: unit.h
rdunit.o: unit.h
read_confc_rows.CPP.o: unit.h
read_confc_columns.CPP.o: unit.h
readco.o: unit.h
readtp.o: cpropar.h parst.h
recpr.o: parst.h cpropar.h
repl_coord.o: unit.h
repl_tpg.o: parst.h cpropar.h
restim.o: resid.h unit.h
rmbend.o: parst.h
rotb.o: unit.h
rattle_verlet.o: parst.h
rattle_correc.o: parst.h
read_input.CPP.o: read_input.CPP.f parst.h cpropar.h unit.h
read_run.o: parst.h cpropar.h unit.h
read_inout.o: parst.h cpropar.h unit.h
read_simulation.o: parst.h cpropar.h unit.h Module_Stress.o
read_potential.o: parst.h cpropar.h unit.h
read_setup.o: parst.h cpropar.h unit.h
read_integrator.o: parst.h cpropar.h unit.h
read_solute.o: parst.h cpropar.h unit.h
read_properties.o: parst.h cpropar.h unit.h HYDRATION_Mod.CPP.o DENSITY_Mod.o\
	PDBs_Mod.o RMS_Subtract_Mod.o
read_solvent.o: parst.h cpropar.h unit.h
read_analysis.o: parst.h cpropar.h unit.h CENTER_SOL_Mod.o 
PDBs_Mod.o: PDBs_Read.f90 INPUT_Mod.o 
DENSITY_Mod.o: DENSITY_Read.f90 INPUT_Mod.o
RMS_Subtract_Mod.o: RMS_Subtract_Read.f90
read_parameters.o: parst.h parameters.h 
read_parallel.CPP.o: parst.h cpropar.h unit.h
rigfit.o: parst.h
run_minimize.CPP.o: run_minimize.CPP.f cpropar.h unit.h parst.h\
	pme.h lc_list.h parallel.h
savech.o: unit.h
setuns.o: unit.h
set_tempp.o: unit.h Module_Stress.o
set_tempt.o: unit.h
ssitep.o: parst.h
start.CPP.o: start.CPP.f unit.h cpropar.h parst.h
strpco.o: pbc.h
timcor.o: unit.h vacf.h
timcor0.o: vacf.h
tmphvr.o: unit.h
tmppr.o: unit.h
tmppr_iso.o: unit.h
tr_inbox.o: pbc.h
mts_umbrl.o: pbc.h parst.h
update.o: parst.h unit.h pbc.h
vstart.o: config.h
verify_input.CPP.f: pressure.h
verify_input.CPP.o: parst.h cpropar.h unit.h pme.h
verify_variables.o: parst.h cpropar.h unit.h
write_fragm_dist: pbc.h
write_cavities_vor.CPP.o: write_cavities_vor.CPP.f
write_cavities_vor.CPP.o: parst.h
write_gyr.o: parst.h
write_gofrp.o: unit.h
write_gofrw.o: unit.h
write_hbonds.o: parst.h pbc.h
write_hbonds_vor.o: parst.h voronoi.h
write_hhisto.o: unit.h
write_bonds.o: parst.h
write_bends.o: parst.h
write_confc.CPP.o: unit.h
write_pot_bond.o: unit.h parst.h
write_pot_nbond.o: unit.h parst.h
write_tors.o: parst.h
write_diffusion.o: unit.h
write_electric_field.o: parst.h cpropar.h unit.h
write_vacf_phi.o: unit.h
write_xrms.o:parst.h
xerror.o: unit.h
zero_voronoi.o: parst.h voronoi.h
M_memory.CPP.o: M_memory.CPP.f
M_memory_s.CPP.o: M_memory_s.CPP.f
write_confc.CPP.o: write_confc.CPP.f
read_confc_rows.CPP.o: read_confc_rows.CPP.f
openf.CPP.o:openf.CPP.f
timer.CPP.0: timer.CPP.f
P_atoms_split_intra.CPP.o: P_atoms_split_intra.CPP.f
P_setup_decomp.CPP.o: P_setup_decomp.CPP.f
M_get_length.CPP.o: M_get_length.CPP.f
get_length_iobuffer.CPP.o:get_length_iobuffer.CPP.f
fft_back.CPP.o: fft_back.CPP.f
fft_forward.CPP.o: fft_forward.CPP.f
gpfft.CPP.o: gpfft.CPP.f
do_pmesh_kspace.CPP.o: do_pmesh_kspace.CPP.f do_rfft3d.CPP.o
do_pmesh_kspace_dipole.CPP.o: do_pmesh_kspace_dipole.CPP.f do_rfft3d.CPP.o
do_pmesh_kspace2.CPP.o: do_pmesh_kspace2.CPP.f
free_mem.CPP.o: free_mem.CPP.f
memory.CPP.o: memory.CPP.f
P_Change_Decomposition.CPP.o: P_Change_Decomposition.CPP.f90
P_Comm_Intra.CPP.o: P_Comm_Intra.CPP.f90 mpi_size.h
do_rfft3d.o: xerror_mod.o
get_memory_iobuffer.CPP.o: get_memory_iobuffer.CPP.f iobuffer.h
P_adjust_decomp.CPP.o: P_adjust_decomp.CPP.f mpi_size.h
P_expand_i.o: mpi_size.h
P_expand_r8.o: mpi_size.h
P_expand_r8x3.o: mpi_size.h
P_fold_i.o: mpi_size.h
P_fold_r8.o: mpi_size.h
P_fold_r8x3.o: mpi_size.h
P_get_errmsg.o: mpi_size.h
P_get_iret.o: mpi_size.h
P_merge_i.o: mpi_size.h
P_merge_r8.o: mpi_size.h
P_merge_veci.o: mpi_size.h
P_merge_vecr8.o: mpi_size.h
P_open_input.o:  mpi_size.h
P_Reduce_Forces.o: mpi_size.h
P_Write_neighbors.o: mpi_size.h
P_broadcast_i.o: mpi_size.h
do_rfft3d.CPP.o: do_rfft3d.CPP.f90 
Pme_init.o: do_rfft3d.CPP.o
P_atoms_split_intra.CPP.o: P_atoms_split_intra.CPP.f
Gauss_Overlap_Mod.o: Neighbors_Mod.o
Polarization_Forces.o: Gauss_Overlap_Mod.o
correc3x3.o: Module_Stress.o
find_thermos.o: Module_Stress.o
rattle_verlet_co.o: Module_Stress.o
rattle_correc_co.o: Module_Stress.o
kinetic.o: Module_Stress.o
Module_Extra_Forces.o: unit.h
HYDRATION_Mod.CPP.o: HYDRATION_Mod.CPP.f90 mpi_size.h
