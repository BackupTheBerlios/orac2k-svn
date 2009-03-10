# Makefile for ORAC. Works for sure with the intel compiler ifort and gfortran 
# supposedly with pgf90 
# 
#
#   "make" or "make MODE=PARALLEL" for parallel systems
#   "make MODE=SEQUENTIAL" for single processor runs
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
DIST_DIR = tmp_dir/orac2k/source
LIBRARIES = NONE

MODE=PARALLEL
MEMORY=DYNAMIC_MEM
MPI_TYPE = OPENMPI
PARA_MODE = MPI

ifneq ($(MODE), PARALLEL)
PARA_MODE = NOMPI
endif 

export MODE
export PARA_MODE
export MEMORY 
export MPI_TYPE

##
##  Machine dependent include files
##  change fftw and lapack/blas library to suit
##  your system
##

#include Make.Mac_g95_dbg
#include Make.Mac_g95
#include Make.Mac_xlf
#include Make.Mac_intel
#include Make.Lx_pgi_64
#include Make.Lx_ifc_argent
#include Make.Lx_ifc_intel
#include Make.Lx_ifc_intel.err

include Make.Lx_ifc_64_breizh

#include Make.Lx_gfortran_64_breizh

#include Make.Lx_ifc_64
#include Make.Lx_ifc_platine
#include Make.Lx_ifc_64_dbg
#include Make.OSF_f90


#==========================================================================
# Define default rules
#==========================================================================

.f.o :
	$(FC) $(FFLAGS) -D$(TARGET0) -D$(TARGET2) -D$(MODE) -D$(MEMORY) -D$(PARA_MODE) -c $<
.SUFFIXES: .o .f90
.f90.o :
	$(F90) $(FFLAGS) -D$(TARGET0) -D$(TARGET2) -D$(MODE) -D$(MEMORY) -D$(PARA_MODE) -c $<

#==========================================================================
# Object modules
#==========================================================================
OBJ_F90 = CalcDistConf.o \
	CENTER_SOL_Mod.o Class_Bending.o Class_Connect.o Class_Gauss_Param.o	\
	Class_Gcharges.o COMM_Mod.o CompElecPotentialOnGrid.o Correct_Ewald.o	\
	DCopy_V3.o DENSITY_Mod.o DENSITY_Read.o Do_Electric_Fields.o		\
	DOMAINS_Mod.o do_pfft3d.o do_rfft3d.o ElecPotential_Mod.o EUL_Mod.o	\
	Extra_Forces.o FIELD_Mod.o FIELD__Read.o FOURIER_Mod.o FOURIER_Read.o	\
	Gauss_Overlap_Mod.o GEOM_groups_Mod.o GEOM_groups_Readit.o		\
	GetHardSphere.o GetKineticTest.o Get_vp_without_vc.o Get_vp_with_vc.o	\
	GROUPS_Mod.o GROUPS_Readit.o HEATING_Mod.o HYDRATION_Mod.o		\
	HYDYNAMICS_Mod.o Init_Fixed_Molecule.o INPUT_Mod.o IOBUFFER_Mod.o	\
	LinkedCellIndex.o LoadBspModuli.o Make_Pol_Models.o MD_Mod.o		\
	Mod_Fourier.o Module_Extra_Forces.o Module_Fourier.o			\
	Module_Neighbors.o Module_Node.o			\
	Module_Orthogonal_Stress.o Module_Stress.o Module_Thole.o		\
	Neighbors.o Neighbors_Mod.o Open_Input.o		\
	ORAC_Finalize.o PBC_Mod.o P_Buffer.o P_Buffer_Mod.o			\
	P_Change_Decomposition.o P_Comm_Intra.o PDB.o PDB_FRAGM_Mod.o		\
	PDB_FRAGM__Read.o PDB__Read.o PDBs_Mod.o PDBs_Read.o PDB_utils.o	\
	P_expand_i.o P_expand_r8.o P_expand_r8x3b.o P_expand_r8x3.o P_fold_i.o	\
	P_fold_r8.o P_fold_r8x3.o P_isend.o P_isend_i.o Pme_init.o P_merge_i.o	\
	P_merge_r8.o P_merge_veci.o P_merge_vecr8.o Poisson_Boltzmann.o		\
	Polarization_Forces.o POLAR_Mod.o POLAR__Read.o P_pme_init.o		\
	P_recv_i.o P_Reduce_Forces.o P_Write_neighbors.o REDUCE_Mod.o		\
	ResetCofM.o reset_nbond.o RIGID_Mod.o RMS_General.o RMS_General_Read.o	\
	RMS_Matrix_Mod.o RMS_Subtract_Mod.o RMS_Subtract_Read.o scalar_sum.o	\
	ScalarSum.o SCALARSUM_Mod.o self_dipole.o Setup_fake_decompa.o		\
	Setup_fake_decomp.o test.o Tetrahedra.o VORONOI_Mod.o	\
	VORONOI_Read.o write_mol_dipole.o WSC_Mod.o WSC__Read.o xerror_mod.o

OBJ_F = add1h.o add2h.o add3h.o add_char1.o add_concta.o		\
	add_energies.o addh2o.o addh.o add_int1.o add_int2.o addlnk.o		\
	add_mapnl.o add_real1.o add_real2.o add_solvent_coord.o			\
	add_solvent_tpg.o add_str.o addtor.o addtpg.o adjust_bonds.o		\
	analyse_voronoi.o anneal.o appatm.o appbkn.o appbou.o appct.o appgrp.o	\
	asng_xrmsw.o asstpe.o atoms.o atthd.o bcnstp.o bend.o bldbox.o		\
	blkdta.o bond.o brent.o brent_total.o brot.o calc_avg2_str.o		\
	calc_avg_str.o calc_avg_xrms.o calc_cofm.o calc_gofr.o calc_hyd.o	\
	calc_inst_fit.o calc_lda_avg.o calc_lda.o calc_lda_flu.o		\
	calc_lda_hyd.o calc_lda_rest.o calc_RotMat.o calc_sk_eet.o calc_sk.o	\
	calc_xrms.o CalcXrmsSecStruct.o calc_xslt.o caldr.o cc.o		\
	change_coord_ex.o change_coord_inner.o change_frame.o change_origin.o	\
	change_step_sim.o change_tors.o change_virtual_potential.o chder_dyn.o	\
	chebev_dyn.o chebft_dyn.o check_length_sim.o check_read_columns.o	\
	check_topology.o check_zero_coord.o chk14.o chkewl.o chlab.o chmass.o	\
	chnpr.o close_files.o clpcm.o clt14.o collision.o comd0.o		\
	comp_abmd_fcryst.o comp_abmd_fdiss.o comp_abmd_fnative.o		\
	comp_abmd_tors.o comp_ang_vel.o comp_cavities_vor.o comp_cell_center.o	\
	comp_concth.o comp_contacts.o comp_dip.o comp_dynamic_mat.o comp_fcm.o	\
	comp_forcep.o comp_gmgp.o comp_molmass.o comp_native.o			\
	comp_neigh_vor.o comp_nmol.o comp_noneq_free.o comp_rmsq.o		\
	comp_stress_bnd.o comp_stress_conf.o comp_stress_kinetic.o		\
	comp_thermos_energy.o comp_thermos_forces_cm.o comp_thermos_forces.o	\
	comp_total_energy.o comp_vcm.o comp_vel_labframe.o comp_voronoi.o	\
	coordinate_spline.o coordinate_spline_init.o copy_protl.o correc3x3.o	\
	correc_etap.o correc_exp_scale.o correc.o correc_gmgp.o correc_matr.o	\
	correc_scale_2.o correc_scale_co.o correc_scale.o correc_stress.o	\
	correc_stress_n0.o covbd.o covchg.o covfod.o covnb.o covpr.o		\
	cov_thermos.o covuen.o crdrtf.o cself_dipole.o cself.o cself_phi.o	\
	csort.o dblas.o dcopy3.o delrfg.o dipole_interaction.o dist_ijk.o	\
	dlineq.o do_pmesh_kspace2.o do_pmesh_kspace_dipole.o do_pmesh_kspace.o	\
	drive_analysis.o dspat.o dumprs.o dumptp.o duni.o dupl_char1.o		\
	dupl_concta.o dupl_int1.o dupl_int2.o dupl_mapnl.o dupl_real1.o		\
	dupl_real2.o dupl_reali2.o eigrs.o equi_bend.o equi_bond.o		\
	equi_impro.o equi_tors.o erf_corr_cutoff.o erfc_spline.o		\
	erfc_spline_init.o f1dim_adjust_bonds.o f1dim_der.o f1dim_total.o	\
	ferrf.o ferrf_tag.o ffct4.o fft3-public.o fft4rc.o fft_back.o		\
	fft_default.o fft_forward.o fftpack.o fft_pme2.o fft_pme_dipole.o	\
	fft_pme.o fft_pme_init.o fft_setup.o fill_charge_grid.o			\
	fill_dipole_grid.o find_diss_mol.o find_igint_vor.o find_length_run.o	\
	find_thermos.o fit.o fitors.o fix_molecules.o fnb14_2.o fnb14.o		\
	fnbgrp2.o fnbgrp.o fndfmt.o fndgrp.o force_ele_nospline.o		\
	force_ele_nospline_sw.o force_ele_spline.o force_ele_spline_sw.o	\
	force_noele.o force_noele_sw.o fpbend.o fpbond_adjust_bonds.o fpbond.o	\
	fptors.o free_mem.o furippa.o gauss.o gen_abmd_kvect.o			\
	generate_rnd_rot.o generate_slv.o gen_stress_tag_bnd.o get_atres.o	\
	get_bspline_coeffs2.o get_bspline_coeffs.o getchr.o			\
	get_electric_field.o get_memory_iobuffer.o get_native.o			\
	get_pointer_protl.o get_prot_cnstr.o get_sgr.o GetSimpleCO.o		\
	get_spectra_vacf.o get_total_energy.o get_type_slv.o GetVariableUnit.o	\
	get_velocities.o get_virtual_energy.o gpfft.o grad_sum2.o		\
	grad_sum_dipole.o grad_sum.o hapacp.o hapdon.o iceil.o icopy.o igmap.o	\
	iindexx.o indexx.o inicmp.o initialize_cavities.o insert_solute.o	\
	int_corr_erf_spline.o interpol_dyna.o join.o kinetic.o latt_func.o	\
	L-BFGS-B.o lc_index.o lc_list2.o lc_list3.o lc_list.o link.o		\
	linmin_adjust_bonds.o linmin_total.o low_up.o lsqqtn.o mapnl_divide.o	\
	mappa.o matinv.o M_cavities.o memory.o M_free.o M_free_s.o		\
	M_get_length.o min_pack.o minv.o M_memory.o M_memory_s.o mnbrakf.o	\
	mnbrakf_total.o M_newlen.o mts_forces2.o mts_forces.o mts_forpp2.o	\
	mts_forpp.o mts_furier2.o mts_furier.o mts_intra_n0.o mts_intra_n1.o	\
	mtsmd_avg_inc.o mtsmd_dump_inc.o mtsmd.o mts_plot_fragm.o mts_plotp.o	\
	mts_test.o mts_umbrl.o M_voronoi.o near0.o normal.o openf.o		\
	open_file_dump.o orac.o P_adjust_decomp.o parse.o parse_numbers.o	\
	P_atoms_split_inter.o P_atoms_split_intra.o P_barrier.o pbend.o		\
	pbond.o P_broadcast_i.o P_Buffer_attach.o P_Buffer_detach.o		\
	P_check_decomp.o pckben.o pckfrg.o pckitor.o pcktor.o P_cm_split.o	\
	perc_bondh.o P_get_errmsg.o P_get_gofr.o P_get_intra.o P_get_iret.o	\
	P_get_neighbors.o P_get_prot.o phapp.o pitors.o plgndr.o		\
	P_load_balance.o plot_center.o plotc.o plotd.o plot_fragm.o		\
	plot_fragm_pdb.o pme_code.o P_npread.o P_nread.o pntres.o P_nwrite.o	\
	P_omit_intra.o P_open_input.o P_pme_assign_node.o pptors.o		\
	print_title_analysis.o prtacc.o prtat.o prtba.o prtbnd.o prtcn.o	\
	prtfrc.o prtite.o prtite_min.o prtit.o prtpt.o prtsq.o prttors.o	\
	P_setup_decompa.o P_setup_decomp.o P_split_constr.o P_split_intra.o	\
	P_split_scalar0.o P_split_scalar.o P_whoami.o ranf.o ranvel.o		\
	rattle_correc_co.o rattle_correc.o rattle_verlet_co.o rattle_verlet.o	\
	rdaprm.o rdatpg.o rdcmac.o rdrgd.o read_analysis.o readco.o		\
	read_confc_columns.o read_confc_rows.o read_inout.o read_input.o	\
	read_integrator.o read_parallel.o read_parameters.o read_potential.o	\
	read_properties.o readrs.o readrs_old.o read_run.o read_setup.o		\
	read_simulation.o read_solute.o read_solvent.o read_space_group.o	\
	readtp.o recpr.o redefine_slv.o redmss.o remove_mv.o repl_coord.o	\
	repl_tpg.o restim.o rigfit.o rmbend.o rmpbond.o rotb.o			\
	RotCoordColumn.o RotCoordRow.o run_minimize.o sbend.o sbond.o		\
	scalar_sum.o scalar_sum_normal.o scalar_sum_transp.o scale_charges.o	\
	sctab.o search_clsth.o serchd.o set_const_co.o set_hbonds_masks.o	\
	set_hhisto.o set_ss_array.o set_tempp.o set_tempt.o setuns.o		\
	setup_box.o setup_fake_map.o setup_skin_shell.o spec14.o start.o	\
	starting_verlet.o strblk.o test.o \
	time_correlation.o timer.o tors.o trans_center.o tr_inbox.o tr_wsc.o	\
	tr_wsc_new.o twofft.o update.o up_low.o uscrpl.o usort.o		\
	verify_input.o verify_variables.o verlet.o verlet_free_eta.o		\
	verlet_free.o virtual_extract.o windows.o wrenc.o write_bends.o		\
	write_bonds.o write_cavities_vor.o write_confc.o write_contacts.o	\
	write_diffusion.o write_electric_field.o write_fragm_dist.o		\
	write_gofrp.o write_gofrw.o write_gyr.o write_hbonds.o			\
	write_hbonds_vor.o write_hhisto.o write_pot_bond.o write_pot_nbond.o	\
	write_rms.o write_sk.o write_tors.o write_vacf_phi.o			\
	write_voronoi_header.o write_xrms_atm.o write_xrms.o xerror.o xfit.o	\
	zero0.o zero3.o zero3x3.o zeroa.o zero.o zero_gofr.o zeroi.o		\
	zero_voronoi.o


OBJ_NP = Tetrahedra.o xerror_mod.o Module_Neighbors.o RMS_Matrix_Mod.o	\
	Neighbors.o EUL_Mod.o ElecPotential_Mod.o Gauss_Overlap_Mod.o		\
	PDB_utils.o PDB.o Neighbors_Mod.o PBC_Mod.o DENSITY_Mod.o		\
	Class_Connect.o Class_Bending.o Class_Gauss_Param.o Module_Fourier.o	\
	Class_Gcharges.o POLAR_Mod.o Module_Extra_Forces.o reset_nbond.o	\
	Poisson_Boltzmann.o Module_Stress.o Module_Thole.o CENTER_SOL_Mod.o	\
	RIGID_Mod.o INPUT_Mod.o PDBs_Mod.o RMS_Subtract_Mod.o GROUPS_Mod.o	\
	GEOM_groups_Mod.o WSC_Mod.o HEATING_Mod.o IOBUFFER_Mod.o		\
	PDB_FRAGM_Mod.o add_energies.o add_int1.o add_int2.o add_real1.o	\
	add_real2.o add_char1.o add_concta.o add_mapnl.o add_solvent_tpg.o	\
	add_solvent_coord.o add1h.o add2h.o add3h.o add_str.o addh.o addh2o.o	\
	addlnk.o addtor.o addtpg.o adjust_bonds.o anneal.o appatm.o appbkn.o	\
	appbou.o appct.o appgrp.o asng_xrmsw.o asstpe.o atoms.o atthd.o		\
	bcnstp.o bend.o blkdta.o bond.o brent.o brot.o CalcDistConf.o		\
	calc_avg_str.o calc_avg2_str.o calc_avg_xrms.o calc_gofr.o calc_sk.o	\
	calc_sk_eet.o calc_xslt.o calc_inst_fit.o calc_RotMat.o calc_cofm.o	\
	calc_lda.o calc_lda_avg.o calc_lda_flu.o calc_lda_hyd.o calc_hyd.o	\
	calc_xrms.o CalcXrmsSecStruct.o caldr.o cc.o change_coord_ex.o		\
	change_coord_inner.o change_virtual_potential.o virtual_extract.o	\
	calc_lda_rest.o change_frame.o change_origin.o change_tors.o		\
	change_step_sim.o check_length_sim.o check_read_columns.o		\
	check_topology.o check_zero_coord.o chk14.o chkewl.o chlab.o chmass.o	\
	chnpr.o close_files.o clpcm.o clt14.o collision.o comd0.o		\
	comp_abmd_fdiss.o comp_abmd_fcryst.o comp_ang_vel.o			\
	comp_cavities_vor.o comp_abmd_tors.o comp_cell_center.o comp_concth.o	\
	comp_contacts.o comp_dip.o comp_dynamic_mat.o comp_fcm.o comp_forcep.o	\
	comp_gmgp.o comp_molmass.o comp_native.o comp_neigh_vor.o		\
	comp_noneq_free.o comp_nmol.o comp_stress_conf.o comp_stress_kinetic.o	\
	comp_stress_bnd.o comp_thermos_energy.o comp_rmsq.o			\
	comp_thermos_forces.o comp_thermos_forces_cm.o comp_vcm.o		\
	comp_vel_labframe.o correc.o correc_etap.o Correct_Ewald.o		\
	coordinate_spline.o coordinate_spline_init.o copy_protl.o correc3x3.o	\
	correc_gmgp.o correc_matr.o correc_exp_scale.o correc_scale.o		\
	correc_scale_2.o correc_scale_co.o cov_thermos.o covbd.o covchg.o	\
	covfod.o covnb.o covpr.o covuen.o crdrtf.o cself.o cself_dipole.o	\
	cself_phi.o csort.o dcopy3.o delrfg.o dist_ijk.o dlineq.o DCopy_V3.o	\
	dspat.o dumptp.o duni.o dupl_char1.o dupl_concta.o dupl_int1.o		\
	dupl_int2.o dupl_mapnl.o dupl_real1.o dupl_real2.o dupl_reali2.o	\
	erf_corr_cutoff.o erfc_spline.o erfc_spline_init.o equi_bond.o		\
	equi_bend.o equi_tors.o equi_impro.o f1dim_adjust_bonds.o		\
	f1dim_total.o ferrf.o fit.o ferrf_tag.o fftpack.o fft3-public.o		\
	fft_default.o fft_pme_init.o fill_charge_grid.o fill_dipole_grid.o	\
	find_diss_mol.o find_igint_vor.o find_length_run.o find_thermos.o	\
	fitors.o fix_molecules.o fnb14.o fnb14_2.o fnbgrp.o fnbgrp2.o fndfmt.o	\
	fndgrp.o fpbend.o fpbond.o fpbond_adjust_bonds.o fptors.o gauss.o	\
	gen_abmd_kvect.o gen_stress_tag_bnd.o generate_slv.o			\
	generate_rnd_rot.o get_atres.o get_bspline_coeffs.o			\
	get_bspline_coeffs2.o GetKineticTest.o GetHardSphere.o			\
	get_virtual_energy.o Get_vp_with_vc.o Get_vp_without_vc.o		\
	get_pointer_protl.o get_native.o get_prot_cnstr.o get_sgr.o		\
	get_spectra_vacf.o get_type_slv.o iceil.o get_velocities.o getchr.o	\
	GetSimpleCO.o grad_sum.o grad_sum2.o grad_sum_dipole.o hapacp.o		\
	hapdon.o igmap.o icopy.o indexx.o iindexx.o inicmp.o			\
	initialize_cavities.o Init_Fixed_Molecule.o int_corr_erf_spline.o	\
	interpol_dyna.o chebft_dyn.o chder_dyn.o chebev_dyn.o join.o kinetic.o	\
	L-BFGS-B.o latt_func.o lc_index.o link.o LinkedCellIndex.o		\
	linmin_adjust_bonds.o LoadBspModuli.o low_up.o lsqqtn.o mappa.o		\
	mapnl_divide.o matinv.o min_pack.o minv.o mnbrakf.o Make_Pol_Models.o	\
	mts_plot_fragm.o mts_plotp.o mts_test.o mts_umbrl.o near0.o normal.o	\
	Open_Input.o parse.o parse_numbers.o pbend.o pbond.o pckben.o pckfrg.o	\
	pckitor.o pcktor.o phapp.o pitors.o plgndr.o plot_center.o plotc.o	\
	plot_fragm.o plotd.o pntres.o pptors.o prtacc.o prtat.o prtba.o		\
	prtbnd.o prtcn.o prtfrc.o prtit.o prtite.o prtite_min.o			\
	print_title_analysis.o prtpt.o prtsq.o prttors.o ranf.o ranvel.o	\
	rattle_verlet.o rattle_verlet_co.o rattle_correc.o rattle_correc_co.o	\
	rdaprm.o rdatpg.o rdcmac.o rdrgd.o read_analysis.o read_inout.o		\
	read_integrator.o read_parameters.o read_potential.o read_properties.o	\
	read_run.o read_setup.o read_simulation.o read_solute.o read_solvent.o	\
	read_space_group.o readco.o readtp.o redefine_slv.o redmss.o		\
	remove_mv.o repl_coord.o repl_tpg.o ResetCofM.o restim.o rigfit.o	\
	rmbend.o rmpbond.o rotb.o RotCoordColumn.o RotCoordRow.o sbend.o	\
	sbond.o ScalarSum.o scalar_sum_normal.o scalar_sum_transp.o		\
	scalar_sum.o scale_charges.o sctab.o search_clsth.o serchd.o		\
	set_const_co.o set_hbonds_masks.o set_hhisto.o set_ss_array.o		\
	set_tempp.o set_tempt.o setuns.o setup_fake_map.o Setup_fake_decompa.o	\
	setup_skin_shell.o spec14.o starting_verlet.o strblk.o			\
	time_correlation.o tors.o tr_inbox.o tr_wsc.o trans_center.o twofft.o	\
	up_low.o update.o uscrpl.o usort.o verify_variables.o verlet.o		\
	verlet_free.o verlet_free_eta.o windows.o wrenc.o write_contacts.o	\
	write_fragm_dist.o write_pot_bond.o write_pot_nbond.o write_gyr.o	\
	write_gofrp.o write_gofrw.o write_hbonds.o write_hbonds_vor.o		\
	write_hhisto.o write_bends.o write_bonds.o write_tors.o			\
	write_vacf_phi.o write_diffusion.o write_voronoi_header.o write_rms.o	\
	write_sk.o write_xrms.o write_xrms_atm.o xfit.o zero0.o zero.o zero3.o	\
	zero3x3.o zero_gofr.o zero_voronoi.o zeroa.o zeroi.o			\
	write_mol_dipole.o write_electric_field.o P_cm_split.o P_get_prot.o	\
	P_get_intra.o Polarization_Forces.o Do_Electric_Fields.o		\
	ORAC_Finalize.o fft_pme.o VORONOI_Mod.o HYDRATION_Mod.o			\
	P_Change_Decomposition.o P_Comm_Intra.o do_rfft3d.o setup_box.o		\
	pme_code.o analyse_voronoi.o bldbox.o start.o furippa.o			\
	insert_solute.o mts_forpp.o comp_voronoi.o comp_abmd_fnative.o		\
	write_cavities_vor.o timer.o verify_input.o open_file_dump.o orac.o	\
	read_input.o free_mem.o M_free.o M_free_s.o run_minimize.o		\
	get_electric_field.o get_memory_iobuffer.o get_total_energy.o		\
	f1dim_der.o comp_total_energy.o linmin_total.o mnbrakf_total.o		\
	brent_total.o mts_intra_n0.o mts_intra_n1.o mts_furier.o mts_furier2.o	\
	do_pmesh_kspace.o do_pmesh_kspace_dipole.o write_confc.o		\
	read_confc_rows.o read_confc_columns.o M_get_length.o M_newlen.o	\
	M_memory.o memory.o M_memory_s.o fft_pme_dipole.o openf.o fft_setup.o	\
	fft_forward.o fft_back.o correc_stress.o correc_stress_n0.o lc_list.o	\
	lc_list2.o lc_list3.o mts_forces.o mts_forces2.o mts_forpp2.o		\
	read_parallel.o P_setup_decomp.o P_setup_decompa.o			\
	P_atoms_split_intra.o P_split_intra.o P_split_constr.o			\
	P_atoms_split_inter.o xerror.o HYDYNAMICS_Mod.o FIELD_Mod.o Pme_init.o	\
	REDUCE_Mod.o


OBJ_PP = P_omit_intra.o P_broadcast_i.o P_check_decomp.o P_expand_r8.o	\
	P_expand_r8x3.o P_expand_i.o P_fold_r8.o P_fold_r8x3.o P_fold_i.o	\
	P_get_errmsg.o P_get_gofr.o P_load_balance.o P_merge_i.o P_merge_r8.o	\
	P_merge_vecr8.o P_merge_veci.o P_split_scalar0.o P_split_scalar.o	\
	P_get_iret.o P_adjust_decomp.o P_Buffer_attach.o P_Buffer_detach.o	\
	P_get_neighbors.o P_Write_neighbors.o P_whoami.o P_open_input.o		\
	P_nwrite.o P_nread.o P_npread.o P_barrier.o P_pme_assign_node.o



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


ifeq ($(TARGET0), Linux)
OBJREC = recpr_linux.o perc_bondh_linux.o
OBJREC1 = recpr.o perc_bondh.o
else
OBJREC = recpr.o perc_bondh.o
OBJREC1 = recpr_linux.o perc_bondh_linux.o
endif

OBJ_LOWOPT = eigrs.o mtsmd.o drive_analysis.o 

OBJ_C = vstart.c

INCLUDE = cpropar.inc parst.inc orac.inc fourier.inc unit.inc \
	resid.inc vacf.inc pbc.inc analysis.inc mts_spectra.inc\
	mts_vacf.inc pme.inc parallel.inc parameters.inc\
	lc_list.inc voronoi.inc iobuffer.inc mpi_size.inc

INCPRE = $(patsubst %.inc, %.h, $(INCLUDE))

#==========================================================================
# Uncomment the folling if an old restart has to be used
# 
#DUMP = dumprs_old.o readrs_old.o
#==========================================================================

DUMP = dumprs.o readrs.o readrs_old.o

SRC = $(patsubst %.o, %.f90, $(OBJ_F90)) $(patsubst %.o, %.f, $(OBJ_F)) \
	$(patsubst %.o, %.c, $(OBJ_C))

#==========================================================================
# List of sources 
#==========================================================================


LIST_FS  = $(SRC)
LIST_CS  = vstart.c 
LIST_EXT = Makefile Make.* config.h mtsmd_avg_inc.f mtsmd_dump_inc.f $(INCLUDE)


#==========================================================================
# Make Default 
#==========================================================================

all: makedepf90 $(INCPRE)
	${MAKE} $(TARGET0)


makedepf90:
	(cd $(PWD)/makedepf90-2.8.8/; ./configure; ${MAKE} ; cp makedepf90 ../)



#==========================================================================
# Make include
#==========================================================================

iobuffer.h: iobuffer.inc
	$(CPP) -D$(TARGET0) iobuffer.inc | awk '!/^ *$$/ && !/^#/' > temp
	./python/include_f90.py temp
	file=iobuffer;\
	test -f $$file.h && ( cmp -s $$file.h temp || \
	cp temp $$file.h) || cp temp $$file.h ;\
	rm -f temp temp2 >/dev/null

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
# Decide what to make
#==========================================================================


PGI :
	${MAKE} obj_PGI "FFLAGS = -Mpreprocess -O2 "
	${MAKE} lowopt_PGI "FFLAGS = -Mpreprocess -O0"
	${MAKE} recur_PGI "FFLAGS = -Mpreprocess  "
	${MAKE} orac_PGI

IFC :
	${MAKE} obj_IFC "FFLAGS = -fpp -O3" 
	${MAKE} lowopt_IFC "FFLAGS = -fpp -O0"
	${MAKE} recur_IFC "FFLAGS = -fpp  "
	${MAKE} orac_IFC 

IFC_DBG :
	${MAKE} obj_IFC "FFLAGS = -fpp -w95 -g -traceback -check all" 
	${MAKE} lowopt_IFC "FFLAGS = -fpp -w95 -g -traceback -check all"
	${MAKE} recur_IFC "FFLAGS = -fpp -w95 -g -traceback -check all"
	${MAKE} orac_IFC 

XLF :
	${MAKE} obj_XLF "FFLAGS = -O5 -qnosave -qhot -qipa" 
	${MAKE} lowopt_XLF "FFLAGS = -O2 -qsave -qipa "
	${MAKE} recur_XLF "FFLAGS = -O4 -qnosave -qipa"
	${MAKE} orac_XLF 

GFORT :
	${MAKE} obj_GFT "FFLAGS = -x f95-cpp-input  -O3" 
	${MAKE} lowopt_GFT "FFLAGS = -x f95-cpp-input  -O0"
	${MAKE} recur_GFT "FFLAGS = -x f95-cpp-input   "
	${MAKE} orac_GFT

Do_Nothing:


#========================================================================

DISTRIBUTION : 
	@rm -r -f $(CURDIR)/tmp_dir
	mkdir $(CURDIR)/tmp_dir ;\
	mkdir $(CURDIR)/tmp_dir/orac2k ;\
	mkdir $(CURDIR)/tmp_dir/orac2k/source ;\
	mkdir $(CURDIR)/tmp_dir/orac2k/source/python ;\
	mkdir $(CURDIR)/tmp_dir/orac2k/source/makedepf90-2.8.8 ;\
	for j in $(LIST_FS);\
	do cp -f $(CURDIR)/$$j $(CURDIR)/$(DIST_DIR)/;done;\
	for j in $(LIST_CS);\
	do cp -f $(CURDIR)/$$j $(CURDIR)/$(DIST_DIR)/;done;\
	for j in $(LIST_EXT);\
	do cp -f $(CURDIR)/$$j $(CURDIR)/$(DIST_DIR)/;done;\
	cp python/include_f90.py $(CURDIR)/$(DIST_DIR)/python/;
	cp -r makedepf90-2.8.8 $(CURDIR)/$(DIST_DIR)/;
	(cd $(CURDIR)/$(DIST_DIR)/makedepf90-2.8.8 ; make clean);

#==========================================================================
# Compile and load ORAC on OSF   
#==========================================================================


obj_PGI : $(OBJ) $(DUMP) $(OBJTIM)

recur_PGI: $(OBJREC)

lowopt_PGI: $(OBJ_LOWOPT)

orac_PGI : $(OBJ) $(OBJREC) $(DUMP) $(OBJ_LOWOPT) $(OBJTIM)
	$(LD) -Mvect -Mipa=fast -fastsse  -o orac_PGI $(OBJ) $(OBJREC) $(DUMP) \
	$(OBJ_LOWOPT) $(OBJTIM)\
	-L$(FFTW_LIB) -L$(OTHER_LIBS) $(LLFLAGS) -lm


#==========================================================================
# Compile and load ORAC on INTEL Compiler   
#==========================================================================


obj_IFC : $(OBJ) $(DUMP) $(OBJTIM)

recur_IFC: $(OBJREC)

lowopt_IFC: $(OBJ_LOWOPT)


orac_IFC : $(OBJ) $(OBJREC) $(DUMP) $(OBJ_LOWOPT) $(OBJTIM)
	$(LD) -o orac_IFC  $(OBJTIM) $(OBJ) $(OBJREC) $(DUMP) $(OBJ_LOWOPT)  \
	-L$(FFTW_LIB) -L$(OTHER_LIBS) $(LLFLAGS) 


#==========================================================================
# Compile and load ORAC on gfortran
#==========================================================================


obj_GFT : $(OBJ) $(DUMP) $(OBJTIM)

recur_GFT: $(OBJREC)

lowopt_GFT: $(OBJ_LOWOPT)


orac_GFT : $(OBJ) $(OBJREC) $(DUMP) $(OBJ_LOWOPT) $(OBJTIM)
	$(LD) -o orac_GFT  $(OBJTIM) $(OBJ) $(OBJREC) $(DUMP) $(OBJ_LOWOPT)  \
	-L$(FFTW_LIB) -L$(OTHER_LIBS) $(LLFLAGS) 


#==========================================================================
# Compile and load ORAC on XLF compiler
#==========================================================================

obj_XLF : $(OBJ) $(DUMP) $(OBJTIM)

recur_XLF: $(OBJREC)

lowopt_XLF: $(OBJ_LOWOPT)

orac_XLF: $(OBJ) $(OBJREC) $(DUMP) $(OBJ_LOWOPT) $(OBJTIM)
	$(LD) -O5 -qnosave -qhot -qipa -o orac_XLF $(OBJ) $(OBJREC)\
	$(DUMP) $(OBJ_LOWOPT) $(OBJTIM)\
	-L$(FFTW_LIB) -L$(OTHER_LIBS) $(LLFLAGS)

#==========================================================================
# Clean Objects
#==========================================================================

clean :
	rm -f *.o *.mod orac_*

#==========================================================================
# Clean Objects and executable
#==========================================================================

extra_clean :
	rm -f *.o *.mod orac_* makedepf90
	(cd $(CURDIR)/makedepf90-2.8.8 ; make clean)

#==========================================================================
# Get all *.f files,  checkin all files, tar 
#==========================================================================

tar: DISTRIBUTION
	(cd $(CURDIR)/tmp_dir; tar zcvf orac.tar.gz orac2k)

clean_distribution:
	rm -r -f tmp_dir;

#==========================================================================
# Dependencies 
#==========================================================================

# Dependencies
Dependencies.h: $(shell ./makedepf90 $(SRC) > Dependencies.h)
include Dependencies.h
##$(shell rm Dependencies.h)
