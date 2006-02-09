#define _REAL_ real*4
      SUBROUTINE drive_analysis(mapnl,xp0,yp0,zp0,xpg,ypg,zpg,eta,xpcm
     &     ,ypcm,zpcm,node,nodex,nodey,nodez,ictxt,npy,npz,nprocs,ncube)

************************************************************************
*   Time-stamp: <2006-02-09 14:06:55 marchi>                             *
*                                                                      *
*     drive_analysis analize a trajectory file written by mtsmd        *
*     In addition to that file also a binary topology file must        *
*     be provided                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed Jul  2 1997 -                                     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*  drive_analysis externals:   	       	       	       	       	       *
*       analyse_voronoi appbou asng_xrmsw calc_avg2_str calc_avg_str   *
*       calc_avg_xrms calc_gofr calc_xrms change_frame change_step_sim *
*       check_length_sim check_read_columns check_topology 	       *
*       comp_cell_center comp_dip comp_molmass comp_neigh_vor 	       *
*       comp_rmsq comp_voronoi coordinate_spline coordinate_spline_init*
*       daxpy dcffti dcopy dscal find_igint_vor find_length_run fndgrp *
*       get_spectra_vacf get_type_slv get_velocities inicmp matinv     *
*       mts_plot_fragm plotc plotd print_title_analysis		       *
*       prtat prtba prtbnd prtit prtpt prtsq			       *
*       read_confc_columns read_confc_rows set_hbonds_masks timer      *
*       time_correlation tr_inbox update write_bends write_bonds       *
*       write_diffusion write_fragm_dist write_gofrp write_gofrw       *
*       write_gyr write_hbonds write_rms write_tors write_vacf_phi     *
*       write_voronoi_header write_xrms write_xrms_atm xerror	       *
*       zero0 zeroa zero_gofr zero_voronoi			       *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================= DECLARATIONS ==================================

      USE HYDRATION_Mod, ONLY: hydration,HYD_n_neighbors=>n_neighbors,
     &     HYD_Initialize_P=>Initialize_P,
     &     HYD_Initialize_Array=>Initialize_Array,
     &     HYD_Compute=>Compute,
     &     HYD_Compute_Neighbors=>Compute_Neighbors,
     &     HYD_n_write=>n_write,HYD_write_it=>write_it

      USE RMS_Matrix_Mod, ONLY: RMS_Initialize=>Initialize,
     &     RMS_Compute=>Compute,RMS_Write_it=>Write_it, RMS_Compute_avg
     &     =>Compute_avg,RMS_Write_it_avg=>Write_it_avg, rms_matrix
     &     ,rms_matrix_avg, rms_matrix_plot,RMS_Rotate=>Rotate, RMS_plot
     &     =>plot
      USE EUL_Mod, ONLY: EUL_Initialize=>Initialize, eul_angles,
     &     EUL_Compute=>Compute
      USE DENSITY_Mod, ONLY: DEN_write=>n_write,Density_Calc
     &     ,DEN_Initialize=>Initialize,DEN_Initialize_Ar=>Initialize_Ar
     &     ,DEN_compute=>Compute,DEN_write_it=>write_it
      USE CENTER_SOL_Mod, ONLY: CEN_Center=>Center_Object, CEN_Compute
     &     =>Compute,CEN_Initialize=>Initialize
      USE PDBs_Mod, ONLY: PDB_pdbs=>PDBs, PDB_Initialize=>Initialize
     &     ,PDB_compute=>Compute,PDB_write=>Write_it,PDB_nwrite=>n_write
     &     ,PDB_ncompute=>n_compute
      USE RMS_Subtract_Mod, ONLY: SUB_Initialize=>Initialize
     &     ,SUB_Subtract=>RMS_Subtract,SUB_Compute=>Compute,SUB_write
     &     =>n_write,SUB_Write_it=>Write_it,SUB_Rotate=>Rotate,SUB_Start
     &     =>Start
      USE GROUPS_Mod, ONLY: GR_groups=>groups,GR_Init=>Init
      USE GEOM_groups_Mod, ONLY: GE_Groups=>Geom_groups, GE_init
     &     =>Init, GE_compute=>Compute,GE_Write=>n_write,GE_output
     &     =>Write_it 
      
      IMPLICIT none
      
      include 'parst.h'

      INTEGER ma,mb,mh,ms,mf
      PARAMETER (ma=tsites,mb=tsitep,mf=mb*4)

*----------------------- ARGUMENTS -------------------------------------

      INTEGER node,ncube,nprocs,mapnl(*)
      INTEGER nodex,nodey,nodez,ictxt,npy,npz
      REAL*8  xp0(*),yp0(*),zp0(*),xpg(*),ypg(*),zpg(*),eta(*),xpcm(*)
     &     ,ypcm(*),zpcm(*)

*-------------------- VARIABLES IN COMMONS -----------------------------

      INCLUDE 'cpropar.h'
      include 'fourier.h'
      include 'pme.h'
      INCLUDE 'lc_list.h'
      INCLUDE 'iobuffer.h'
      INCLUDE 'voronoi.h'
      INCLUDE 'analysis.h'
      INCLUDE 'unit.h'
      
      
*-------------------- DEFINITION OF AN EXTERNAL FUNCTION ---------------

      EXTERNAL  near0
      LOGICAL near0

*-------------------- LOCAL VARIABLES ----------------------------------

      CHARACTER*80 errmsg

      INTEGER    l2
      INTEGER numatoms,descQ(12)
      INTEGER    nato_slt,nbun_slt,iter_avg,iter_avg2
      REAL*8  fact

      INTEGER mapnl0(m1),ingrpp0,ingrp0(2,3*m11),ingrp0_x(3*m11)

      INTEGER ntot_fragm,fragm_1,fragm_2
      INTEGER i,j,nstep,iret,istep,ia,ib
      LOGICAL end,err,Regular

*----------- LOCAL WORK ARRAYS FOR THE RUN -----------------------------

      INTEGER, DIMENSION(:), POINTER ::  indxi,indxj,indxk,nnlpp0
      INTEGER nnlppf(sitslu+1)

      INTEGER mapdn(mf),nmapdn(m1),hhisto_list(3,hhisto_maxbin)
     &     ,hhisto_dim,hhisto_count,indxyz,ind_a,nind,npp

      INTEGER nstart_h,nend_h,nlocal_h,nstart_ah,nend_ah,nlocal_ah
     &     ,nstart_uh,nend_uh,nlocal_uh,nend_lda

      REAL*8  qt0(4)
      INTEGER worka(m1,m10)
      LOGICAL a_mask(m1),d_mask(m1)
      REAL*8  work(mspline)      

*--- DYNAM is a scratch common block: here to save storage; not passed
*    to any of the external

      
      INTEGER type_slv(slvatm)
      INTEGER(4), POINTER :: krdf(:)
#ifdef PARALLEL 
#define _KRDF_   krdfa
      INTEGER(4), POINTER ::  krdfa(:)
#else
#define _KRDF_   krdf
#endif
      INTEGER ntype_slv,offset_slv,nbetab_slv,nstep_start_o
      PARAMETER(nbetab_slv=90)
      CHARACTER*1 betab_slv(nbetab_slv),P_shell

      INTEGER list(3,m1*5),hlist(3,m1*5),rlist(3,m1*5),nlist


      INTEGER itype_slv(nbetab_slv),npp_u,mpp_a,M_get_length,npp_m
      _REAL_  xau(mb),yau(mb),zau(mb)

      REAL*8 ffstep,fstep,fstep1,fstep2,fstep_o,xt_cm,yt_cm,zt_cm
      REAL*8  tmass(npm),vfcp,tfcp,gcpu
      REAL(8), DIMENSION (:), POINTER :: xp_avg,yp_avg,zp_avg,xp_avg2
     &     ,yp_avg2,zp_avg2
      REAL(8), DIMENSION (:), POINTER :: errca,errhe,errbc,erral,drpca
     &     ,drpbc,drphe,drpal,wca,whe,wbc

      REAL(8), POINTER ::  xpo(:),ypo(:),zpo(:)
      REAL(8), POINTER ::  xpa(:),ypa(:),zpa(:),xpga(:),ypga(:),zpga(:)
     &     ,xpcma(:),ypcma(:),zpcma(:)

      REAL*8 dips(3,2),vol_gofr,sum_volume,elapse,aux

      REAL*8  sofk_data(maxsk,2),nsofk_data(maxsk,2),Dist
      REAL*8  ucns_p,ucos_p,virs_p,virsp_p,ucnp_p,ucop_p,ucnsp_p,ucosp_p
     &     ,fpx_p,fpy_p,fpz_p,stressd_p,cpu_h,gr,gra,uconf,enout,fpx(m1)
     &     ,fpy(m1),fpz(m1),ffwork(2),InstRotMat(3,3),Ixpcc,Iypcc,Izpcc
      INTEGER offset,nnstep,start_time,end_time,length_run,length_tot
     &     ,length_fft,iatom,iatom0,niatom,naux,i_old,i_start
     &     ,i_end,typei,nato1,nato2,jmax_cav,mqq,mma,ncpu_h,mmb
      INTEGER nat_listp,nat_cntactp,nat_list(2,nores),nat_cntact(nores)

      COMMON /dynam/ w1,w2,mapdn,nmapdn,worka,nnlpp0
     &     ,a_mask,d_mask,xpc,ypc,zpc,vxpc,vypc,vzpc,vxpd,vypd,vzpd,list
     &     ,hlist,rlist,wsave1,spline_x,spline_y

      COMMON /rag2/ vacf_data,rms_disp,tot_rms_disp,xpb,ypb,zpb,xpcc
     &     ,ypcc,zpcc,xpo,ypo,zpo,RotMat,xau,yau,zau,xpa,ypa,zpa,xpga
     &     ,ypga,zpga,xpcma,ypcma,zpcma,wca,whe,wbc,errca,errhe,errbc
     &     ,erral,drpca,drpbc,drphe,drpal,xp_avg,yp_avg,zp_avg,xp_avg2
     &     ,yp_avg2,zp_avg2,tmass

*==================== EXECUTABLE STATEMENTS ============================
c$$$====================================================================
c$$$-- Allocate something 
c$$$====================================================================

      ALLOCATE(xpo(ntap),ypo(ntap),zpo(ntap))
      ALLOCATE(xpa(ntap),ypa(ntap),zpa(ntap))
      ALLOCATE(xpga(ngrp),ypga(ngrp),zpga(ngrp))
      ALLOCATE(xpcma(nprot),ypcma(nprot),zpcma(nprot))

c$$$====================================================================
c$$$-- 
c$$$====================================================================

      npp_u=mpp
      indxyz=indmax
      end=.FALSE.

#ifdef PARALLEL      
      CALL P_setup_decompa(nstart_h,nend_h,nlocal_h,nstart_ah,nend_ah
     &     ,nlocal_ah,nstart_uh,nend_uh,nlocal_uh,ss_index,nbun,ngrp
     &     ,grppt,mres,node,nprocs,ncube,nbyte)
#else
      CALL Setup_fake_decompa(nstart_h,nend_h,nlocal_h,nstart_ah,nend_ah
     &     ,nlocal_ah,nstart_uh,nend_uh,nlocal_uh,ntap,ngrp,nbun)
#endif


      IF(linked_cell) THEN
         indxyz=(2*(ncx-1)+1)*(2*(ncy-1)+1)*(2*(ncz-1)+1)
         ALLOCATE(indxi(indxyz),indxj(indxyz),indxk(indxyz))
      END IF

      IF(nprocs .EQ. 1) THEN
         npp_u=mpp
      ELSE
         npp_u=(mpp/nprocs)*2
      END IF
      ALLOCATE(nnlpp0(npp_u))
      IF(nprocs .EQ. 1) THEN
         mmb=m1
      ELSE
         mmb=m1/(nprocs-1)
      END IF
      IF(hydration) THEN
         CALL HYD_Initialize_P(node,nprocs,ngrp,nbun)
         CALL HYD_Initialize_Array(grppt,mres,resg)
      END IF
      IF(Density_Calc) THEN
         CALL DEN_Initialize(prsymb,nres(1,2),mass,ntap)
         CALL DEN_Initialize_Ar
      END IF
      IF(CEN_Center) THEN
         CALL CEN_Initialize(prsymb,beta,nres(1,2),ss_index,ntap)
      END IF
      IF(PDB_pdbs) THEN
         CALL PDB_Initialize(chrge,prsymb,beta,nres(1,1),nres(1,2),ntap)
      END IF

*===  Check if the dimension of the work array are sufficient 

      IF(mb.LT.ntap) THEN
         errmsg=' While running Analysis: PARAMETER MB dimensions the'
     &/ /' work arrays in sufficiently. Abort. '
         CALL xerror(errmsg,80,1,2)
      END IF
    
*=======================================================================
*----- Initialize some stuff -------------------------------------------
*=======================================================================

      rcuth=rspoff
      rtolh=0.0D0
      rneih=rspcut

*===  set few variable to zero

      iter_avg=0
      iter_avg2=0

      sum_volume=0.0D0
      iret=0
      errmsg=' '
      nato_slt=ntap-nmol*nato_slv

      if(prot_lda) then
         nend_lda = ngrp-nmol-nlda_mol*15
      end if
*=======================================================================
*---- Divide mapnl -----------------------------------------------------
*=======================================================================
 
      IF(efield) THEN
 
         call zero(ecc6,ecc12,m5,m5)
         mapnl0(1:m1)=0
         CALL igmap(ngrp,grppt,ingrpp0,ingrp0,3*m11,mapnl0,errmsg,iret)
         IF(iret.EQ.1) CALL xerror(errmsg,80,1,2)
         CALL mapnl_divide(node,nstart_h,nend_h,grppt,mapnl0)

         if(polar) then

#if defined PARALLEL
            CALL P_split_intra(node,nprocs,ncube,nstart_ah,nend_ah,worka
     &           ,iret,errmsg)
            IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)

            CALL P_get_intra(ingrp0,2,ingrpp0,ingrp0_x,nstart_ah,nend_ah
     &           ,iret)
            CALL P_get_iret(iret,node,nprocs,ncube,nbyte)
#endif
            IF(iret .EQ. 1) THEN
               errmsg=
     &  ' Shared ingrp0 bond larger than expected. Change n_h'
     & / /' in MTSMD '
               RETURN
            END IF

#if defined PARALLEL
            CALL P_omit_intra(node,nprocs,ncube,nbyte,ingrp0_x,worka
     &           ,iret)
#endif
         end if

      END IF

      CALL mapnl_divide(node,nstart_h,nend_h,grppt,mapnl)

*========================================================================
*==== Calls init routine for conventional kspace Ewald or PME -----------
*========================================================================

*=======================================================================
*-------- Find out the first and last group of each protein ------------
*=======================================================================

      CALL fndgrp(nprot,ngrp,protl,grppt,atomg,protg,groupp,atomp,npm
     &     ,iret,errmsg)
               
      CALL comp_molmass(nprot,protl,mass,tmass)

      IF(native .OR. check_native) THEN
         CALL get_atres(atres,nres(1,1),nato_slt,nbun_slt)
      END IF
      IF(check_native) THEN
         CALL get_native(knative_tpl,nat_listp,nat_list,nores,atres,beta
     &        ,atres_map1,atres_map2,iret,errmsg)
         IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
      END IF

      IF(voronoi) THEN
         IF(nprocs .EQ. 1) THEN
            mma=ntap
         ELSE
            mma=ntap/(nprocs-1)+10
         END IF
         ALLOCATE(nnlpp_vor(pnnlpp_vor,mma),ig_nnl(pig_nnl,mma))
         ALLOCATE(area_vor(pnnlpp_vor,mma), volume_vor(ntap))
      END IF
      IF(cavities) THEN
         IF(nprocs .EQ. 1) THEN
            mma=ntap
         ELSE
            mma=ntap/(nprocs-1)+10
         END IF
         mqq=maxcav_atom*mma
         ALLOCATE(cavity_n(2,mma),cavity_h(maxcav_bin,maxcav_nres)
     &        ,cavity_r(mqq))
      END IF
      IF(voronoi) THEN
         IF(node .EQ. 0) CALL write_voronoi_header(kvoronoi)
      END IF
      IF(cavities) THEN
         CALL initialize_cavities(ss_index,ntap,nres(1,1),bin_size_cav
     &        ,rmax_size_cav,cavity_h,maxcav_bin,maxcav_nres,jmax_cav
     &        ,iret,errmsg)
         IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
      END IF
      IF(hbonds_tot .OR. hbonds_res .OR. hhisto) THEN
         hrcut_update=rcut_hb+hrcut_update
      END IF

      IF(hhisto) THEN
         CALL set_hhisto(hhisto_list,hhisto_bin,hhisto_dim,hhisto_maxbin
     &        ,iret,errmsg)
         hhisto_count=0
         IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
      END IF

      IF(hbonds_vor) THEN
         CALL set_hbonds_masks(ntap,lacc,llacc,ldon,lldon,a_mask,d_mask)
      END IF
      IF(hbonds_tot .OR. hbonds_res .OR. hhisto) THEN
         acut_hb=DCOS(acut_hb*pi/180.0D0)
         a2cut_hb=DCOS(a2cut_hb*pi/180.0D0)
      END IF

      ALLOCATE(wca(ntap),whe(ntap),wbc(ntap))
      wca=0.0D0
      whe=0.0D0
      wbc=0.0D0
      CALL asng_xrmsw(ss_point,m1+1,wca,whe,wbc,beta,mback,nbone)
      IF(rms_matrix) THEN
         CALL RMS_Initialize(wca,xpt0,ypt0,zpt0)
      END IF
      IF(EUL_angles) THEN
         CALL EUL_Initialize(wca,xpt0,ypt0,zpt0)
      END IF
      IF(SUB_Subtract) THEN
         CALL get_atres(atres,nres(1,1),nato_slt,nbun_slt)
         CALL SUB_Initialize(atres,nbun_slt,nres(1,2),prsymb,beta
     &        ,chrge,DSQRT(unitc),nato_slt)
      END IF

      IF(GR_Groups) THEN
         CALL GR_Init(ntap,nbun,mres)
      END IF
      IF(GE_Groups) THEN
         CALL GE_Init(ntap)
      END IF

      IF(anxrms) THEN
         ALLOCATE(errca(nprot),errhe(nprot),errbc(nprot),erral(nprot)
     &        ,drpca(ntap),drpbc(ntap),drphe(ntap),drpal(ntap))
         errca=0.0D0
         errhe=0.0D0
         errbc=0.0D0
         erral=0.0D0
         drpca=0.0D0
         drpbc=0.0D0
         drphe=0.0D0
         drpal=0.0D0
      END IF
      IF(anxrms_cell) THEN
         anprot=.TRUE.
         annpro=annpro+1
         anpoint(1,annpro)=1
         anpoint(2,annpro)=ntap
      END IF

      IF(gofr) THEN
         CALL get_type_slv(nato_slt,nato_slv,beta,betab_slv
     &        ,nbetab_slv,ntype_slv,type_slv,itype_slv,offset_slv
     &        ,types_gofr,iret,errmsg)
         IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
      END IF

      IF(avg_str .OR. avg_rms .OR. (diffusion .AND. slt_exist)) THEN
         ALLOCATE(xp_avg(ntap),yp_avg(ntap),zp_avg(ntap))
         CALL zeroa(xp_avg,yp_avg,zp_avg,ntap,1)
      END IF
      IF(avg_rms) THEN
         ALLOCATE(xp_avg2(ntap),yp_avg2(ntap),zp_avg2(ntap))
         CALL zeroa(xp_avg2,yp_avg2,zp_avg2,ntap,1)
      END IF
      IF(inst_fit .OR. lda_hyd) CALL zeroa(xpo,ypo,zpo,ntap,1)

*=======================================================================
*---- Initialize G of Rs -----------------------------------------------
*=======================================================================

      IF(gofr) THEN
         ALLOCATE(krdf(maxint*g1))
#ifdef PARALLEL
         ALLOCATE(krdfa(maxint*g1))
#endif
         CALL zero_gofr(maxint,krdf,ngrdon,offset_slv)
      END IF

      IF(nplot_fragm .GT. 0) THEN
         do i=1,nfragm
            IF ( fragm(1,i).gt.ntap.or.fragm(2,i).gt.ntap) THEN 
               WRITE(kprint,10260) i
10260          FORMAT(/ /  
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
*---  set fake neighbor lists
*=======================================================================

      CALL setup_fake_map(ngrp,grppt,nnlppf,pfix,mass)

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
*----- Find out what was the timestep ----------------------------------
*=======================================================================

      IF(dmprnd_i) THEN
         CALL read_confc_rows(co,xp0,yp0,zp0,xau,yau,zau,ntap,fstep1,1
     &        ,end,err,divide_records,atom_record)
         IF(end) THEN
            errmsg=' Trajectory file seems to be empty.'
            CALL xerror(errmsg,80,1,2)
         END IF
         CALL read_confc_rows(co,xp0,yp0,zp0,xau,yau,zau,ntap,fstep2,2
     &        ,end,err,divide_records,atom_record)
         IF(end) THEN
            errmsg=' Trajectory file contains only one configuration.'
            CALL xerror(errmsg,80,1,1)
            fstep=fstep1
            stop_anl=1
         ELSE
            fstep=fstep2-fstep1
            IF(dmprnd_o) THEN
               fstep_o=fstep
               nstep_start_o=maxstp
            END IF

*=======================================================================
*----- Find out the length of the run ----------------------------------
*=======================================================================

            naux=stop_anl
            CALL find_length_run(stop_anl)
            WRITE(kprint,2000) stop_anl
            IF(naux .LT. stop_anl) stop_anl=naux
         END IF
      END IF
*========== Compute Native contacts ====================================
                  
      IF(native) THEN
         CALL GetSimpleCO(co,oc)
         CALL change_frame(co,oc,-1,nato_slt,xpt0,ypt0,zpt0,xpa,ypa
     &        ,zpa)
         CALL comp_native(knative,co,xpa,ypa,zpa,beta
     &        ,nbun_slt,atres,native_dist,nres(1,1))
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

      IF(prttpg) THEN
         WRITE(kprint,1000)
         IF(prtseq) CALL prtsq(nbun,mend,prsymb)
         IF(prtatl) CALL prtat(ss_index,ntap,beta,betb,nres,m1,prsymb
     &        ,chrge,mass)
         IF(prtbndl) CALL prtbnd(beta,nres(1,1),lbnd,lbond,potbo
     &        ,m9)
*         IF(prtcnl) CALL prtcn(dssp,lcnstr,lconstr,beta,nres(1,1))
         IF(prtbal) CALL prtba(beta,nres(1,1),lbndg,lbend,potbe,m2)
         IF(prtptl) CALL prtpt(beta,nres(1,1),ltor,ltors,potto,m3)
         IF(prtitl) CALL prtit(beta,nres(1,1),litr,litor,potit,m4)
         WRITE(kprint,'(/ / / /77(''=''))')
      END IF


*=======================================================================
*----- Print titles for the run ----------------------------------------
*=======================================================================

      CALL print_title_analysis(fstep)
      ffstep=fstep
      IF((.NOT. dmprnd_i)) STOP

*=======================================================================
*----- Decide where to start to read -----------------------------------
*=======================================================================

      nstep=0
      IF(.NOT. pdb_read .AND. start_anl .EQ. 0) start_anl=1
      IF(start_anl .NE. 0) nstep=start_anl-1

      CALL timer(vfcp,tfcp,elapse)
      gcpu=tfcp

      WRITE(*,*) start_anl,stop_anl
      IF(ndipole .GT. 0) THEN
         CALL get_memory_iobuffer(start_anl,stop_anl,start_time
     &        ,end_time,length_run,atom_record,length_tot,length_fft)

         ALLOCATE(vxpd(buffer_fft),vypd(buffer_fft),vzpd(buffer_fft))

         CALL dcffti(length_fft,wsave1)
      END IF

      IF(time_corr) THEN

*=======================================================================
*--- Ask for memory ----------------------------------------------------
*=======================================================================

         CALL get_memory_iobuffer(start_anl,stop_anl,start_time
     &        ,end_time,length_run,atom_record,length_tot,length_fft)

         IF(vacf) THEN
            WRITE(kprint,2500)
         END IF
         IF(diffusion) THEN
            WRITE(kprint,2501)
         END IF

*=======================================================================
*--- Check the length of the simulation --------------------------------
*=======================================================================

         IF(start_anl .NE. 0) CALL check_length_sim(start_anl,stop_anl
     &        ,buffer_time,buffer_fft,start_time,end_time,length_run
     &        ,length_tot,length_fft,iret,errmsg)
         IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)

*=======================================================================
*--- Check if buffer is sufficiently large for simulation parameters ---
*=======================================================================

         CALL check_read_columns(nbuffer,length_tot,atom_record,iret
     &        ,errmsg)
         IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)

*=======================================================================
*--- Initialize cubic spline arrays to compute velocities --------------
*=======================================================================

         CALL coordinate_spline_init(spline_x,fstep,length_tot)

*=======================================================================
*--- Reinitialize length of simulation and fft's if necessary ----------
*=======================================================================

         IF(vacf .AND. divide_spline .GT. 1) THEN
            CALL change_step_sim(buffer_time,buffer_fft,fstep
     &           ,length_tot,length_fft,divide_spline,iret,errmsg)
            IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
         END IF

*=======================================================================
*--- Initialize FFT routines -------------------------------------------
*=======================================================================

         CALL dcffti(length_fft,wsave1)

*=======================================================================
*--- Zero correlation arrays -------------------------------------------
*=======================================================================

         DO j=1,2
            CALL zero0(vacf_data,3*length_tot)
            CALL zero0(tot_rms_disp,2*length_tot)
         END DO

         IF(diffusion) THEN
            CALL zero0(xpcc,length_tot)
            CALL zero0(ypcc,length_tot)
            CALL zero0(zpcc,length_tot)
            CALL zero0(RotMat,length_tot*9)

            WRITE(kprint,5000)
            nnstep=0
            DO istep=start_time,end_time
               nnstep=nnstep+1

*=======================================================================
*----- Read trajectory file by row -------------------------------------
*=======================================================================
            
               CALL read_confc_rows(co,xp0,yp0,zp0,xau,yau,zau,ntap
     &              ,fstep1,istep,end,err,divide_records,atom_record)

               CALL calc_RotMat(wca,xpt0,ypt0,zpt0,xp0,yp0,zp0
     &              ,xpcc(nnstep),ypcc(nnstep),zpcc(nnstep)
     &              ,RotMat(1,1,nnstep),xt_cm,yt_cm,zt_cm,nato_slt)

               WRITE(*,*) 'Nstep =',nnstep,istep
            END DO
         END IF

*=======================================================================
*----  Set atom number of solute and solvent to zero and decide --------
*----  for which atom to compute correlation ---------------------------
*=======================================================================

         nato1=0
         nato2=0

         IF(corr_atoms(1) .EQ. 0) THEN
            niatom=ntap
         ELSE
            niatom=corr_atoms(1)
         END IF

*=======================================================================
*---- Beginning of time correlation loop -------------------------------
*=======================================================================

         i_old=0
         DO iatom0=1,niatom
            IF(corr_atoms(1) .EQ. 0) THEN
               iatom=iatom0
            ELSE
               iatom=corr_atoms(1+iatom0)
            END IF
            typei=ss_index(iatom)
            IF(beta(iatom)(1:2) .NE. 'o1') GOTO 12000
*=======================================================================
*---- Read coordinates of each atom as a function of time --------------
*=======================================================================

            CALL read_confc_columns(xpc,ypc,zpc,xpb,ypb,zpb,ntap,iatom
     &           ,length_run,start_time,end_time,divide_records
     &           ,atom_record,iret,errmsg)
            IF(vacf) THEN

*=======================================================================
*---- Compute velocity autocorrelation ---------------------------------
*=======================================================================

*---- Do x component

               CALL coordinate_spline(spline_x,spline_y,xpc,length_tot
     &              ,divide_spline)
               CALL get_velocities(spline_x,spline_y,vxpc,length_tot
     &              ,fstep)
               CALL time_correlation(length_fft,vxpc,vxpc,wsave1
     &              ,vacf_data((typei-1)*buffer_time),fstep,w1,w2)
               
*---- Do y component
               
               CALL coordinate_spline(spline_x,spline_y,ypc,length_tot
     &              ,divide_spline)
               CALL get_velocities(spline_x,spline_y,vypc,length_tot
     &              ,fstep)
               CALL time_correlation(length_fft,vypc,vypc,wsave1
     &              ,vacf_data((typei-1)*buffer_time),fstep,w1,w2)
               
*---- Do z component
               
               CALL coordinate_spline(spline_x,spline_y,zpc,length_tot
     &              ,divide_spline)
               CALL get_velocities(spline_x,spline_y,vzpc,length_tot
     &              ,fstep)
               CALL time_correlation(length_fft,vzpc,vzpc,wsave1
     &              ,vacf_data((typei-1)*buffer_time),fstep,w1,w2)
               
               IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)

            END IF
            IF(diffusion) THEN
               IF(typei .EQ. 1) THEN
                  nato1=nato1+1
               ELSE IF(typei .EQ. 2) THEN
                  nato2=nato2+1
               END IF

*=======================================================================
*---- Compute rms displacement -----------------------------------------
*=======================================================================

               CALL zero0(rms_disp((typei-1)*buffer_time),length_tot)

*---- Do x component
               
               CALL RotCoordColumn(length_tot,xpc,ypc,zpc,xpcc,ypcc,zpcc
     &              ,xt_cm,yt_cm,zt_cm,RotMat)
               CALL dcopy(length_tot,xpc,1,vxpc,1)
               CALL dcopy(length_tot,ypc,1,vypc,1)
               CALL dcopy(length_tot,zpc,1,vzpc,1)
               CALL time_correlation(length_fft,vxpc,vxpc,wsave1
     &              ,rms_disp((typei-1)*buffer_time),fstep,w1,w2)
               
*---- Do y component
               
               CALL time_correlation(length_fft,vypc,vypc,wsave1
     &              ,rms_disp((typei-1)*buffer_time),fstep,w1,w2)
               
*---- Do z component
               
               CALL time_correlation(length_fft,vzpc,vzpc,wsave1
     &              ,rms_disp((typei-1)*buffer_time),fstep,w1,w2)
               
*--- Compute the displacement for the current atom and accumulate
               
               CALL comp_rmsq(beta(iatom),iatom,length_tot
     &              ,rms_disp((typei-1)*buffer_time),xpc,ypc,zpc
     &              ,tot_rms_disp((typei-1)*buffer_time),w1,w2)
            END IF
12000       CONTINUE
            naux=MOD(iatom,atom_record)
            IF(naux .EQ. 0) THEN
               naux=iatom/atom_record
            ELSE
               naux=iatom/atom_record+1
            END IF
            i_start=(naux-1)*atom_record+1
            i_end=naux*atom_record
            IF(i_start .NE. i_old) THEN
               IF(i_end .GT. ntap) i_end=ntap
               WRITE(kprint,3000) i_start,i_end
            END IF
            i_old=i_start
         END DO

         IF(vacf) THEN
            CALL get_spectra_vacf(buffer_time,length_tot,length_tot/2
     &           ,vacf_data,vacf_spectra,wsave1)
            WRITE(kprint,3500)
            CALL write_vacf_phi(buffer_time,length_tot,fstep,vacf_data
     &           ,vacf_spectra)
         END IF
         IF(diffusion) THEN
            WRITE(kprint,4500)
            CALL write_diffusion(buffer_time,length_tot,fstep
     &           ,tot_rms_disp,nato1,nato2)
         END IF
      ELSE IF(.NOT. dmprnd_o) THEN

*=======================================================================
*----- Analyse interparticle properties --------------------------------
*=======================================================================
               
         IF(voronoi .OR. gofr .OR. hbonds_tot .OR. hbonds_res .OR.
     &        cavities .OR. native  .OR. check_native .OR. hhisto .OR. 
     &        efield .OR. prot_hyd .OR. prot_lda .OR. EPotential .OR.
     &        scan_traj) THEN
            nnstep=0
            DO WHILE(.NOT. end .AND. nstep .LT. stop_anl)
               nstep=nstep+1
               nnstep=nnstep+1
*=======================================================================
*----- Read trajectory file by row -------------------------------------
*=======================================================================
               
               IF(start_anl .NE. 0 .AND. dmprnd_i) THEN
                  CALL read_confc_rows(co,xp0,yp0,zp0,xau,yau,zau,ntap
     &                 ,fstep,nstep,end,err,divide_records,atom_record)

                  IF(template) THEN
                     CALL calc_RotMat(whe,xpt0,ypt0,zpt0,xp0,yp0,zp0
     &                    ,Ixpcc,Iypcc,Izpcc,InstRotMat(1,1),xt_cm,yt_cm
     &                    ,zt_cm,nato_slt)
                     CALL RotCoordRow(ntap,xp0,yp0,zp0,Ixpcc,Iypcc,Izpcc
     &                    ,xt_cm,yt_cm,zt_cm,InstRotMat)
                     CALL CalcDistConf(xp0,yp0,zp0,xpt0,ypt0,zpt0,whe
     &                    ,nato_slt,Dist)
                     IF(node .EQ. 0) THEN
                        WRITE(*,90200) DSQRT(Dist)
                     END IF
                  END IF
               ELSE IF(start_anl .NE. 0 .AND. pdb_read) THEN
                  CALL rdcmac(kconf,kprint,nres(1,1),beta,xp0,yp0,zp0
     &                 ,ntap,.FALSE.,1,1,iret,errmsg)
*=======================================================================
*----- Compute the positions of hydrogens if required ------------------
*=======================================================================

                  CALL serchd(xp0,yp0,zp0,ntap,concta,m1,list,hlist
     &                 ,rlist,nlist,iret,errmsg)
                  IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
                  CALL atthd(list,hlist,rlist,nlist,xp0,yp0,zp0,beta
     &                 ,ntap,xpa,ypa,zpa,iret,errmsg)
                  IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
                  CALL brot(read_co,aaxis,baxis,caxis,alf,bet,gam,icl
     &                 ,icm,icn,co,oc,volume)
                  IF(iret .NE. 0) call xerror(errmsg,80,1,21)
                  WRITE(kprint,'(/ /)')
               END IF
               
               IF(.NOT. end) THEN
                  Regular=.TRUE.
                  IF(err) Regular=.FALSE.
*=======================================================================
*----- Check if all coordinates are non zeroes -------------------------
*=======================================================================
                  
                  CALL check_zero_coord(xp0,yp0,zp0,ntap,iret,errmsg)

*=======================================================================
*--- Compute the volume of the system ----------------------------------
*=======================================================================
                  
                  CALL matinv(3,3,co,oc,volume)
                  volume=volume*boxl**3
                  sum_volume=sum_volume+volume

*=======================================================================
*----- Write information about the processed trajectory ----------------
*=======================================================================

                  IF(fstep .NE. ffstep*DFLOAT(nstep)) THEN
                     Regular=.FALSE.
                     fstep=ffstep*DFLOAT(nstep)
                  END IF
                  IF(volume .EQ. 0.0D0 .OR. iret .NE. 0) THEN
                     Regular=.FALSE.
                     fstep=ffstep*DFLOAT(nstep)
                  END IF
                  IF((gofr .OR. voronoi) .AND. update_anl .NE. 0) THEN
                     WRITE(kprint,90000) fstep
                  ELSE IF(hbonds_tot .OR. hbonds_res) THEN
                     IF(update_anl .NE. 0) THEN 
                        IF(MOD(nnstep-1,update_anl) .NE. 0) THEN
                           WRITE(kprint,90000) fstep
                        END IF
                     ELSE 
                        WRITE(kprint,90000) fstep
                     END IF
                  ELSE
                     WRITE(kprint,90000) fstep
                  END IF

*=======================================================================
*----- Scan trajectory -------------------------------------------------
*=======================================================================                  

                  IF(scan_traj) THEN
                     IF(Regular) THEN
                        WRITE(*,'(''---> Configuration is regular'')') 
                        WRITE(99,*) nstep,.TRUE.
                     ELSE
                        IF(volume .EQ. 0.0D0) THEN
                           WRITE(*,'(''---> CO-Matrix is zero!'')') 
                        END IF
                        IF(iret .NE. 0) THEN
                           WRITE(*,'(''---> Coordinates are zero!'')') 
                        END IF
                        IF(volume .NE. 0.0D0 .AND. iret .EQ. 0) THEN
                           WRITE(*,'(''---> Configuration is''
     &                          ,'' not regular'')') 
                        END IF
                        WRITE(99,*) nstep,.FALSE.
                     END IF
                     CYCLE
                  END IF
                  IF(iret .NE. 0) THEN
                     IF(.NOT. skip_step) THEN
                        call xerror(errmsg,80,2,21)
                     ELSE
                        errmsg='Now Skip: Found CO-matrix '
     &                       / /'singular.'
                        call xerror(errmsg,80,1,21)
                     END IF
                  END If
                  IF(volume .EQ. 0.0D0) THEN
                     IF(.NOT. skip_step) THEN
                        errmsg='Now Stop: Found zeros in trajectory'
     &                       / /' file.'
                        call xerror(errmsg,80,2,21)
                     ELSE
                        errmsg='Now Skip: Found zeros in trajectory'
     &                       / /' file.'
                        call xerror(errmsg,80,2,21)
                     END IF
                  END IF
                  IF(.NOT. Regular) THEN
                     IF(skip_step) THEN
                        CYCLE
                     ELSE
                        errmsg='Now Stop: Trajectory'
     &                       / /' file is not regular.'
                        call xerror(errmsg,80,2,21)                        
                     END IF
                  END IF

                  IF(CEN_Center) THEN
                     CALL CEN_Compute(xp0,yp0,zp0,xpa,ypa,zpa,co,oc
     &                    ,mass,ntap)
                  END IF

*=======================================================================
*-------- Calculate group position  ------------------------------------
*=======================================================================

                  CALL appbou(xp0,yp0,zp0,xpg,ypg,zpg,pmass,1,ngrp,grppt
     &                 )
                  
*=======================================================================
*-------- Calculate solute center of mass ------------------------------
*=======================================================================
                  
                  CALL inicmp(ss_index,xp0,yp0,zp0,xpcm,ypcm,zpcm,mass
     &                 ,nprot,protl)
                  
*=======================================================================
*--- Change frame to get xpa, ypa, zpa etc in box fractions ------------
*=======================================================================
                  
                  CALL change_frame(co,oc,-1,ntap,xp0,yp0,zp0,xpa,ypa
     &                 ,zpa)
                  CALL change_frame(co,oc,-1,ngrp,xpg,ypg,zpg,xpga,ypga
     &                 ,zpga)
                  CALL change_frame(co,oc,-1,nprot,xpcm,ypcm,zpcm,xpcma
     &                 ,ypcma,zpcma)
                  
*=======================================================================
*--- Compute neighbor lists --------------------------------------------
*=======================================================================

                  IF(gofr .OR. voronoi .OR. native.OR. check_native .OR.
     &                 efield .OR. prot_hyd  .OR. prot_lda .OR.
     &                 EPotential) THEN
                     IF(update_anl .NE. 0) THEN
                        IF(MOD(nnstep-1,update_anl) .EQ. 0) THEN
                           aux=rcuth+rtolh+rneih
                           IF(linked_cell) THEN
                              CALL lc_index(indxyz,ncx,ncy,ncz,nind
     &                             ,indxi,indxj,indxk,aux,co)
                              IF(voronoi) THEN
                                 CALL lc_list(ncx,ncy,ncz,nind,indxi
     &                                ,indxj,indxk,aux,co,xpga,ypga,zpga
     &                                ,ngrp,nstart_h,nend_h,node,nprocs
     &                                ,ncube,nnlpp0,npp_u,npp,worka
     &                                ,kprint,.FALSE.)
                              ELSE IF(prot_hyd) THEN
                                 CALL lc_list2(ncx,ncy,ncz,nind,indxi
     &                                ,indxj,indxk,aux,co,xpga,ypga,zpga
     &                                ,ngrp,nstart_h,ngrp-nmol
     &                                ,node,nprocs
     &                                ,ncube,nnlpp0,npp_u,npp,worka
     &                                ,kprint,ngrp-nmol)
                              ELSE IF(prot_lda) THEN
                                 CALL lc_list3(ncx,ncy,ncz,nind,indxi
     &                                ,indxj,indxk,aux,co,xpga,ypga,zpga
     &                                ,ngrp,nstart_h,nend_lda
     &                                ,node,nprocs
     &                                ,ncube,nnlpp0,npp_u,npp,worka
     &                                ,kprint,ngrp-nmol)
                              ELSE
                                 CALL lc_list(ncx,ncy,ncz,nind,indxi
     &                                ,indxj,indxk,aux,co,xpga,ypga,zpga
     &                                ,ngrp,nstart_h,nend_h,node,nprocs
     &                                ,ncube,nnlpp0,npp_u,npp,worka
     &                                ,kprint,.TRUE.)
                              END IF 
                           ELSE
*----- Compute neighbors
                              P_shell='l'
                              IF(voronoi) THEN
                                 CALL mts_forces('v',xpa,ypa,zpa,xpga
     &                                ,ypga,zpga,xpcma,ypcma,zpcma,mapnl
     &                                ,mapdn,nmapdn,ucns_p,ucos_p,virs_p
     &                                ,virsp_p,ucnp_p,ucop_p,ucnsp_p
     &                                ,ucosp_p,fpx_p,fpy_p,fpz_p
     &                                ,stressd_p,nnlppf,nnlpp0,npp_u,npp
     &                                ,worka,cpu_h,ncpu_h,nstart_h
     &                                ,nend_h,nstart_ah,nend_ah
     &                                ,nlocal_ah,node,nprocs,ncube
     &                                ,P_shell)
                              ELSE
                                 CALL mts_forces('u',xpa,ypa,zpa,xpga
     &                                ,ypga,zpga,xpcma,ypcma,zpcma,mapnl
     &                                ,mapdn,nmapdn,ucns_p,ucos_p,virs_p
     &                                ,virsp_p,ucnp_p,ucop_p,ucnsp_p
     &                                ,ucosp_p,fpx_p,fpy_p,fpz_p
     &                                ,stressd_p,nnlppf,nnlpp0,npp_u,npp
     &                                ,worka,cpu_h,ncpu_h,nstart_h
     &                                ,nend_h,nstart_ah,nend_ah
     &                                ,nlocal_ah,node,nprocs,ncube
     &                                ,P_shell)
                              END IF
                           END IF

                        END IF
                     END IF
                  ELSE IF(hbonds_tot .OR. hbonds_res .OR. hhisto) THEN
                     IF(update_anl .NE. 0) THEN
                        IF(MOD(nnstep-1,update_anl) .EQ. 0) THEN
                           CALL update(co,xpa,ypa,zpa,llacc,lacc,lldon
     &                          ,ldon,hrcut_update,nmapdn,mapdn,mf)
                           IF(pdb_read) THEN
                              WRITE(kprint,90001) 
                           ELSE
                              WRITE(kprint,90000) fstep
                           END IF
                        END IF
                     END IF
                  END IF

*========== Compute Native contacts ====================================
                  
                  IF(native) THEN
                     IF(MOD(nstep,nnative) .EQ. 0) THEN
                        CALL comp_native(knative,co,xpa,ypa,zpa,beta
     &                       ,nbun_slt,atres,native_dist,nres(1,1))
                     END IF
                  END IF

*========== Compute Contact for the trajectory =========================
                  
                  IF(check_native) THEN
                     IF(MOD(nstep,nnative) .EQ. 0) THEN
                        CALL comp_contacts(co,xpa,ypa,zpa,beta
     &                       ,native_theta,atres,nat_listp,nat_list
     &                       ,native_dist,nat_cntactp,nat_cntact)
                        CALL write_contacts(knative,fstep,nat_listp
     &                       ,nat_list,nat_cntactp,nat_cntact)
                     END IF
                  END IF

*========== Compute Voronoi volumes and areas ==========================
                  
                  IF(voronoi) THEN
                     IF(MOD(nstep,nvoronoi) .EQ. 0) THEN
                        WRITE(*,93000) 
                        CALL zero_voronoi
                        CALL comp_neigh_vor(nstart_h,nend_h,nstart_ah
     &                       ,nend_ah,heavy_vor,beta,ntap,ngrp,grppt
     &                       ,nnlpp0,cutoff_vor,xpa,ypa,zpa,xpga,ypga
     &                       ,zpga,co,iret,errmsg)
#ifdef PARALLEL
                        CALL P_get_errmsg(iret,errmsg,80,node,nprocs
     &                       ,ncube,nbyte)
#else
                        IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
#endif
                        CALL comp_voronoi(nstart_ah,nend_ah,ntap,xpga
     &                       ,ypga,zpga,atomg,xpa,ypa,zpa,co,iret,errmsg
     &                       )
                        IF(voronoi_volume .OR. voronoi_access .OR.
     &                       voronoi_contact .OR. voronoi_neighbor) THEN
                           CALL analyse_voronoi(nstart_ah,nend_ah
     &                          ,nlocal_ah,nstart_uh,nend_uh,nlocal_uh
     &                          ,node,nprocs,ncube,fstep,volume
     &                          ,voronoi_volume,voronoi_access
     &                          ,voronoi_contact,voronoi_neighbor
     &                          ,ncontact_slt,contact_slt,voronoi_res
     &                          ,kvoronoi,ss_index,ss_point(1,1),grppt
     &                          ,mend,protl,nprot,ntap,nbun,beta,atomp
     &                          ,nres(1,1),nres(1,2),mres,prsymb,iret
     &                          ,errmsg)
                           IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
                        END IF
                     END IF
                  END IF
                  
*========== Compute Cavity distribution ================================
                  IF(cavities) THEN
                     IF(MOD(nstep,nvoronoi) .EQ. 0) THEN
                        CALL comp_cavities_vor(nstart_ah,nend_ah
     &                       ,ss_index,nx_cav,ny_cav,nz_cav
     &                       ,size_atom_cav,bin_size_cav,rmax_size_cav
     &                       ,mqq,jmax_cav,xpa,ypa,zpa,xpga,ypga,zpga
     &                       ,ntap,atomg,nbtype,nres(1,1),mback,nbone
     &                       ,nbun,pnbd1,co,iret,errmsg)
#ifdef PARALLEL
                        CALL P_get_errmsg(iret,errmsg,80,node,nprocs
     &                       ,ncube,nbyte)
#else
                        IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
#endif
                     END IF
                     IF(MOD(nstep,ncavities) .EQ. 0) THEN
                        WRITE(kprint,6000) 
                        CALL write_cavities_vor(kcavities,cavities_file
     &                       ,prsymb,mend,jmax_cav,bin_size_cav
     &                       ,rmax_size_cav,cavity_h,maxcav_bin,node
     &                       ,nprocs,ncube)
                     END IF
                  END IF
*========== Compute Hydrogen bonds =====================================
                  
                  IF(hbonds_tot .OR. hbonds_res .OR. hhisto) THEN
                     IF(nhbonds .NE. 0) THEN
                        IF(MOD(nstep,nhbonds) .EQ. 0)  THEN
                           CALL write_hbonds(update_anl,hbonds_tot
     &                          ,hbonds_res,fstep,khbonds,ss_index,nbun
     &                          ,nres,m1,atomp,co,xpa,ypa,zpa,lacc,ldon
     &                          ,llacc,lldon,rcut_hb,acut_hb,a2cut_hb
     &                          ,nmapdn,mapdn,hhisto,hhisto_list
     &                          ,hhisto_bin,hhisto_dim,hhisto_count,beta
     &                          ,prsymb,pdb_read)
                        END IF 
                     END IF
                     IF(hhisto) THEN
                        IF(nhhisto .NE. 0) THEN
                           IF(MOD(nstep,nhhisto) .EQ. 0) THEN
                              CALL write_hhisto(hhisto_list,hhisto_bin
     &                             ,hhisto_dim,hhisto_count)
                           END IF
                        END IF
                     END IF
                  END IF
                  
*========== Compute G of R's of solute and solvent =====================
                  
                  IF(gofr) THEN
                     IF(MOD(nstep,gofr_ncomp) .EQ. 0) THEN
                        l2=maxint
                        CALL calc_gofr(nato_slt,nato_slv,type_slv
     &                       ,ss_index,atomp,gofr_neighbor,gofr_intra,co
     &                       ,xpa,ypa,zpa,xpga,ypga,zpga,wca,whe,delrg
     &                       ,nstart_h,nend_h,nstart_ah,nend_ah,ntap
     &                       ,ngrp,grppt,l2,krdf,ngrdon,gofr_cut)
                        IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
                     END IF
                  END IF
                  
*========== Compute s of k's of solute only ============================
                  
                  IF(s_of_k) THEN
                     IF(MOD(nstep,sofk_ncomp) .EQ. 0 .AND. nstep .GT.
     &                    nrject)THEN
                        CALL calc_sk(ss_index,co,oc,xpa,ypa,zpa,whe,ntap
     &                       ,sofk_cut,sofk_delta,sofk_data,nsofk_data
     &                       ,maxsk,node,nprocs)
c$$$                        CALL calc_sk_eet(xpa,ypa,zpa,whe,nstart_ah
c$$$     &                       ,nend_ah,sofk_cut,sofk_delta,sofk_data
c$$$     &                       ,nsofk_data,maxsk)
                     END IF
                     IF(MOD(nstep,sofk_nprint) .EQ. 0 .AND. nstep .GT.
     &                    nrject)THEN
                        CALL write_sk(ksofk,fstep,sofk_cut,sofk_delta
     &                       ,sofk_data,nsofk_data,maxsk)
                     END IF
                  END IF
            
*========== Print G of R's of solute and solvent =======================
                  
                  IF(gofr) THEN
                     IF(MOD(nstep,gofr_nprint) .EQ. 0)THEN
c                  vol_gofr=sum_volume/DFLOAT(nstep)
                        vol_gofr=volume
#ifdef PARALLEL 
                        CALL P_get_gofr(krdf,krdfa,worka,maxint
     &                       ,offset_slv,node,nprocs,ncube,nbyte)
#endif
                        IF(node .EQ. 0) THEN
                           IF(slt_exist) THEN
                              offset=0
                              CALL write_gofrp(.NOT.gofr_avg,fstep
     &                             ,_KRDF_,maxint,0,wca,whe,nato_slt
     &                             ,delrg,gofr_cut,ngrdon,iret,errmsg)
                              IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
                           END IF
                           IF(slv_exist .AND. slt_exist) THEN
                              CALL write_gofrw(.FALSE.,fstep,_KRDF_
     &                             ,maxint,g1,delrg,gofr_cut,itype_slv
     &                             ,betab_slv,vol_gofr,ntype_slv,nmol
     &                             ,ngrdon,4,iret,errmsg)
                              IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
                           ELSE IF(slv_exist .AND. (.NOT. slt_exist))
     &                             THEN
                              CALL write_gofrw(.TRUE.,fstep,_KRDF_
     &                             ,maxint,g1,delrg,gofr_cut,itype_slv
     &                             ,betab_slv,vol_gofr,ntype_slv,nmol
     &                             ,ngrdon,4,iret,errmsg)
                              IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
                           END IF
                        END IF
                     END IF
                  END IF
                  
*=======================================================================
*---- Initialize G of Rs when required ---------------------------------
*=======================================================================
                  
                  IF(gofr .AND. gofr_avg) THEN
                     IF(MOD(nstep,gofr_navg) .EQ. 0) THEN
                        CALL zero_gofr(maxint,krdf,ngrdon,offset_slv)
                     END IF
                  END IF
 
*=======================================================================
*---- Calculate Electric field when required ---------------------------
*=======================================================================

*=======================================================================
*---- Calculate Electric field when required ---------------------------
*=======================================================================


*=======================================================================
*---- Compute hydration for the protein residue ------------------------
*---- NO parallelization
*=======================================================================
                     IF(prot_hyd) THEN


                        CALL calc_hyd(kprot_hyd,nprot_hyd,nstep,fstep
     &                       ,nres,nbun,nbtype,pnbd1
     &                       ,coeff_hyd,min_hyd,max_hyd,nnlpp0
     &                       ,nstart_h,nend_h-nmol,grppt,co
     &                       ,xpa,ypa,zpa,nato_slt,nmol
     &                       ,res_time,kprot_rest)
                     END IF

*=======================================================================
*---- Compute lda-solvation for the protein residue --------------------
*---- NO parallelization
*=======================================================================
                     IF(prot_lda) THEN
                        CALL calc_lda_rest(kprot_lda,nprot_lda,nstep
     &                       ,fstep,nres,nbun,nbtype,pnbd1
     &                       ,coeff_lda,min_lda,max_lda,nnlpp0
     &                       ,nstart_h,nend_lda,grppt,co
     &                       ,xpa,ypa,zpa,res_time,klda_rest)
                     END IF

c--------------------------
               ELSE
                  nstep=nstep-1
               END IF
            END DO
         ELSE
            nnstep=0
            DO WHILE(.NOT. end .AND. nstep .LT. stop_anl)
               nstep=nstep+1
               nnstep=nnstep+1
            
*=======================================================================
*----- Read trajectory file by row -------------------------------------
*=======================================================================
            
               IF(start_anl .NE. 0 .AND. dmprnd_i) THEN
                  CALL read_confc_rows(co,xp0,yp0,zp0,xau
     &                 ,yau,zau,ntap,fstep,nstep,end,err,divide_records
     &                 ,atom_record)
               ELSE IF(pdb_read .AND. start_anl .NE. 0) THEN
                  CALL rdcmac(kconf,kprint,nres(1,1),beta,xp0,yp0,zp0
     &                 ,ntap,.FALSE.,1,1,iret,errmsg)
                  CALL brot(read_co,aaxis,baxis,caxis,alf,bet,gam,icl
     &                 ,icm,icn,co,oc,volume)
                  IF(iret .NE. 0) call xerror(errmsg,80,1,21)
                  WRITE(kprint,'(/ /)')
               END IF
               
               IF(.NOT. end) THEN
                  
*=======================================================================
*----- Check if all coordinates are non zeroes -------------------------
*=======================================================================
                  
                  CALL check_zero_coord(xp0,yp0,zp0,ntap,iret,errmsg)
                  IF(iret .NE. 0) call xerror(errmsg,80,1,21)
                  IF(iret .EQ. 1) CYCLE

*=======================================================================
*--- Compute the volume of the system ----------------------------------
*=======================================================================
                  
                  CALL matinv(3,3,co,oc,volume)
                  volume=volume*boxl**3
                  sum_volume=sum_volume+volume
                  
*=======================================================================
*----- Write information about the processed trajectory ----------------
*=======================================================================
                  
                  WRITE(kprint,90000) fstep
*=======================================================================
*-------- Calculate group position  ------------------------------------
*=======================================================================
                  
                  IF(CEN_Center) THEN
                     CALL CEN_Compute(xp0,yp0,zp0,xpa,ypa,zpa,co,oc
     &                    ,mass,ntap)
                  END IF

                  CALL appbou(xp0,yp0,zp0,xpg,ypg,zpg,pmass,1,ngrp,grppt
     &                 )
                  
*=======================================================================
*-------- Calculate solute center of mass ------------------------------
*=======================================================================
                  
                  CALL inicmp(ss_index,xp0,yp0,zp0,xpcm,ypcm,zpcm,mass
     &                 ,nprot,protl)
                  
*=======================================================================
*--- Change frame to get xpa, ypa, zpa etc in box fractions ------------
*=======================================================================
                  
                  CALL change_frame(co,oc,-1,ntap,xp0,yp0,zp0,xpa,ypa
     &                 ,zpa)
                  CALL change_frame(co,oc,-1,ngrp,xpg,ypg,zpg,xpga,ypga
     &                 ,zpga)
                  CALL change_frame(co,oc,-1,nprot,xpcm,ypcm,zpcm,xpcma
     &                 ,ypcma,zpcma)
                  
*=======================================================================
*---- Compute instantaneous X-rms --------------------------------------
*=======================================================================

                  IF(hydration) THEN
                     IF(MOD(nstep,HYD_n_neighbors) == 0) THEN
                        CALL HYD_Compute_Neighbors(xpga,ypga,zpga,co)
                     END IF
                     CALL HYD_Compute(xpa,ypa,zpa,co,nbtype,pnbd1)
                     IF(MOD(nstep,HYD_n_write) == 0) THEN
                        CALL HYD_Write_it(fstep)
                     END IF
                  END IF

                  IF(node .EQ. 0) THEN
                     IF(Density_Calc) THEN                     
                        CALL DEN_Compute(xpa,ypa,zpa,co)
                        IF(MOD(nstep,DEN_write) == 0) THEN
                           CALL DEN_write_it(Volume,co)
                        END IF
                     END IF
                  END IF
                  IF(node .EQ. 0) THEN
                     IF(GE_Groups) THEN
                        CALL GE_Compute(xpa,ypa,zpa,mass,co)
                        IF(MOD(nstep,GE_write) == 0) THEN
                           CALL GE_output(fstep)
                        END IF
                     END IF
                  END IF
                  IF(node .EQ. 0) THEN
                     IF(SUB_Subtract) THEN
                        IF(nnstep == 1) CALL SUB_Start(xp0,yp0,zp0)
                        CALL SUB_Rotate(xp0,yp0,zp0)
                        CALL SUB_Compute(co,oc,xp0,yp0,zp0)
                        IF(MOD(nstep,SUB_write) == 0) THEN
                           CALL SUB_write_it(fstep)
                        END IF
                     END IF
                  END IF
                  IF(node .EQ. 0) THEN
                     IF(PDB_pdbs) THEN
                        IF(MOD(nstep,PDB_ncompute) == 0) THEN
                           CALL PDB_Compute(xp0,yp0,zp0)
                        END IF
                        IF(MOD(nstep,PDB_nwrite) == 0) THEN
                           CALL PDB_write(fstep)
                        END IF
                     END IF
                  END IF

                  IF(node .EQ. 0) THEN
                     IF(rms_matrix) THEN
                        CALL RMS_Rotate(xp0,yp0,zp0)
                        IF(rms_matrix_avg) THEN
                           CALL RMS_Compute_avg(xp0,yp0,zp0)
                           CALL RMS_Write_it_avg
                        ELSE
                           CALL RMS_Compute(xp0,yp0,zp0)
                           CALL RMS_Write_it
                           IF(rms_matrix_plot) CALL RMS_plot
                        END IF
                     END IF
                     IF(eul_angles) THEN
                        CALL EUL_Compute(xp0,yp0,zp0)
                     END IF
                     IF(anxrms) THEN
                        CALL CalcXrmsSecStruct(anxca,anxbc,anxhe,anxal
     &                       ,SecStructure,SecStructTotal,SecPointer
     &                       ,grppt,mres,protl,wca,whe,wbc,xp0,yp0,zp0
     &                       ,xpt0,ypt0,zpt0,nato_slt,errca,errhe,errbc
     &                       ,erral,drpca,drpbc,drphe,drpal)
                        
                        IF(nxrms .NE. 0 .AND.
     &                       MOD(nstep,nxrms) .EQ. 0) THEN
                           IF(SecStructure) THEN
                              WRITE(kxrms,50000) fstep
                              IF(anxca) CALL write_xrms(kxrms
     &                             ,SecStructTotal,'CA',errca)
                              IF(anxbc) CALL write_xrms(kxrms
     &                             ,SecStructTotal,'BC',errbc)
                              IF(anxhe) CALL write_xrms(kxrms
     &                             ,SecStructTotal,'HE',errhe)
                              IF(anxal) CALL write_xrms(kxrms
     &                             ,SecStructTotal,'AL',erral)
                           ELSE
                              WRITE(kxrms,50000) fstep
                              IF(anxca) CALL write_xrms(kxrms,nprot,'CA'
     &                             ,errca)
                              IF(anxbc) CALL write_xrms(kxrms,nprot,'BC'
     &                             ,errbc)
                              IF(anxhe) CALL write_xrms(kxrms,nprot,'HE'
     &                             ,errhe)
                              IF(anxal) CALL write_xrms(kxrms,nprot,'AL'
     &                             ,erral)
                           END IF
                           REWIND kxrms_atm
                           IF(anxca) CALL write_xrms_atm(kxrms_atm,ntap
     &                          ,'CA',drpca,nnstep,fstep,ngrp,grppt
     &                          ,nres(1,1))
                           IF(anxbc) CALL write_xrms_atm(kxrms_atm,ntap
     &                          ,'BC',drpbc,nnstep,fstep,ngrp,grppt
     &                          ,nres(1,1))
                           IF(anxhe) CALL write_xrms_atm(kxrms_atm,ntap
     &                          ,'HE',drphe,nnstep,fstep,ngrp,grppt
     &                          ,nres(1,1))
                           IF(anxal) CALL write_xrms_atm(kxrms_atm,ntap
     &                          ,'AL',drpal,nnstep,fstep,ngrp,grppt
     &                          ,nres(1,1))
                        END IF

                        IF(nxslt .NE. 0) THEN
                           IF(MOD(nstep,nxslt) .EQ. 0) THEN
                              
                              CALL calc_xslt(anxca,anxbc,anxhe,anxal
     &                             ,anprot,annpro,anpoint,protl,wca,whe
     &                             ,wbc,xp0,yp0,zp0,xpt0,ypt0,zpt0
     &                             ,nato_slt,errca,errhe,errbc,erral)

                              IF(anprot) THEN
                                 WRITE(kxslt,50000) fstep
                                 IF(anxca) CALL write_xrms(kxslt,annpro
     &                                ,'CA',errca)
                                 IF(anxbc) CALL write_xrms(kxslt,annpro
     &                                ,'BC',errbc)
                                 IF(anxhe) CALL write_xrms(kxslt,annpro
     &                                ,'HE',errhe)
                                 IF(anxal) CALL write_xrms(kxslt,annpro
     &                                ,'AL',erral)
                              ELSE
                                 WRITE(kxslt,50000) fstep
                                 IF(anxca) CALL write_xrms(kxslt,nprot,'
     &                                CA',errca)
                                 IF(anxbc) CALL write_xrms(kxslt,nprot,'
     &                                BC',errbc)
                                 IF(anxhe) CALL write_xrms(kxslt,nprot,'
     &                                HE',errhe)
                                 IF(anxal) CALL write_xrms(kxslt,nprot,'
     &                                AL',erral)
                              END IF
                           END IF
                        END IF
                     END IF
                     
*=======================================================================
*---- Compute properties (cofm)  ---------------------------------------
*=======================================================================

                     IF(ncalc_cofm .NE. 0) THEN
                        IF(MOD(nstep,ncalc_cofm) .EQ.0) THEN
                           
                           IF(avg_ca) CALL calccofm(kcalc_cofm,latms
     &                          ,natms,patms,wca,xpt0,ypt0,zpt0,qt0,xp0
     &                          ,yp0,zp0,mass,nato_slt,fstep)
                           IF(avg_bc) CALL calccofm(kcalc_cofm,latms
     &                          ,natms,patms,wbc,xpt0,ypt0,zpt0,qt0,xp0
     &                          ,yp0,zp0,mass,nato_slt,fstep)
                           IF(avg_he) CALL calccofm(kcalc_cofm,latms
     &                          ,natms,patms,whe,xpt0,ypt0,zpt0,qt0,xp0
     &                          ,yp0,zp0,mass,nato_slt,fstep)
                        END IF
                     END IF
                     
*=======================================================================
*---- Compute lda properties -------------------------------------------
*=======================================================================
                     
                     IF(ninst_lda.NE. 0) THEN
                        IF(MOD(nstep,ninst_lda) .EQ.0) THEN
                           
                           IF(avg_ca) CALL calc_lda(klda_inst
     &                          ,nlda_mol,nlda_atm,nlda_zero
     &                          ,wca,xpt0,ypt0,zpt0,qt0,xp0,yp0,zp0,mass
     &                          ,nato_slt,fstep)
                           IF(avg_bc) CALL calc_lda(klda_inst
     &                          ,nlda_mol,nlda_atm,nlda_zero
     &                          ,wbc,xpt0,ypt0,zpt0,qt0,xp0,yp0,zp0,mass
     &                          ,nato_slt,fstep)
                           IF(avg_he) CALL calc_lda(klda_inst
     &                          ,nlda_mol,nlda_atm,nlda_zero
     &                          ,whe,xpt0,ypt0,zpt0,qt0,xp0,yp0,zp0,mass
     &                          ,nato_slt,fstep)
                           IF(avg_lda) THEN
                              
                              IF(avg_ca) CALL calc_lda_avg(klda_rmin
     &                             ,klda_eend,nlda_mol,nlda_atm
     &                             ,nlda_zero,wca,xpt0,ypt0,zpt0,qt0,xp0
     &                             ,yp0,zp0,nato_slt,nstep,ninst_lda
     &                             ,fstep)
                              IF(avg_bc) CALL calc_lda_avg(klda_rmin
     &                             ,klda_eend,nlda_mol,nlda_atm
     &                             ,nlda_zero,wbc,xpt0,ypt0,zpt0,qt0,xp0
     &                             ,yp0,zp0,nato_slt,nstep,ninst_lda
     &                             ,fstep)
                              IF(avg_he) CALL calc_lda_avg(klda_rmin
     &                             ,klda_eend,nlda_mol,nlda_atm
     &                             ,nlda_zero,whe,xpt0,ypt0,zpt0,qt0,xp0
     &                             ,yp0,zp0,nato_slt,nstep,ninst_lda
     &                             ,fstep)
                           END IF
                           
                        END IF
                        
                     END IF
c----------
                     IF(lda_hyd) THEN
                        CALL calc_lda_hyd(klda_hyd,nlda_hyd,nstep,fstep
     &                       ,nlda_mol,nlda_atm,nlda_zero,co,xpa,ypa,zpa
     &                       ,nato_slt,nato_slv,nmol)
                     END IF
                     
                     IF(lda_flu) THEN
                        CALL calc_lda_flu(klda_flu,nlda_flu,nstep,fstep
     &                       ,nlda_mol,nlda_atm,nlda_zero,xp0,yp0,zp0)
                     END IF
*=======================================================================
*---- Compute instantaneous fit  ---------------------------------------
*=======================================================================
                     
                     IF(ninst_fit .NE. 0) THEN
                        IF(MOD(nstep,ninst_fit) .EQ.0) THEN
                           
                           IF(avg_ca) WRITE(kfit,80000)
                           IF(avg_he) WRITE(kfit,80100)
                           IF(avg_bc) WRITE(kfit,80101)
                           fact=0.0D0
                           CALL dscal(nato_slt,fact,xpo,1)
                           CALL dscal(nato_slt,fact,ypo,1)
                           CALL dscal(nato_slt,fact,zpo,1)
                           IF(avg_ca) CALL calc_inst_fit(anprot,annpro
     &                          ,anpoint,protl,wca,xpt0,ypt0,zpt0,xpo
     &                          ,ypo,zpo,qt0,xp0,yp0,zp0,nato_slt)
                           IF(avg_he) CALL calc_inst_fit(anprot,annpro
     &                          ,anpoint,protl,whe,xpt0,ypt0,zpt0,xpo
     &                          ,ypo,zpo,qt0,xp0,yp0,zp0,nato_slt)
                           IF(avg_bc) CALL calc_inst_fit(anprot,annpro
     &                          ,anpoint,protl,wbc,xpt0,ypt0,zpt0,xpo
     &                          ,ypo,zpo,qt0,xp0,yp0,zp0,nato_slt)
                           CALL plotd(fstep,kfit,beta,xpt0,ypt0,zpt0,xpo
     &                          ,ypo,zpo,nato_slt,nres,m1,prsymb)
                        END IF
                     END IF
                     
*=======================================================================
*---- Compute averaged structure ---------------------------------------
*=======================================================================
                     
                     IF(avg_str .OR. avg_rms) THEN
                        IF(avg_ca) CALL calc_avg_str(anprot,annpro
     &                       ,anpoint,protl,wca,xpt0,ypt0,zpt0,xp_avg
     &                       ,yp_avg,zp_avg
     &                       ,qt0,xp0,yp0,zp0,nato_slt,iter_avg)
                        IF(avg_he) CALL calc_avg_str(anprot,annpro
     &                       ,anpoint,protl,whe,xpt0,ypt0,zpt0,xp_avg
     &                       ,yp_avg,zp_avg
     &                       ,qt0,xp0,yp0,zp0,nato_slt,iter_avg)
                        IF(avg_bc) CALL calc_avg_str(anprot,annpro
     &                       ,anpoint,protl,wbc,xpt0,ypt0,zpt0,xp_avg
     &                       ,yp_avg,zp_avg
     &                       ,qt0,xp0,yp0,zp0,nato_slt,iter_avg)
                        
                     END IF
                     
                     IF(navg_str .NE. 0) THEN
                        IF(MOD(nstep,navg_str) .EQ.0) THEN
                           CALL dcopy(nato_slt,xp_avg,1,xpo,1)
                           CALL dcopy(nato_slt,yp_avg,1,ypo,1)
                           CALL dcopy(nato_slt,zp_avg,1,zpo,1)
                           fact=1.0D0/DFLOAT(iter_avg)
                           CALL dscal(nato_slt,fact,xpo,1)
                           CALL dscal(nato_slt,fact,ypo,1)
                           CALL dscal(nato_slt,fact,zpo,1)
                           IF(avg_ca) WRITE(kavg,80000) 
                           IF(avg_he) WRITE(kavg,80100)
                           IF(avg_bc) WRITE(kavg,80101) 
                           CALL plotd(fstep,kavg,beta,xpt0,ypt0,zpt0,xpo
     &                          ,ypo,zpo,nato_slt,nres,m1,prsymb)
                        END IF
                     END IF
                     
                     IF(navg_str_xrms .NE. 0) THEN
                        IF(MOD(nstep,navg_str_xrms) .EQ.0) THEN
                           CALL dcopy(nato_slt,xp_avg,1,xpo,1)
                           CALL dcopy(nato_slt,yp_avg,1,ypo,1)
                           CALL dcopy(nato_slt,zp_avg,1,zpo,1)
                           fact=1.0D0/DFLOAT(iter_avg)
                           CALL dscal(nato_slt,fact,xpo,1)
                           CALL dscal(nato_slt,fact,ypo,1)
                           CALL dscal(nato_slt,fact,zpo,1)
                           CALL calc_avg_xrms(avg_ca,avg_he,avg_bc,fstep
     &                          ,kavg_xrms,xpt0,ypt0,zpt0,xpo,ypo,zpo
     &                          ,wca,whe,wbc,protl
     &                          ,anprot,annpro,anpoint,nato_slt)
                        END IF
                     END IF
                     
*=======================================================================
*---- Compute average RMS ----------------------------------------------
*=======================================================================
                     
                     IF(avg_rms) THEN
                        IF(avg_ca) CALL calc_avg2_str(anprot,annpro
     &                       ,anpoint,protl,wca,xpt0,ypt0,zpt0,xp_avg2
     &                       ,yp_avg2,zp_avg2
     &                       ,qt0,xp0,yp0,zp0,nato_slt,iter_avg2)
                        IF(avg_he) CALL calc_avg2_str(anprot,annpro
     &                       ,anpoint,protl,whe,xpt0,ypt0,zpt0,xp_avg2
     &                       ,yp_avg2,zp_avg2
     &                       ,qt0,xp0,yp0,zp0,nato_slt,iter_avg2)
                        IF(avg_bc) CALL calc_avg2_str(anprot,annpro
     &                       ,anpoint,protl,wbc,xpt0,ypt0,zpt0,xp_avg2
     &                       ,yp_avg2,zp_avg2
     &                       ,qt0,xp0,yp0,zp0,nato_slt,iter_avg2)
                     END IF
                     IF(nrms .NE. 0) THEN
                        IF(MOD(nstep,nrms) .EQ.0) THEN
                           REWIND krms
                           CALL write_rms(krms,xp_avg,yp_avg,zp_avg
     &                          ,xp_avg2,yp_avg2,zp_avg2,iter_avg,fstep
     &                          ,ngrp,grppt,nres(1,1),beta)
                        END IF 
                     END IF
                     
*=========== Compute dipole ============================================
                     
                     IF(ndipole.gt.0) THEN
                        IF(MOD(nstep,ndipole).EQ.0)THEN
                           CALL comp_dip(ss_index,co,xpga,ypga,zpga,xpa
     &                          ,ypa,zpa,chrge,dips,ntap,ngrp,grppt)
c                           aux=elechg*unitl/3.336D-30
                           aux=1.0D0
                           vxpc(nstep)=dips(1,1)*aux
                           vypc(nstep)=dips(2,1)*aux
                           vzpc(nstep)=dips(3,1)*aux
                           vxpd(nstep) =dips(1,2)*aux
                           vypd(nstep) =dips(2,2)*aux
                           vzpd(nstep) =dips(3,2)*aux
                           WRITE(kdipole,106) fstep, (dips(j,1)*aux
     &                          ,dips(j,2)*aux,j=1,3)
106                        FORMAT(' Dip. ',f12.3,3e15.5,4x,3e15.5)
                        END IF
                     END IF
                     
                     
*=========== Plot .pdb file ============================================
                     
                     IF(nascii .NE. 0) THEN
                        IF(MOD(nstep,nascii) .EQ. 0) THEN
                           IF(ascii_nocell) THEN
                              CALL change_frame(co,oc,1,ntap,xpa,ypa,zpa
     &                             ,xpo,ypo,zpo)
                           ELSE IF(ascii_wsc) THEN
                              CALL tr_wsc(co,xpa,ypa,zpa,xpo,ypo,zpo
     &                             ,mass,nprot,protl)
                              CALL change_frame(co,oc,1,ntap,xpo,ypo,zpo
     &                             ,xpo,ypo,zpo)
                           ELSE
                              CALL tr_inbox(xpa,ypa,zpa,xpo,ypo,zpo,mass
     &                             ,nprot,protl)
                              CALL change_frame(co,oc,1,ntap,xpo,ypo,zpo
     &                             ,xpo,ypo,zpo)
                           END IF
                           CALL plotc(co,abmd,gr,gra,fstep,beta,xpo,ypo
     &                          ,zpo,ntap,nres,m1,prsymb,chrge)
                        END IF
                     END IF
                     
*=========== plot a fragment ===========================================
                     
                     IF(nplot_fragm .GT. 0) THEN
                        IF(MOD(nstep,nplot_fragm).EQ.0 ) THEN
                           CALL tr_wsc(co,xpa,ypa,zpa,xpo,ypo,zpo
     &                          ,mass,nprot,protl)
                           CALL change_frame(co,oc,1,ntap,xpo,ypo,zpo
     &                          ,xpo,ypo,zpo)
                           CALL plot_fragm(fstep,beta,nfragm,fragm,xpo
     &                          ,ypo,zpo,ntap,nres,m1,prsymb,chrge)
                        END IF
                     END IF
                     
*=======================================================================
*----- Write topology --------------------------------------------------
*=======================================================================
                     
                     IF(prttopl) THEN
                        IF(MOD(nstep,ntop_print) .EQ. 0) THEN
                           IF(top_bonds(1) .GT. 0) CALL
     &                          write_bonds(ktopol,fstep,top_bonds,lbnd
     &                          ,lbond,xp0,yp0,zp0)
                           IF(top_bendings(1) .GT. 0) CALL
     &                          write_bends(ktopol,fstep,top_bendings
     &                          ,lbndg,lbend,xp0,yp0,zp0)
                           IF(top_ptors(1) .GT. 0) CALL write_tors('P'
     &                          ,ktopol,fstep,top_ptors,ltor,ltors,xp0
     &                          ,yp0,zp0)
                           IF(top_itors(1) .GT. 0) CALL write_tors('I'
     &                          ,ktopol,fstep,top_itors,litr,litor,xp0
     &                          ,yp0,zp0)
                        END IF
                     END IF
                     
*=======================================================================
*----- Write distance between fragments --------------------------------
*=======================================================================
                     
                     IF(fragm_dist) THEN
                        IF(MOD(nstep,nfragm_dist) .EQ. 0) THEN
                           CALL write_fragm_dist(fstep,kfragm_dist,fragm
     &                          ,nfragm,xpa,ypa,zpa,co)
                        END IF
                     END IF
                     
*=======================================================================
*----- Compute and write gyration ratio --------------------------------
*=======================================================================
                     
                     IF(wrtgyr) THEN
                        CALL write_gyr(kgyr,fstep,protl,nprot,ss_index
     &                       ,xp0,yp0,zp0)
                     END IF
                     
c$$$                  ELSE
c$$$                     nstep=nstep-1
                  END IF
               END IF
            END DO
         END IF
      ELSE IF(dmprnd_o) THEN
         nnstep=nstep_start_o
         DO WHILE(.NOT. end .AND. nstep .LT. stop_anl)
            nstep=nstep+1
            nnstep=nnstep+1
            
*=======================================================================
*----- Read trajectory file by row -------------------------------------
*=======================================================================
               
            IF(start_anl .NE. 0 .AND. dmprnd_i) THEN
               CALL read_confc_rows(co,xp0,yp0,zp0,xau,yau,zau,ntap
     &              ,fstep,nstep,end,err,divide_records,atom_record)

               WRITE(kprint,90000) fstep
               fstep=fstep_o*nnstep
               CALL write_confc(co,xp0,yp0,zp0,ntap,fstep,nnstep
     &              ,nconf,divide_records,atom_record)
               WRITE(kprint,90100) fstep
               WRITE(kprint,*)

            END IF
         END DO

      END IF

      IF(ndipole .GT. 0) THEN

*=======================================================================
*--- Zero correlation arrays -------------------------------------------
*=======================================================================

         CALL zero0(vacf_data,length_fft)

         CALL time_correlation(length_fft,vxpc,vxpc,wsave1,vacf_data
     &        ,ffstep,w1,w2)
         CALL time_correlation(length_fft,vypc,vypc,wsave1,vacf_data
     &        ,ffstep,w1,w2)
         CALL time_correlation(length_fft,vzpc,vzpc,wsave1,vacf_data
     &        ,ffstep,w1,w2)
         WRITE(kdipole,'('' Correlation Solute-Solute'')')
         WRITE(kdipole,'(f10.2,e16.6)') (i*ffstep,vacf_data(i)
     &        /DFLOAT(nstep-i),i=0,nstep-1)
         CALL zero0(vacf_data,length_fft)

         CALL time_correlation(length_fft,vxpc,vxpd,wsave1,vacf_data
     &        ,ffstep,w1,w2)
         CALL time_correlation(length_fft,vypc,vypd,wsave1,vacf_data
     &        ,ffstep,w1,w2)
         CALL time_correlation(length_fft,vzpc,vzpd,wsave1,vacf_data
     &        ,ffstep,w1,w2)

         WRITE(kdipole,'('' Correlation Solute-Solvent'')')
         WRITE(kdipole,'(f10.2,e16.6)') (i*ffstep,2.0D0*(vacf_data(i)
     &        /DFLOAT(nstep-i)),i=0,nstep-1)
         CALL zero0(vacf_data,length_fft)

         CALL time_correlation(length_fft,vxpd,vxpd,wsave1,vacf_data
     &        ,ffstep,w1,w2)
         CALL time_correlation(length_fft,vypd,vypd,wsave1,vacf_data
     &        ,ffstep,w1,w2)
         CALL time_correlation(length_fft,vzpd,vzpd,wsave1,vacf_data
     &        ,ffstep,w1,w2)

         WRITE(kdipole,'('' Correlation Solvent-Solvent'')')
         WRITE(kdipole,'(f10.2,e16.6)') (i*ffstep,vacf_data(i)
     &        /DFLOAT(nstep-i),i=0,nstep-1)
      END IF

      CALL timer(vfcp,tfcp,elapse)
      gcpu=-gcpu + tfcp

*=======================================================================
*     Write timing
*=======================================================================

      WRITE(kprint,*)
      WRITE(kprint,60030)
      WRITE(kprint,17000) gcpu
      IF(time_corr) THEN
         WRITE(kprint,60300) gcpu/DFLOAT(ntap)
      ELSE
         WRITE(kprint,60200) gcpu/DFLOAT(nnstep)
      END IF
      IF(polar) THEN
         WRITE(kprint,60400) polar_counter/DFLOAT(nnstep)
      END IF
      WRITE(kprint,60030)

*================= END OF EXECUTABLE STATEMENTS ========================

1000  FORMAT(/ / / /20('*'),'  M o l e c u l a r   T o p o l o g y  ',
     x     20('*')/ / / / /)
1100  FORMAT('*',13(' '),
     x     ' I n p u t   O p e r a t i o n s   C o m p l e t e d ',
     x     12(' '),'*'/'*',78(' '),'*')
2000  FORMAT(/ /72('*')/'*',70(' '),'*'/
     &     '*     T R A J E C T O R Y   files contains  ',i6,
     &     '   configurations    *'/'*',70(' '),'*'/72('*')/)
2500  FORMAT(/ /72('*')/'*',70(' '),'*'/
     &     '*       C O M P U T I N G   velocity autocorrelation',
     &     ' function          *'/
     &     '*',70(' '),'*'/72('*')/ /)
2501  FORMAT(/ /72('*')/'*',70(' '),'*'/
     &     '*             C O M P U T I N G   atomic rms ',
     &     'displacements             *'/
     &     '*',70(' '),'*'/72('*')/ /)
3000  FORMAT(
     &     '                    ================================'/
     &     '                    =                              ='/
     &     '             -----> =  Processing atoms no. ',i5,'  ='/
     &     '                    =           through no. ',i5,'  ='/
     &     '                    =                              ='/
     &     '                    ================================'/ /)
3500  FORMAT(/ /'Writing v.a.c.f to file  --------> '/ /)
4500  FORMAT(/ /'Writing rms displacement to file  --------> '/ /)
5000  FORMAT(
     &     12x,'================================================='/
     &     12x,'=                                               ='/
     &     12x,'=                                               ='/
     &     12x,'=        SCANNING trajectory file ....          ='/
     &     12x,'=                                               ='/
     &     12x,'=                                               ='/
     &     12x,'================================================='/
     &     )
6000  FORMAT(10x,'<------- Writing Cavities Distribution  --------> ')
60030 FORMAT(/10x,'==========================================='/)
1200  FORMAT(80('*'))
1300  FORMAT('*',78(' '),'*')
17000 FORMAT(  10x,' Total CPU time for data processing  = ',f10.3)
60200 FORMAT(  10x,' Averaged CPU time per step          = ',3x,f7.3)
60300 FORMAT(  10x,' Averaged CPU time per atomic corr.  = ',3x,f7.3)
60400 FORMAT(  10x,' Averaged (polar)iteration per step  = ',3x,f7.3)
90000 FORMAT(' Tstep ', f11.2,' fs successfully processed.'
     &     ,' Program continues... ')
90100 FORMAT(' Tstep ', f11.2,' fs successfully written.'
     &     ,' Program continues... ')
90001 FORMAT(' PDB Configuration successfully processed.')
50000 FORMAT(12x,' Tstep =',f10.1)

80000 FORMAT('REMARK   Rigid body fit on CA atoms')
80100 FORMAT('REMARK   Rigid body fit on heavy atoms')
80101 FORMAT('REMARK   Rigid body fit on backbone atoms')
78410 FORMAT
     &     (/ /' *******ERROR: MAXT for PME too small. INCREASE MAXT.'/
     &     /)
90200 FORMAT('Template and current configuration are at ',f10.3
     &     ,' Angs distance')

93000 FORMAT('Compute Voronoi Volumes--->...')
      RETURN
      END
