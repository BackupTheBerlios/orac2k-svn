      SUBROUTINE read_properties(fxrms,favg,favg_xrms,fvi,gofr_fprint
     &     ,fvoronoi,fcavities,ftop_print,gofr_favg,gofr_fcomp,fprtvaf
     &     ,ftotvaf,fnovaf,fdipole,fnative,ffragm_dist,fhbonds,fhhisto
     &     ,frms,fgyr,freq_ef,freq_dp,fxslt,finst_fit,fcalc_cofm
     &     ,finst_lda,flda_flu,flda_hyd,fprot_hyd,fprot_lda
     &     ,sofk_fprint,sofk_fcomp,err_open,err_args,err_unr
     &     ,err_end)

************************************************************************
*   Time-stamp: <99/09/27 19:04:27 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Nov 18 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      USE HYDRATION_Mod, ONLY: HYD_coeff=>coeff, HYD_cutoff=>cutoff_max,
     &     HYD_Initialize=>Initialize, hydration, HYD_n_neighbors
     &     =>n_neighbors,HYD_ncx=>ncx,HYD_ncy=>ncy,HYD_ncz=>ncz,
     &     HYD_n_write=>n_write
      USE RMS_Matrix_Mod, ONLY: krms_matrix, rms_matrix, Write_Freq,
     &     rms_matrix_avg,RMS_Navg=>Navg,rms_matrix_plot
     &     ,krms_matrix_plot
      USE EUL_Mod, ONLY: EUL_Create=> Create
     &     ,eul_angles,EUL_index_l=>index_l, EUL_add=>Add, EUL_kdump
     &     =>kdump
      USE DENSITY_Mod, ONLY: DEN_input=>Read_it
      USE PDBs_Mod, ONLY: PDB_input=>Read_it

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER iret
      CHARACTER*80 errmsg
      REAL*8  fxrms,fvi,gofr_fprint,gofr_favg,gofr_fcomp,fprtvaf,ftotvaf
     &     ,fdipole,fnative,fnovaf,favg,favg_xrms,fvoronoi,fcavities
     &     ,ftop_print,ffragm_dist,fhbonds,fhhisto,frms,fgyr,sofk_fprint
     &     ,sofk_fcomp
      REAL*8 freq_ef,freq_dp,finst_fit,fcalc_cofm,finst_lda
     &     ,fxslt,flda_flu,flda_hyd,fprot_hyd,fprot_lda,dummy
      CHARACTER*22 err_open
      CHARACTER*37 err_args(*)
      CHARACTER*20 err_end 
      CHARACTER*27 err_unr(4)

*----------------------- VARIABLES IN COMMON --------------------------*
      
      INCLUDE 'parst.h'
      INCLUDE 'cpropar.h'
      INCLUDE 'unit.h'

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER nword,i,j,m,nsevere,nwarning,n_res,n_bonds,NSec
     &     ,n_bendings,n_ptors,n_itors,n_atoms,naux1,naux2,pmol_efp
     &     ,n_st_dummy,n_sv_dummy,read_err
      INTEGER, DIMENSION(:), ALLOCATABLE :: st_dummy,sv_dummy
      CHARACTER*80 line,strngs(40),lined
      CHARACTER*8 fmt
      CHARACTER*1 sep(2),comm(2)
      LOGICAL  exist,lprint
      DATA sep/' ',','/comm/'(',')'/

*----------------------- EXECUTABLE STATEMENTS ------------------------*


c=======================================================================
c     Environment parser starts here 
c=======================================================================

      nfragm = 0
      nsevere=0
      nwarning=0
      n_res=0
      n_bonds=0
      n_bendings=0
      n_ptors=0
      n_itors=0
      n_atoms=0
      n_st_dummy=0
      n_sv_dummy=0
      NSecStruct=1
      SecStructTotal=0

      DO i=1,80
         lined(i:i)=' '
      END DO
      line(79:80)='  '
 100  READ(knlist,'(a78)',END=600) line(1:78)
      CALL wrenc(kprint,line)
      IF(line(1:1) .EQ. '#') GOTO 100 
      CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
      IF(iret.EQ.1) THEN 
         errmsg='while parsing line: toomany strings'
         CALL xerror(errmsg,80,1,2)
         nsevere = nsevere + 1
         go to 100
      END IF

c==== Command  TIME_CORRELATIONS =====================================

      IF(strngs(1).EQ. 'TIME_CORRELATIONS') THEN
         time_corr=.TRUE.
1200     READ(knlist,'(a78)',END=600) line(1:78)
         CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 1200
         CALL parse(line,sep,2,comm,strngs,40,nword,
     x        iret,errmsg)

c----- subcommand vacf

         IF(strngs(1) .EQ. 'vacf' ) THEN
            vacf=.TRUE.
            IF(strngs(2) .EQ. 'OPEN') THEN
               CALL uscrpl(strngs(3),80)
               INQUIRE(FILE=strngs(3),EXIST=exist)
               IF(exist) THEN
                  CALL openf(kvaf,strngs(3),'FORMATTED','OLD',0)
               ELSE
                  CALL openf(kvaf,strngs(3),'FORMATTED','NEW',0)
               END IF
            ELSE
               errmsg=err_open
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF
         
c---- subcommand divide_step
c---  provide number of intepolated points between data points

         ELSE IF(strngs(1).EQ. 'divide_step') THEN
            IF(nword .NE. 2) THEN
               errmsg=err_args(3)//'1'
               CALL xerror(errmsg,80,1,30)
               nsevere=nsevere+1
            ELSE
               CALL fndfmt(1,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) divide_spline
            END IF

c----- subcommand diffusion

         ELSE IF(strngs(1).EQ. 'diffusion') THEN
            diffusion=.TRUE.
            IF(strngs(2) .EQ. 'OPEN') THEN
               CALL uscrpl(strngs(3),80)
               INQUIRE(FILE=strngs(3),EXIST=exist)
               IF(exist) THEN
                  CALL openf(kdiff,strngs(3),'FORMATTED','OLD',0)
               ELSE
                  CALL openf(kdiff,strngs(3),'FORMATTED','NEW',0)
               END IF
            ELSE
               errmsg=err_open
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. 'atoms') THEN
            CALL parse_numbers(err_unr,strngs,nword,corr_atoms
     &           ,n_atoms,nsevere)

         ELSE IF(strngs(1) .EQ. ' ') THEN
            GOTO 1200

         ELSE IF(strngs(1).EQ. 'END' ) THEN
            GOTO 100
         ELSE
            errmsg=err_unr(3)//strngs(2)//err_end(1:14)//err_end(16:20)
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         GOTO 1200

c==== Command  PRINT_DIPOLE===========================================

      ELSE IF(strngs(1).EQ. 'PRINT_DIPOLE') THEN
         not_time_corr=.TRUE.
         IF(strngs(3) .EQ. 'OPEN') THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) fdipole
            CALL uscrpl(strngs(4),80)
            INQUIRE(FILE=strngs(4),EXIST=exist)
            IF(exist) THEN
               CALL openf(kdipole,strngs(4),'FORMATTED','OLD',0)
            ELSE
               CALL openf(kdipole,strngs(4),'FORMATTED','NEW',0)
            END IF
         ELSE
            errmsg=err_unr(3)//strngs(3)
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF

c==== Command  NATIVE ================================================

      ELSE IF(strngs(1).EQ. 'NATIVE') THEN
         not_time_corr=.TRUE.
1500     READ(knlist,'(a78)',END=600) line(1:78)
         CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 1500
         CALL parse(line,sep,2,comm,strngs,40,nword,
     x        iret,errmsg)

                  
         IF(strngs(1) .EQ. 'native' ) THEN
            native=.TRUE.
            IF(strngs(3) .EQ. 'OPEN') THEN
               CALL fndfmt(2,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) fnative
               CALL uscrpl(strngs(2),80)
               INQUIRE(FILE=strngs(4),EXIST=exist)
               IF(exist) THEN
                  CALL openf(knative,strngs(4),'FORMATTED','OLD',0)
               ELSE
                  CALL openf(knative,strngs(4),'FORMATTED','NEW',0)
               END IF
            ELSE
               errmsg=err_unr(3)//strngs(3)
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

c==== Subcommand  cutoff ============================================

         ELSE IF(strngs(1).EQ. 'cutoff') THEN
            not_time_corr=.TRUE.
            IF(nword .NE. 2) THEN
               errmsg=err_args(1)//'1'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            ELSE
               CALL fndfmt(2,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) native_dist
               not_time_corr=.TRUE.
               native=.TRUE.
            END IF

c==== Subcommand  theta ============================================

         ELSE IF(strngs(1).EQ. 'theta') THEN
            not_time_corr=.TRUE.
            IF(nword .NE. 2) THEN
               errmsg=err_args(1)//'1'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            ELSE
               CALL fndfmt(2,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) native_theta
            END IF

c==== Subcommand  check_native =====================================

         ELSE IF(strngs(1).EQ. 'check_native') THEN
            not_time_corr=.TRUE.
            check_native=.TRUE.
            native=.FALSE.
            IF(strngs(3) .EQ. 'OPEN') THEN
               CALL fndfmt(2,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) fnative
               CALL uscrpl(strngs(2),80)
               INQUIRE(FILE=strngs(4),EXIST=exist)
               IF(exist) THEN
                  CALL openf(knative,strngs(4),'FORMATTED','OLD',0)
               ELSE
                  CALL openf(knative,strngs(4),'FORMATTED','NEW',0)
               END IF
            ELSE
               errmsg=err_unr(3)//strngs(3)
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF


c==== Subcommand  check_native =====================================

         ELSE IF(strngs(1).EQ. 'template') THEN
            not_time_corr=.TRUE.
            check_native=.TRUE.
            native=.FALSE.
            IF(nword .EQ. 2) THEN
               INQUIRE(FILE=strngs(2),EXIST=exist)
               IF(exist) THEN
                  CALL openf(knative_tpl,strngs(2),'FORMATTED','OLD',0)
               ELSE
                  errmsg=
     &                 'Native template File does not exist.'
                  CALL xerror(errmsg,80,1,30)
                  nsevere=nsevere+1
               END IF
            ELSE
               errmsg=err_args(1)//'1'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. ' ') THEN
            GOTO 1500

         ELSE IF(strngs(1).EQ. 'END' ) THEN
            GOTO 100
         ELSE
            errmsg=err_unr(3)//strngs(2)//err_end(1:14)//err_end(16:20)
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         GOTO 1500

c==== Command  DIST_FRAGMENT =====================================

      ELSE IF(strngs(1).EQ. 'DIST_FRAGMENT' ) THEN
         not_time_corr=.TRUE.
         fragm_dist=.TRUE.
         IF(strngs(3) .EQ. 'OPEN') THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) ffragm_dist
            CALL uscrpl(strngs(4),80)
            INQUIRE(FILE=strngs(4),EXIST=exist)
            IF(exist) THEN
               CALL openf(kfragm_dist,strngs(4),'FORMATTED','OLD',0)
            ELSE
               CALL openf(kfragm_dist,strngs(4),'FORMATTED','NEW',0)
            END IF
         ELSE
            errmsg=err_unr(3)//strngs(3)
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF

c==== Command  DEF_SECONDARY =============================================

      ELSE IF(strngs(1).EQ. 'DEF_SECONDARY' ) THEN
         SecStructure=.TRUE.
         Nsec=0
         CALL parse_numbers(err_unr,strngs,nword,SecPointer(NSecStruct)
     &        ,NSec,nsevere)
         SecPointer(NSecStruct)=NSec
         NSecStruct=NSecStruct+NSec+1
         SecStructTotal=SecStructTotal+1

c==== Command  DEF_FRAGMMENT =========================================

      ELSE IF(strngs(1).EQ. 'DEF_FRAGMENT' ) THEN
         not_time_corr=.TRUE.
         nfragm = nfragm+1
         IF(nfragm.le.20) THEN
            IF(nword .NE. 3) THEN
               errmsg=err_args(1)//'2'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            ELSE
               CALL fndfmt(1,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) fragm(1,nfragm) 
               CALL fndfmt(1,strngs(3),fmt)
               READ(strngs(3),fmt,err=20) fragm(2,nfragm) 
            END IF
         ELSE
            errmsg='Toomany fragments defined; Max is 20'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF

c==== Command  GOFR===================================================

      ELSE IF(strngs(1) .EQ. 'GOFR') THEN
         not_time_corr=.TRUE.
         gofr=.TRUE.
 700     READ(knlist,'(a78)',END=600) line(1:78)
         CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 700
         CALL parse(line,sep,2,comm,strngs,40,nword,
     x        iret,errmsg)
         IF(strngs(1) .EQ. 'print' ) THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) gofr_fprint
            IF(strngs(3) .EQ. 'OPEN') THEN
               CALL uscrpl(strngs(4),80)
               INQUIRE(FILE=strngs(4),EXIST=exist)
               IF(exist) THEN
                  CALL openf(kgofr_sk,strngs(4),'FORMATTED','OLD',0)
               ELSE
                  CALL openf(kgofr_sk,strngs(4),'FORMATTED','NEW',0)
               END IF
            ELSE
               errmsg='OPEN keyword not found'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. 'delta') THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) delrg

         ELSE IF(strngs(1) .EQ. 'average') THEN
            gofr_avg=.TRUE.
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) gofr_favg

         ELSE IF(strngs(1) .EQ. 'compute') THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) gofr_fcomp

         ELSE IF(strngs(1) .EQ. 'use_neighbor') THEN
            gofr_neighbor=.TRUE.

         ELSE IF(strngs(1) .EQ. 'cutoff') THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) gofr_cut

         ELSE IF(strngs(1) .EQ. 'intra') THEN
            gofr_intra=.TRUE.

         ELSE IF(strngs(1) .EQ. ' ') THEN
            GOTO 700

         ELSE IF(strngs(1).EQ. 'END' ) THEN
            GOTO 100
         ELSE
            errmsg=err_unr(3)//strngs(2)//err_end(1:14)//err_end(16:20)
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         GOTO 700

c==== Command  SOFK ==================================================

      ELSE IF(strngs(1) .EQ. 'SOFK') THEN
         not_time_corr=.TRUE.
         s_of_k=.TRUE.
1400     READ(knlist,'(a78)',END=600) line(1:78)
         CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 1400
         CALL parse(line,sep,2,comm,strngs,40,nword,
     x        iret,errmsg)
         IF(strngs(1) .EQ. 'print' ) THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) sofk_fprint
            IF(strngs(3) .EQ. 'OPEN') THEN
               CALL uscrpl(strngs(4),80)
               INQUIRE(FILE=strngs(4),EXIST=exist)
               IF(exist) THEN
                  CALL openf(ksofk,strngs(4),'FORMATTED','OLD',0)
               ELSE
                  CALL openf(ksofk,strngs(4),'FORMATTED','NEW',0)
               END IF
            ELSE
               errmsg='OPEN keyword not found'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. 'delta') THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) sofk_delta

         ELSE IF(strngs(1) .EQ. 'compute') THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) sofk_fcomp

         ELSE IF(strngs(1) .EQ. 'cutoff') THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) sofk_cut

         ELSE IF(strngs(1) .EQ. ' ') THEN
            GOTO 1400

         ELSE IF(strngs(1).EQ. 'END' ) THEN
            GOTO 100
         ELSE
            errmsg=err_unr(3)//strngs(2)//err_end(1:14)//err_end(16:20)
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         GOTO 1400


c==== Command  HBONDS ================================================

      ELSE IF(strngs(1).EQ. 'HBONDS' ) THEN
         not_time_corr=.TRUE.
1000     READ(knlist,'(a78)',END=600) line(1:78)
         CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 1000
         CALL parse(line,sep,2,comm,strngs,40,nword,
     x        iret,errmsg)

         IF(strngs(1) .EQ. 'total') THEN
            hbonds_tot=.TRUE.

         ELSE IF(strngs(1) .EQ. 'residues') THEN
            hbonds_res=.TRUE.

         ELSE IF(strngs(1) .EQ. 'radial_cutoff') THEN
            IF(nword .EQ. 2) THEN
               CALL fndfmt(2,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) rcut_hb
            ELSE
               errmsg=err_args(2)//'1'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. 'angular_cutoff') THEN
            IF(nword .EQ. 3) THEN
               CALL fndfmt(2,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) acut_hb
               CALL fndfmt(2,strngs(3),fmt)
               READ(strngs(3),fmt,err=20) a2cut_hb
            ELSE IF(nword .EQ. 2) THEN
               CALL fndfmt(2,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) acut_hb
            ELSE
               errmsg=err_args(1)//'1'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. 'print_hyd' ) THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) fhbonds
            IF(strngs(3) .EQ. 'OPEN') THEN
               CALL uscrpl(strngs(4),80)
               INQUIRE(FILE=strngs(4),EXIST=exist)
               IF(exist) THEN
                  CALL openf(khbonds,strngs(4),'FORMATTED','OLD',0)
               ELSE
                  CALL openf(khbonds,strngs(4),'FORMATTED','NEW',0)
               END IF
            ELSE
               errmsg='OPEN keyword not found'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. 'histogram' ) THEN
            hhisto=.TRUE.
            IF(nword .EQ. 2) THEN
               CALL fndfmt(2,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) hhisto_bin
            END IF

         ELSE IF(strngs(1) .EQ. 'print_histo' ) THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) fhhisto
            IF(strngs(3) .EQ. 'OPEN') THEN
               CALL uscrpl(strngs(4),80)
               INQUIRE(FILE=strngs(4),EXIST=exist)
               IF(exist) THEN
                  CALL openf(khhisto,strngs(4),'FORMATTED','OLD',0)
               ELSE
                  CALL openf(khhisto,strngs(4),'FORMATTED','NEW',0)
               END IF
            ELSE
               errmsg='OPEN keyword not found'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. 'use_neighbors' ) THEN
            IF(nword .EQ. 3) THEN
               CALL fndfmt(1,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) update_anl
               CALL fndfmt(2,strngs(3),fmt)
               READ(strngs(3),fmt,err=20) hrcut_update
            ELSE
               nsevere = nsevere + 1
               errmsg=err_args(1) //'2'
               CALL xerror(errmsg,80,1,30)
            END IF

         ELSE IF(strngs(1) .EQ. ' ') THEN
            CONTINUE

         ELSE IF(strngs(1).EQ. 'END' ) THEN
            GOTO 100
         ELSE
            errmsg=err_unr(2)//strngs(2)//err_end(1:14)//err_end(16:20)
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         GOTO 1000

c==== Command  VORONOI==================================================

      ELSE IF(strngs(1) .EQ. 'VORONOI' ) THEN
         not_time_corr=.TRUE.
         voronoi=.TRUE.
         heavy_vor=.FALSE.
800      READ(knlist,'(a78)',END=600) line(1:78)
         CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 800
         CALL parse(line,sep,2,comm,strngs,40,nword,
     x        iret,errmsg)

         IF(strngs(1) .EQ. 'cutoff' ) THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt) cutoff_vor

         ELSE IF(strngs(1) .EQ. 'heavy_atoms') THEN
            heavy_vor=.TRUE.

         ELSE IF(strngs(1) .EQ. 'compute') THEN
            IF(strngs(2) .EQ. 'contact_solute') THEN
               IF(nword .NE. 4) THEN
                  errmsg=err_args(1)//'3'
                  CALL xerror(errmsg,80,1,30)
                  nsevere=nsevere+1
               END IF
               ncontact_slt=ncontact_slt+1
               IF(ncontact_slt .GT. pcontact_slt) THEN
                  errmsg=
     &  ' Number of solute molecules for contact surface exceeded.'
                  CALL xerror(errmsg,80,1,30)
                  nsevere=nsevere+1
               END IF
               voronoi_contact=.TRUE.
               CALL fndfmt(1,strngs(3),fmt)
               READ(strngs(3),fmt,err=20) contact_slt(1,ncontact_slt)
               CALL fndfmt(1,strngs(4),fmt)
               READ(strngs(4),fmt,err=20) contact_slt(2,ncontact_slt)

            ELSE IF(strngs(2) .EQ. 'accessibility') THEN
               voronoi_access=.TRUE.

            ELSE IF(strngs(2) .EQ. 'volume') THEN
               voronoi_volume=.TRUE.

            ELSE IF(strngs(2) .EQ. 'neighbors') THEN
               voronoi_neighbor=.TRUE.

            ELSE
               errmsg=err_unr(2)//strngs(2)
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. 'residue') THEN
            CALL parse_numbers(err_unr,strngs,nword,voronoi_res,n_res
     &           ,nsevere)

         ELSE IF(strngs(1) .EQ. 'print' ) THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) fvoronoi
            IF(strngs(3) .EQ. 'OPEN') THEN
               CALL uscrpl(strngs(4),80)
               INQUIRE(FILE=strngs(4),EXIST=exist)
               IF(exist) THEN
                  CALL openf(kvoronoi,strngs(4),'FORMATTED','OLD',0)
               ELSE
                  CALL openf(kvoronoi,strngs(4),'FORMATTED','NEW',0)
               END IF
            ELSE
               errmsg='OPEN keyword not found'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. ' ') THEN
            GOTO 800

         ELSE IF(strngs(1).EQ. 'END' ) THEN
            GOTO 100
         ELSE
            errmsg=err_unr(2)//strngs(2)//err_end(1:14)//err_end(16:20)
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         GOTO 800

c==== Command  CAVITIES ================================================

      ELSE IF(strngs(1) .EQ. 'CAVITIES' ) THEN
         cavities=.TRUE.
1300     READ(knlist,'(a78)',END=600) line(1:78)
         CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 1300
         CALL parse(line,sep,2,comm,strngs,40,nword,
     x        iret,errmsg)

         IF(strngs(1) .EQ. 'grid' ) THEN
            CALL fndfmt(1,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) nx_cav
            CALL fndfmt(1,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) ny_cav
            CALL fndfmt(1,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) nz_cav

         ELSE IF(strngs(1) .EQ. 'size_atom') THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) size_atom_cav

         ELSE IF(strngs(1) .EQ. 'histo_bin') THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) bin_size_cav

         ELSE IF(strngs(1) .EQ. 'histo_rmax') THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) rmax_size_cav

         ELSE IF(strngs(1) .EQ. 'print' ) THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) fcavities

            IF(strngs(3) .EQ. 'OPEN') THEN
               CALL uscrpl(strngs(4),80)
               cavities_file=strngs(4)
               INQUIRE(FILE=strngs(4),EXIST=exist)
               IF(exist) THEN
                  CALL openf(kcavities,strngs(4),'FORMATTED','OLD',0)
               ELSE
                  CALL openf(kcavities,strngs(4),'FORMATTED','NEW',0)
               END IF
            ELSE
               errmsg='OPEN keyword not found'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. ' ') THEN
            GOTO 1300

         ELSE IF(strngs(1).EQ. 'END' ) THEN
            GOTO 100
         ELSE
            errmsg=err_unr(2)//strngs(2)//err_end(1:14)//err_end(16:20)
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         GOTO 1300

c==== Command FIELD_COFACTOR==========================================

      ELSE IF(strngs(1).EQ. 'FIELD_COFACTOR') THEN 
         not_time_corr=.TRUE.
         IF(nword .LT. 2) THEN
            errmsg=' Ivalid synatx: FIELD_COFACTOR r1 [r2 r3 [OPEN fn]]'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF   
         IF(strngs(5) .EQ. 'OPEN') THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) fvi
            CALL fndfmt(2,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) cut_field
            CALL fndfmt(2,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) alphaf
            CALL uscrpl(strngs(6),80)
            INQUIRE(FILE=strngs(6),EXIST=exist)
            IF(exist) THEN
               CALL openf(kvi,strngs(6),'UNFORMATTED',
     &              'OLD',0)
            ELSE
               CALL openf(kvi,strngs(6),'UNFORMATTED'
     &              ,'NEW',0)
            END IF
         ELSE IF(strngs(3) .EQ. 'OPEN') THEN 
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) fvi
            cut_field=-10.d0
            alphaf=-10.d0
            INQUIRE(FILE=strngs(4),EXIST=exist)
            IF(exist) THEN
               CALL openf(kvi,strngs(4),'UNFORMATTED',
     &              'OLD',0)
            ELSE
               CALL openf(kvi,strngs(4),'UNFORMATTED'
     &              ,'NEW',0)
            END IF
         ELSE   
            errmsg=err_open
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF

c==== Command HYDRATION  ==============================================
                                                                       
         ELSE IF(strngs(1).EQ. 'HYDRATION' ) THEN
            hydration=.TRUE.
            strngs(2)='HYDRATION_FILE'
            INQUIRE(FILE=strngs(2),EXIST=exist)
            IF(exist) THEN
               CALL openf(khydration,strngs(2),'FORMATTED',
     &              'OLD',0)
            ELSE
               CALL openf(khydration,strngs(2),'FORMATTED'
     &              ,'NEW',0)
            END IF
            not_time_corr=.TRUE.
            ALLOCATE(st_dummy(nores+1),sv_dummy(nores+1))
4100        READ(knlist,'(a78)',END=600) line(1:78)                    
            CALL wrenc(kprint,line)
            IF(line(1:1) .EQ. '#') GOTO 4100                           
                                                                       
            CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)    
                                                                       
c==== Subcommand hydration coeff ======================================  

            IF(strngs(1).EQ. 'coeff') THEN
               IF(nword .LT. 2) THEN
                  errmsg=err_args(1)//'3'                              
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere=nsevere+1                                    
               ELSE                                                    
                  CALL fndfmt(2,strngs(2),fmt)                         
                  READ(strngs(2),fmt) HYD_coeff
               END IF                                                  

c==== subcommand cutoff =============================================  

            ELSE IF(strngs(1).EQ. 'cutoff' ) THEN
               IF(nword.ne.1) THEN                                     
                  CALL fndfmt(2,strngs(2),fmt)                         
                  READ(strngs(2),fmt,err=20) HYD_cutoff
               END IF                                                  

c==== subcommand solute =============================================  

            ELSE IF(strngs(1).EQ. 'solute') THEN
               CALL parse_numbers(err_unr,strngs,nword,st_dummy
     &              ,n_st_dummy,nsevere)

c==== subcommand solvent ============================================  

            ELSE IF(strngs(1).EQ. 'solvent') THEN
               CALL parse_numbers(err_unr,strngs,nword,sv_dummy
     &              ,n_sv_dummy,nsevere)

c==== subcommand neighbors ============================================  

            ELSE IF(strngs(1).EQ. 'neighbors') THEN
               IF(nword.ne.1) THEN                                     
                  CALL fndfmt(2,strngs(2),fmt)
                  READ(strngs(2),fmt,err=20) dummy
                  HYD_n_neighbors=IDINT(dummy)
               END IF

c==== subcommand box_grid =============================================  

            ELSE IF(strngs(1).EQ. 'box_grid') THEN
               IF(nword .LT. 4) THEN
                  errmsg=err_args(1)//'3'
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere=nsevere+1
               END IF
               CALL fndfmt(1,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) HYD_ncx
               CALL fndfmt(1,strngs(3),fmt)
               READ(strngs(3),fmt,err=20) HYD_ncy
               CALL fndfmt(1,strngs(4),fmt)
               READ(strngs(4),fmt,err=20) HYD_ncz

c==== subcommand write ==============================================  

            ELSE IF(strngs(1).EQ. 'write' ) THEN
               IF(nword.ne.1) THEN                                     
                  CALL fndfmt(2,strngs(2),fmt)
                  READ(strngs(2),fmt,err=20) dummy
                  HYD_n_write=IDINT(dummy)
               END IF                                                  

c--------------------------------------------------------------------  
                                                                       
            ELSE IF(strngs(1) .EQ. ' ') THEN                           
               GOTO 4100                                               
                                                                       
            ELSE IF(strngs(1).EQ. 'END' ) THEN                         
               GOTO 100                                                
                                                                       
            ELSE                                                       
               errmsg=err_unr(2)//strngs(2)//err_end(1:14)/            
     &              /err_end(16:20)                                    
               CALL xerror(errmsg,80,1,30)                             
               nsevere = nsevere + 1                                   
            END IF                                                     
            GOTO 4100                                                  
                                                                       
                                                                       
c==== END Command HYDRATION ===========================================

         ELSE IF(strngs(1).EQ. 'PDB' ) THEN
            not_time_corr=.TRUE.
            CALL PDB_input(knlist,kprint,nsevere,nword,strngs,iret
     &           ,errmsg,read_err)
            IF(read_err == 1) GOTO 20

c==== Command DENSITY  ================================================
                                                                       
         ELSE IF(strngs(1).EQ. 'DENSITY' ) THEN
            not_time_corr=.TRUE.
            CALL DEN_input(knlist,kprint,nsevere,nword,strngs,iret
     &           ,errmsg,read_err)
            IF(read_err == 1) GOTO 20

c==== END Command DENSITY  ============================================

c==== Command POLARIZATION ============================================
                                                                       
         ELSE IF(strngs(1).EQ. 'POLARIZATION' ) THEN                   
            not_time_corr=.TRUE.                                       
            stretch = .TRUE.                                           
            polar=.TRUE.                                               
2200        READ(knlist,'(a78)',END=600) line(1:78)                    
            CALL wrenc(kprint,line)                                    
            IF(line(1:1) .EQ. '#') GOTO 2200                           
                                                                       
            CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)    
                                                                       
            IF(strngs(1).EQ. 'input_file' ) THEN                       
               IF(nword .LT. 2) THEN                                   
                  errmsg=err_args(1)//'2'                              
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere=nsevere+1                                    
               ELSE                                                    
                  CALL uscrpl(strngs(2),80)                            
                  INQUIRE(FILE=strngs(2),EXIST=exist)                  
                  CALL openf(kpol_inp,strngs(2),'FORMATTED','OLD',0)   
               END IF                                                  
                                                                       
c==== Command  QQ-FUDGE=============================================== 
                                                                       
            ELSE IF(strngs(1) .EQ. 'QQ-FUDGE') THEN                    
               if(nword.eq.2) THEN                                     
                  CALL fndfmt(2,strngs(2),fmt)                         
                  READ(strngs(2),fmt,err=20) fudge                     
               ELSE                                                    
                  errmsg=err_args(1)//'1'                              
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere = nsevere + 1                                
               END IF                                                  
c==== Command  EXTERNAL ============================================== 
                                                                       
            ELSE IF(strngs(1) .EQ. 'EXTERNAL') THEN                    
               if(nword.eq.4) THEN                                     
                  CALL fndfmt(2,strngs(2),fmt)                         
                  READ(strngs(2),fmt,err=20) Ext_ef(1)                     
                  CALL fndfmt(2,strngs(3),fmt)                         
                  READ(strngs(3),fmt,err=20) Ext_ef(2)                     
                  CALL fndfmt(2,strngs(4),fmt)                         
                  READ(strngs(4),fmt,err=20) Ext_ef(3)                     
               ELSE                                                    
                  errmsg=err_args(1)//'1'                              
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere = nsevere + 1                                
               END IF                                                  
c==== Command  SCALE=============================================== 
                                                                       
            ELSE IF(strngs(1) .EQ. 'SCALE') THEN                    
               if(nword.eq.2) THEN                                     
                  CALL fndfmt(2,strngs(2),fmt)                         
                  READ(strngs(2),fmt,err=20) polar_scale                     
               ELSE                                                    
                  errmsg=err_args(1)//'1'                              
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere = nsevere + 1                                
               END IF                   
c==== Subcommand print dipole =======================================  
                                                                       
            ELSE IF(strngs(1).EQ. 'print_dipole' ) THEN                       
               IF(nword .LT. 4) THEN                                   
                  errmsg=err_args(1)//'3'                              
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere=nsevere+1                                    
               END IF                                                  
               IF(strngs(3) .EQ. 'OPEN') THEN                          
                  print_dp=.TRUE.
                  freq_dp = 1.0d0
                  CALL fndfmt(2,strngs(2),fmt)                         
                  READ(strngs(2),fmt) freq_dp
                  CALL uscrpl(strngs(4),80)                            
                  INQUIRE(FILE=strngs(4),EXIST=exist)                  
                  IF(exist) THEN                                       
                     CALL openf(kout_dp,strngs(4),'FORMATTED',         
     &                    'OLD',0)                                     
                  ELSE                                                 
                     CALL openf(kout_dp,strngs(4),'FORMATTED'          
     &                    ,'NEW',0)                                    
                  END IF                                               
                                                                       
                  CALL add_str(strngs(4),80,'_sta',4,lined)            
                  INQUIRE(FILE=lined,EXIST=exist)                      
                  IF(exist) THEN                                       
                     CALL openf(kout_dp_sta,lined,'FORMATTED',         
     &                    'OLD',0)                                     
                  ELSE                                                 
                     CALL openf(kout_dp_sta,lined,'FORMATTED'          
     &                    ,'NEW',0)                                    
                  END IF                                               
                                                                       
                  CALL add_str(strngs(4),80,'_tot',4,lined)            
                  INQUIRE(FILE=lined,EXIST=exist)                      
                  IF(exist) THEN                                       
                     CALL openf(kout_dp_tot,lined,'FORMATTED',         
     &                    'OLD',0)                                     
                  ELSE                                                 
                     CALL openf(kout_dp_tot,lined,'FORMATTED'          
     &                    ,'NEW',0)                                    
                  END IF                                               
               ELSE                                                    
                  errmsg=err_open                                      
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere = nsevere + 1                                
               END IF                                                  
                                                                       
c--------------------------------------------------------------------  
                                                                       
            ELSE IF(strngs(1) .EQ. ' ') THEN                           
               GOTO 2200                                               
                                                                       
            ELSE IF(strngs(1).EQ. 'END' ) THEN                         
               GOTO 100                                                
                                                                       
            ELSE                                                       
               errmsg=err_unr(2)//strngs(2)//err_end(1:14)/            
     &              /err_end(16:20)                                    
               CALL xerror(errmsg,80,1,30)                             
               nsevere = nsevere + 1                                   
            END IF                                                     
            GOTO 2200                                                  
                                                                       
c==== Command POTENTIAL ===============================================

c$$$===================================================================
c$$$-- Compute electrostatic potential on a grid
c$$$===================================================================

         ELSE IF(strngs(1).EQ. 'E-POTENTIAL' ) THEN
            not_time_corr=.TRUE.                                       
            EPotential=.TRUE.
2300        READ(knlist,'(a78)',END=600) line(1:78)                    
            CALL wrenc(kprint,line)                                    
            IF(line(1:1) .EQ. '#') GOTO 2300

            CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)    
            
c==== Subcommand EWALD ==============================================  
                                                                       
            IF(strngs(1).EQ. 'EWALD') THEN                             
               clewld=.TRUE.                                           
               IF(strngs(2).EQ. 'PME') THEN                            
                  pme =.TRUE.                                          
                  grpcut=.FALSE.                                       
                  IF(nword .EQ. 7) THEN                                
                     CALL fndfmt(2,strngs(3),fmt)                      
                     READ(strngs(3),fmt,err=20) alphal                 
                     CALL fndfmt(1,strngs(4),fmt)                      
                     READ(strngs(4),fmt,err=20) nfft1                  
                     CALL fndfmt(1,strngs(5),fmt)                      
                     READ(strngs(5),fmt,err=20) nfft2                  
                     CALL fndfmt(1,strngs(6),fmt)                      
                     READ(strngs(6),fmt,err=20) nfft3                  
                     CALL fndfmt(1,strngs(7),fmt)                      
                     READ(strngs(7),fmt,err=20) pme_order              
                  ELSE                                                 
                     nsevere = nsevere + 1                             
                     errmsg=err_args(1) // '4 after keyword "PME"'     
                     CALL xerror(errmsg,80,1,30)                       
                  END IF                                               
                                                                       
               ELSE IF(strngs(2).EQ. 'ON' ) THEN                       
                  clewld=.TRUE.                                        
                  IF(nword .EQ. 4) THEN                                
                     CALL fndfmt(2,strngs(3),fmt)                      
                     READ(strngs(3),fmt,err=20) alphal                 
                     CALL fndfmt(2,strngs(4),fmt)                      
                     READ(strngs(4),fmt,err=20) rkcut                  
                  ELSE IF(nword .EQ. 3) THEN                           
                     CALL fndfmt(2,strngs(3),fmt)                      
                     READ(strngs(3),fmt,err=20) alphal                 
                  ELSE                                                 
                     nsevere = nsevere + 1                             
                     errmsg=err_args(1) // '1 after keyword "on"'      
                     CALL xerror(errmsg,80,1,30)                       
                  END IF                                               
               ELSE IF(strngs(2).EQ. 'OFF' ) THEN                      
                  clewld=.FALSE.                                       
               ELSE                                                    
                  nsevere = nsevere + 1                                
                  errmsg=err_unr(3) // strngs(2)                       
                  CALL xerror(errmsg,80,1,30)                          
               END IF                                                  
                                                                       
c==== Command ERFC_SPLINE==============================================
                                                                       
            ELSE IF(strngs(1).EQ. 'ERFC_SPLINE' ) THEN                 
               erfc_spline=.TRUE.                                      
               IF(nword .NE. 1) THEN                                   
                  CALL fndfmt(2,strngs(2),fmt)                         
                  READ(strngs(2),fmt,err=20) erfc_bin                  
               END IF                                                  
                                                                       
c==== subcommand CUTOFF =============================================  
                                                                       
            ELSE IF(strngs(1).EQ. 'CUTOFF' ) THEN                      
               IF(nword.ne.1) THEN                                     
                  CALL fndfmt(2,strngs(2),fmt)                         
                  READ(strngs(2),fmt,err=20) rspoff                    
               END IF                                                  

c==== Command  FREQUENCY=============================================== 
                                                                       
            ELSE IF(strngs(1) .EQ. 'FREQUENCY') THEN                    
               if(nword.eq.2) THEN                                     
                  CALL fndfmt(1,strngs(2),fmt)                         
                  READ(strngs(2),fmt,err=20) NFreqPotential
               ELSE                                                    
                  errmsg=err_args(1)//'1'                              
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere = nsevere + 1                                
               END IF                                                  

c==== Command  ChargeIon=============================================== 
                                                                       
            ELSE IF(strngs(1) .EQ. 'CHARGEION') THEN                    
               if(nword.eq.2) THEN                                     
                  CALL fndfmt(2,strngs(2),fmt)                         
                  READ(strngs(2),fmt,err=20) ChargeIon
               ELSE                                                    
                  errmsg=err_args(1)//'1'                              
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere = nsevere + 1                                
               END IF                                                  
                                                                       
c==== Command  SigmaIon================================================ 
                                                                       
            ELSE IF(strngs(1) .EQ. 'SIGMAION') THEN                    
               if(nword.eq.2) THEN                                     
                  CALL fndfmt(2,strngs(2),fmt)                         
                  READ(strngs(2),fmt,err=20) SigmaIon
               ELSE                                                    
                  errmsg=err_args(1)//'1'                              
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere = nsevere + 1                                
               END IF                                                  

c==== Command  Smoothfactor============================================ 
                                                                       
            ELSE IF(strngs(1) .EQ. 'SMOOTH') THEN                    
               if(nword.eq.2) THEN                                     
                  CALL fndfmt(2,strngs(2),fmt)                         
                  READ(strngs(2),fmt,err=20) SmoothFactor
               ELSE                                                    
                  errmsg=err_args(1)//'1'                              
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere = nsevere + 1                                
               END IF                                                  

c==== Command  DOFREEENERGY============================================ 
                                                                       
            ELSE IF(strngs(1) .EQ. 'FREEENERGY') THEN
               DoFreeEnergy=.TRUE.

c==== Command  SigmaIon================================================ 
                                                                       
            ELSE IF(strngs(1).EQ. 'PRINT' ) THEN
               IF(nword .LT. 3) THEN                                   
                  errmsg=err_args(1)//'2'                              
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere=nsevere+1                                    
               END IF                                                  
               IF(strngs(2) .EQ. 'OPEN') THEN                          
                  print_ef=.TRUE.                                      
                  CALL uscrpl(strngs(3),80)                            
                  INQUIRE(FILE=strngs(3),EXIST=exist)                  
                  IF(exist) THEN                                       
                     CALL openf(KElecPot,strngs(3),'FORMATTED',         
     &                    'OLD',0)                                     
                  ELSE                                                 
                     errmsg='Auxiliary file for Electrostatic Potential'
     &                    //' does not exist. Abort.'
                     CALL xerror(errmsg,80,1,30)
                     nsevere = nsevere + 1
                  END IF
               ELSE                                                    
                  errmsg=err_open                                      
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere = nsevere + 1                                
               END IF                                                  

            END IF
            GOTO 2300
                                                                       
c==== Command ELECTRIC ================================================
                                                                       
         ELSE IF(strngs(1).EQ. 'ELECTRIC' ) THEN
            not_time_corr=.TRUE.                                       
            efield=.TRUE.                                              
2100        READ(knlist,'(a78)',END=600) line(1:78)                    
            CALL wrenc(kprint,line)                                    
            IF(line(1:1) .EQ. '#') GOTO 2100                           
                                                                       
            CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)    
                                                                       
c==== Subcommand EWALD ==============================================  
                                                                       
            IF(strngs(1).EQ. 'EWALD') THEN                             
               clewld=.TRUE.                                           
               IF(strngs(2).EQ. 'PME') THEN                            
                  pme =.TRUE.                                          
                  grpcut=.FALSE.                                       
                  IF(nword .EQ. 7) THEN                                
                     CALL fndfmt(2,strngs(3),fmt)                      
                     READ(strngs(3),fmt,err=20) alphal                 
                     CALL fndfmt(1,strngs(4),fmt)                      
                     READ(strngs(4),fmt,err=20) nfft1                  
                     CALL fndfmt(1,strngs(5),fmt)                      
                     READ(strngs(5),fmt,err=20) nfft2                  
                     CALL fndfmt(1,strngs(6),fmt)                      
                     READ(strngs(6),fmt,err=20) nfft3                  
                     CALL fndfmt(1,strngs(7),fmt)                      
                     READ(strngs(7),fmt,err=20) pme_order              
                  ELSE                                                 
                     nsevere = nsevere + 1                             
                     errmsg=err_args(1) // '4 after keyword "PME"'     
                     CALL xerror(errmsg,80,1,30)                       
                  END IF                                               
                                                                       
               ELSE IF(strngs(2).EQ. 'ON' ) THEN                       
                  clewld=.TRUE.                                        
                  IF(nword .EQ. 4) THEN                                
                     CALL fndfmt(2,strngs(3),fmt)                      
                     READ(strngs(3),fmt,err=20) alphal                 
                     CALL fndfmt(2,strngs(4),fmt)                      
                     READ(strngs(4),fmt,err=20) rkcut                  
                  ELSE IF(nword .EQ. 3) THEN                           
                     CALL fndfmt(2,strngs(3),fmt)                      
                     READ(strngs(3),fmt,err=20) alphal                 
                  ELSE                                                 
                     nsevere = nsevere + 1                             
                     errmsg=err_args(1) // '1 after keyword "on"'      
                     CALL xerror(errmsg,80,1,30)                       
                  END IF                                               
               ELSE IF(strngs(2).EQ. 'OFF' ) THEN                      
                  clewld=.FALSE.                                       
                                                                       
                                                                       
               ELSE                                                    
                  nsevere = nsevere + 1                                
                  errmsg=err_unr(3) // strngs(2)                       
                  CALL xerror(errmsg,80,1,30)                          
               END IF                                                  
                                                                       
c==== Command ERFC_SPLINE==============================================
                                                                       
            ELSE IF(strngs(1).EQ. 'ERFC_SPLINE' ) THEN                 
               erfc_spline=.TRUE.                                      
               IF(nword .NE. 1) THEN                                   
                  CALL fndfmt(2,strngs(2),fmt)                         
                  READ(strngs(2),fmt,err=20) erfc_bin                  
               END IF                                                  
                                                                       
c==== subcommand CUTOFF =============================================  
                                                                       
            ELSE IF(strngs(1).EQ. 'CUTOFF' ) THEN                      
               IF(nword.ne.1) THEN                                     
                  CALL fndfmt(2,strngs(2),fmt)                         
                  READ(strngs(2),fmt,err=20) rspoff                    
               END IF                                                  

c==== Command  FREQUENCY=============================================== 
                                                                       
            ELSE IF(strngs(1) .EQ. 'FREQUENCY') THEN                    
               if(nword.eq.2) THEN                                     
                  CALL fndfmt(1,strngs(2),fmt)                         
                  READ(strngs(2),fmt,err=20) nfreq_polar                     
               ELSE                                                    
                  errmsg=err_args(1)//'1'                              
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere = nsevere + 1                                
               END IF                                                  
                                                                       
c==== Subcommand sel_mol ============================================  
                                                                       
            ELSE IF(strngs(1).EQ. 'sel_mol') THEN                      
               lmol_ef=.TRUE.                                          
               nmol_ef = nmol_ef + 1                                   
               pmol_efp=0                                              
               CALL parse_numbers(err_unr,strngs,nword,pmol_ef(1,nmol_ef
     &              ),pmol_efp,nsevere)                                
               IF(pmol_efp .GT. f1) THEN                               
                  errmsg=' Length of the chromophores arrays '         
     &                 / /'insufficient. Increase _F1_.'               
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere = nsevere + 1                                
               END IF                                                  
               pmol_ef(1,nmol_ef)=pmol_efp                             
                                                                       
c==== Subcommand print efield =======================================  
                                                                       
            ELSE IF(strngs(1).EQ. 'print' ) THEN                       
               IF(nword .LT. 4) THEN                                   
                  errmsg=err_args(1)//'3'                              
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere=nsevere+1                                    
               END IF                                                  
               IF(strngs(3) .EQ. 'OPEN') THEN                          
                  print_ef=.TRUE.                                      
                  CALL fndfmt(2,strngs(2),fmt)                         
                  READ(strngs(2),fmt) freq_ef                          
                  CALL uscrpl(strngs(4),80)                            
                  INQUIRE(FILE=strngs(4),EXIST=exist)                  
                  IF(exist) THEN                                       
                     CALL openf(kout_ef,strngs(4),'FORMATTED',         
     &                    'OLD',0)                                     
                  ELSE                                                 
                     CALL openf(kout_ef,strngs(4),'FORMATTED'          
     &                    ,'NEW',0)                                    
                  END IF                                               
               ELSE                                                    
                  errmsg=err_open                                      
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere = nsevere + 1                                
               END IF                                                  
                                                                       
c--------------------------------------------------------------------  
                                                                       
            ELSE IF(strngs(1) .EQ. ' ') THEN                           
               GOTO 2100                                               
                                                                       
            ELSE IF(strngs(1).EQ. 'END' ) THEN                         
               GOTO 100                                                
                                                                       
            ELSE                                                       
               errmsg=err_unr(2)//strngs(2)//err_end(1:14)/            
     &              /err_end(16:20)                                    
               CALL xerror(errmsg,80,1,30)                             
               nsevere = nsevere + 1                                   
            END IF                                                     
            GOTO 2100                                                  
                                                                       
                                                                       
c==== END Command ELECTRIC_FIELD ======================================
c==== Command LDA =====================================================
                                                                       
         ELSE IF(strngs(1).EQ. 'LDA' ) THEN                            
            not_time_corr=.TRUE.                                       
3100        READ(knlist,'(a78)',END=600) line(1:78)                    
            CALL wrenc(kprint,line)                                    
            IF(line(1:1) .EQ. '#') GOTO 3100                           
                                                                       
            CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)    
                                                                       
c==== Subcommand sel_lda ===========================================   
                                                                       
            IF(strngs(1).EQ. 'select_lda') THEN                        
               IF(nword .LT. 5) THEN                                   
                  errmsg=err_args(1)//'5'                              
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere=nsevere+1                                    
               ELSE                                                    
                  sel_lda=.TRUE.                                       
                  CALL fndfmt(1,strngs(2),fmt)                         
                  READ(strngs(2),fmt) nlda_mol                         
                  CALL fndfmt(1,strngs(3),fmt)                         
                  READ(strngs(3),fmt) nlda_atm                         
                  CALL fndfmt(1,strngs(4),fmt)                         
                  READ(strngs(4),fmt) nlda_zero                        
                  CALL fndfmt(1,strngs(5),fmt)                         
                  READ(strngs(5),fmt) nlda_grp                         
               END IF                                                  
                                                                       
c==== Subcommand fluidity ==========================================   
                                                                       
            ELSE IF(strngs(1).EQ. 'fluidity')  THEN                    
               IF(nword .LT. 4) THEN                                   
                  errmsg=err_args(1)//'4'                              
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere=nsevere+1                                    
               END IF                                                  
               IF(strngs(3) .EQ. 'OPEN') THEN                          
                  lda_flu=.TRUE.                                       
                  CALL fndfmt(2,strngs(2),fmt)                         
                  READ(strngs(2),fmt) flda_flu                         
                  CALL uscrpl(strngs(4),80)                            
                  INQUIRE(FILE=strngs(4),EXIST=exist)                  
                  IF(exist) THEN                                       
                     CALL openf(klda_flu,strngs(4),'FORMATTED',        
     &                    'OLD',0)                                     
                  ELSE                                                 
                     CALL openf(klda_flu,strngs(4),'FORMATTED'         
     &                    ,'NEW',0)                                    
                  END IF                                               
               ELSE                                                    
                  errmsg=err_open                                      
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere = nsevere + 1                                
               END IF                                                  
                                                                       
c==== Subcommand hydration =========================================   
                                                                       
            ELSE IF(strngs(1).EQ. 'hydration')  THEN                   
               IF(nword .LT. 4) THEN                                   
                  errmsg=err_args(1)//'4'                              
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere=nsevere+1                                    
               END IF                                                  
               IF(strngs(3) .EQ. 'OPEN') THEN                          
                  lda_hyd=.TRUE.                                       
                  CALL fndfmt(2,strngs(2),fmt)                         
                  READ(strngs(2),fmt) flda_hyd                         
                  CALL uscrpl(strngs(4),80)                            
                  INQUIRE(FILE=strngs(4),EXIST=exist)                  
                  IF(exist) THEN                                       
                     CALL openf(klda_hyd,strngs(4),'FORMATTED',        
     &                    'OLD',0)                                     
                  ELSE                                                 
                     CALL openf(klda_hyd,strngs(4),'FORMATTED'         
     &                    ,'NEW',0)                                    
                  END IF                                               
               ELSE                                                    
                  errmsg=err_open                                      
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere = nsevere + 1                                
               END IF                                                  
                                                                       
c==== Subcommand print LDA =========================================   
                                                                       
            ELSE IF(strngs(1).EQ. 'print' .AND.                        
     &                          strngs(2) .EQ. 'LDA')  THEN            
               IF(nword .LT. 5) THEN                                   
                  errmsg=err_args(1)//'4'                              
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere=nsevere+1                                    
               END IF                                                  
                                                                       
               IF(strngs(4) .EQ. 'OPEN') THEN                          
                  avg_lda=.TRUE.                                       
                  CALL fndfmt(2,strngs(3),fmt)                         
                  READ(strngs(3),fmt) finst_lda                        
                  CALL uscrpl(strngs(5),80)                            
                                                                       
                  CALL add_str(strngs(5),80,'_rmin',5,lined)           
                  INQUIRE(FILE=lined,EXIST=exist)                      
                  IF(exist) THEN                                       
                     CALL openf(klda_rmin,lined,'FORMATTED',           
     &                    'OLD',0)                                     
                  ELSE                                                 
                     CALL openf(klda_rmin,lined,'FORMATTED'            
     &                    ,'NEW',0)                                    
                  END IF                                               
                                                                       
                  CALL add_str(strngs(5),80,'_inst',5,lined)           
                  INQUIRE(FILE=lined,EXIST=exist)                      
                  IF(exist) THEN                                       
                     CALL openf(klda_inst,lined,'FORMATTED',           
     &                    'OLD',0)                                     
                  ELSE                                                 
                     CALL openf(klda_inst,lined,'FORMATTED'            
     &                    ,'NEW',0)                                    
                  END IF                                               
                                                                       
                  CALL add_str(strngs(5),80,'_eend',5,lined)           
                  INQUIRE(FILE=lined,EXIST=exist)                      
                  IF(exist) THEN                                       
                     CALL openf(klda_eend,lined,'FORMATTED',           
     &                    'OLD',0)                                     
                  ELSE                                                 
                     CALL openf(klda_eend,lined,'FORMATTED'            
     &                    ,'NEW',0)                                    
                  END IF                                               
               ELSE                                                    
                  errmsg=err_open                                      
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere = nsevere + 1                                
               END IF                                                  
c==== Subcommand solvation (protein) ================================  
            ELSE IF(strngs(1).EQ. 'solvation' .AND.                    
     &         strngs(2) .EQ. 'coeff' ) THEN                           
               IF(nword .LT. 3) THEN                                   
                  errmsg=err_args(1)//'3'                              
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere=nsevere+1                                    
               ELSE                                                    
                  CALL fndfmt(2,strngs(3),fmt)                         
                  READ(strngs(3),fmt) coeff_lda                        
               END IF                                                  
c---------------                                                       
            ELSE IF(strngs(1).EQ. 'solvation' .AND.                    
     &         strngs(2) .EQ. 'residues' ) THEN                        
               IF(nword .LT. 4) THEN                                   
                  errmsg=err_args(1)//'4'                              
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere=nsevere+1                                    
               ELSE                                                    
                  CALL fndfmt(1,strngs(3),fmt)                         
                  READ(strngs(3),fmt) min_lda                          
                  CALL fndfmt(1,strngs(4),fmt)                         
                  READ(strngs(4),fmt) max_lda                          
               END IF                                                  
c==== subcommand CUTOFF =============================================  
            ELSE IF(strngs(1).EQ. 'CUTOFF' ) THEN                      
               IF(nword.ne.1) THEN                                     
                  CALL fndfmt(2,strngs(2),fmt)                         
                  READ(strngs(2),fmt,err=20) rspoff                    
               END IF                                                  
c==== Subcommand print solvation ====================================  
            ELSE IF(strngs(1).EQ. 'residence_times') THEN              
                    res_time=.TRUE.                                    
            ELSE IF(strngs(1).EQ. 'print' .AND.                        
     &         strngs(2) .EQ. 'solvation') THEN                        
               IF(nword .LT. 5) THEN                                   
                  errmsg=err_args(1)//'4'                              
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere=nsevere+1                                    
               END IF                                                  
               IF(strngs(4) .EQ. 'OPEN') THEN                          
                  prot_lda=.TRUE.                                      
                  CALL fndfmt(2,strngs(3),fmt)                         
                  READ(strngs(3),fmt) fprot_lda                        
                  CALL uscrpl(strngs(5),80)                            
                  INQUIRE(FILE=strngs(5),EXIST=exist)                  
                  IF(exist) THEN                                       
                     CALL openf(kprot_lda,strngs(5),'FORMATTED',       
     &                    'OLD',0)                                     
                  ELSE                                                 
                     CALL openf(kprot_lda,strngs(5),'FORMATTED'        
     &                    ,'NEW',0)                                    
                  END IF                                               
                                                                       
                  CALL add_str(strngs(5),80,'_rest',5,lined)           
                  INQUIRE(FILE=lined,EXIST=exist)                      
                  IF(exist) THEN                                       
                     CALL openf(klda_rest,lined,'FORMATTED',           
     &                    'OLD',0)                                     
                  ELSE                                                 
                     CALL openf(klda_rest,lined,'FORMATTED'            
     &                    ,'NEW',0)                                    
                  END IF                                               
                                                                       
               ELSE                                                    
                  errmsg=err_open                                      
                  CALL xerror(errmsg,80,1,30)                          
                  nsevere = nsevere + 1                                
               END IF                                                  
c--------------------------------------------------------------------  
                                                                       
            ELSE IF(strngs(1) .EQ. ' ') THEN                           
               GOTO 3100                                               
                                                                       
            ELSE IF(strngs(1).EQ. 'END' ) THEN                         
               GOTO 100                                                
                                                                       
            ELSE                                                       
               errmsg=err_unr(2)//strngs(2)//err_end(1:14)/            
     &              /err_end(16:20)                                    
               CALL xerror(errmsg,80,1,30)                             
               nsevere = nsevere + 1                                   
            END IF                                                     
            GOTO 3100                                                  
                                                                       
c==== END Command LDA =================================================
c==== Command COMPARE =================================================

         ELSE IF(strngs(1).EQ. 'STRUCTURES' ) THEN
            not_time_corr=.TRUE.
1100        READ(knlist,'(a78)',END=600) line(1:78)
            CALL wrenc(kprint,line)
            IF(line(1:1) .EQ. '#') GOTO 1100

            CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)

c==== Subcommand print inst_fit =====================================

            IF(strngs(1).EQ. 'print' .AND. 
     &              strngs(2) .EQ. 'inst_fit') THEN
               IF(nword .LT. 5) THEN
                  errmsg=err_args(1)//'4'
                  CALL xerror(errmsg,80,1,30)
                  nsevere=nsevere+1
               END IF
               IF(strngs(4) .EQ. 'OPEN') THEN
                  inst_fit=.TRUE.
                  CALL fndfmt(2,strngs(3),fmt)
                  READ(strngs(3),fmt) finst_fit
                  CALL uscrpl(strngs(5),80)
                  INQUIRE(FILE=strngs(5),EXIST=exist)
                  IF(exist) THEN
                     CALL openf(kfit,strngs(5),'FORMATTED',
     &                    'OLD',0)
                  ELSE
                     CALL openf(kfit,strngs(5),'FORMATTED'
     &                    ,'NEW',0)
                  END IF
               ELSE
                  errmsg=err_open
                  CALL xerror(errmsg,80,1,30)
                  nsevere = nsevere + 1
               END IF

c==== Subcommand print COFM=========================================

            ELSE IF(strngs(1).EQ. 'print' .AND.
     &                          strngs(2) .EQ. 'COFM')  THEN
               IF(nword .LT. 5) THEN
                  errmsg=err_args(1)//'4'
                  CALL xerror(errmsg,80,1,30)
                  nsevere=nsevere+1
               END IF
               IF(strngs(4) .EQ. 'OPEN') THEN
                  calc_cofm=.TRUE.
                  CALL fndfmt(2,strngs(3),fmt)
                  READ(strngs(3),fmt) fcalc_cofm
                  CALL uscrpl(strngs(5),80)
                  INQUIRE(FILE=strngs(5),EXIST=exist)
                  IF(exist) THEN

                     CALL openf(kcalc_cofm,strngs(5),'FORMATTED',
     &                    'OLD',0)
                  ELSE
                     CALL openf(kcalc_cofm,strngs(5),'FORMATTED'
     &                    ,'NEW',0)
                  END IF
               ELSE
                  errmsg=err_open
                  CALL xerror(errmsg,80,1,30)
                  nsevere = nsevere + 1
               END IF
            
c==== Subcommand print structure ====================================

            ELSE IF(strngs(1).EQ. 'print' .AND. 
     &              strngs(2) .EQ. 'averaged') THEN
               IF(nword .LT. 5) THEN
                  errmsg=err_args(1)//'4'
                  CALL xerror(errmsg,80,1,30)
                  nsevere=nsevere+1
               END IF
               IF(strngs(4) .EQ. 'OPEN') THEN
                  avg_str=.TRUE.
                  CALL fndfmt(2,strngs(3),fmt)
                  READ(strngs(3),fmt) favg
                  CALL uscrpl(strngs(5),80)
                  INQUIRE(FILE=strngs(5),EXIST=exist)
                  IF(exist) THEN
                     CALL openf(kavg,strngs(5),'FORMATTED',
     &                    'OLD',0)
                  ELSE
                     CALL openf(kavg,strngs(5),'FORMATTED'
     &                    ,'NEW',0)
                  END IF
               ELSE
                  errmsg=err_open
                  CALL xerror(errmsg,80,1,30)
                  nsevere = nsevere + 1
               END IF
            
c==== Subcommand print avg_xrms =====================================

            ELSE IF(strngs(1).EQ. 'print' .AND. strngs(2) .EQ.
     &              'avg_xrms') THEN
               IF(nword .LT. 5) THEN
                  errmsg=err_args(1)//'4'
                  CALL xerror(errmsg,80,1,30)
                  nsevere=nsevere+1
               END IF
               IF(strngs(4) .EQ. 'OPEN') THEN
                  avg_str=.TRUE.
                  CALL fndfmt(2,strngs(3),fmt)
                  READ(strngs(3),fmt) favg_xrms
                  CALL uscrpl(strngs(5),80)
                  INQUIRE(FILE=strngs(5),EXIST=exist)
                  IF(exist) THEN
                     CALL openf(kavg_xrms,strngs(5),'FORMATTED',
     &                    'OLD',0)
                  ELSE
                     CALL openf(kavg_xrms,strngs(5),'FORMATTED'
     &                    ,'NEW',0)
                  END IF
               ELSE
                  errmsg=err_open
                  CALL xerror(errmsg,80,1,30)
                  nsevere = nsevere + 1
               END IF

c==== Subcommand print inst_xslt ====================================
 
            ELSE IF(strngs(1).EQ. 'print' .AND. strngs(2) .EQ.
     &              'inst_xslt') THEN
               IF(nword .LT. 5) THEN
                  errmsg=err_args(1)//'4'
                  CALL xerror(errmsg,80,1,30)
                  nsevere=nsevere+1
               END IF
               IF(strngs(4) .EQ. 'OPEN') THEN
                  CALL fndfmt(2,strngs(3),fmt)
                  READ(strngs(3),fmt,err=20) fxslt
                  CALL uscrpl(strngs(5),80)
                  INQUIRE(FILE=strngs(5),EXIST=exist)
                  IF(exist) THEN
                     CALL openf(kxslt,strngs(5),'FORMATTED',
     &                    'OLD',0)
                  ELSE
                     CALL openf(kxslt,strngs(5),'FORMATTED'
     &                    ,'NEW',0)
                  END IF
               ELSE
                  errmsg=err_open
                  CALL xerror(errmsg,80,1,30)
                  nsevere = nsevere + 1
               END IF

c==== Subcommand print inst_xrms ====================================

            ELSE IF(strngs(1).EQ. 'print' .AND. strngs(2) .EQ.
     &              'inst_xrms') THEN
               IF(nword .LT. 5) THEN
                  errmsg=err_args(1)//'4'
                  CALL xerror(errmsg,80,1,30)
                  nsevere=nsevere+1
               END IF
               IF(strngs(4) .EQ. 'OPEN') THEN
                  CALL fndfmt(2,strngs(3),fmt)
                  READ(strngs(3),fmt,err=20) fxrms
                  CALL uscrpl(strngs(5),80)
                  INQUIRE(FILE=strngs(5),EXIST=exist)
                  IF(exist) THEN
                     CALL openf(kxrms,strngs(5),'FORMATTED',
     &                    'OLD',0)
                  ELSE
                     CALL openf(kxrms,strngs(5),'FORMATTED'
     &                    ,'NEW',0)
                  END IF
                  CALL add_str(strngs(5),80,'_atm',4,lined)
                  INQUIRE(FILE=lined,EXIST=exist)
                  IF(exist) THEN
                     CALL openf(kxrms_atm,lined,'FORMATTED',
     &                    'OLD',0)
                  ELSE
                     CALL openf(kxrms_atm,lined,'FORMATTED'
     &                    ,'NEW',0)
                  END IF
               ELSE
                  errmsg=err_open
                  CALL xerror(errmsg,80,1,30)
                  nsevere = nsevere + 1
               END IF

c==== Subcommand ins_xrms ============================================

            ELSE IF(strngs(1).EQ. 'inst_xrms') THEN
               anxrms=.TRUE.
               DO i=2,nword
                  IF(strngs(i).EQ. 'ca') THEN
                     anxca=.TRUE.
                  ELSE IF(strngs(i) .EQ. 'backbone') THEN
                     anxbc=.TRUE.
                  ELSE IF(strngs(i) .EQ. 'heavy') THEN
                     anxhe=.TRUE.
                  ELSE IF(strngs(i) .EQ. 'allatoms')
     &                    THEN
                     anxal=.TRUE.
                  ELSE IF(strngs(i) .EQ. 'all') THEN
                     anxca=.TRUE.
                     anxbc=.TRUE.
                     anxhe=.TRUE.
                     anxal=.TRUE.
                  ELSE
                     errmsg=err_unr(3) // strngs(i)
                     CALL xerror(errmsg,80,1,30)
                     nsevere = nsevere + 1
                  END IF
               END DO


c==== Subcommand select_cell=========================================

            ELSE IF(strngs(1).EQ. 'select_cell') THEN
               anxrms_cell=.TRUE.

c==== Subcommand select_atms=========================================

            ELSE IF(strngs(1).EQ. 'select_atms') THEN
               latms=.TRUE.
               natms=natms+1
               IF(nword .NE. 3) THEN
                  errmsg=err_args(1)//'2'
                  CALL xerror(errmsg,80,1,30)
                  nsevere = nsevere + 1
               ELSE
                  DO i=2,3
                     CALL fndfmt(1,strngs(i),fmt)
                     READ(strngs(i),fmt,err=20) patms(i-1,natms)
                  END DO
               END IF

c==== Subcommand averaged ===========================================

            ELSE IF(strngs(1).EQ. 'averaged') THEN
               IF(nword .LT. 2) THEN
                  errmsg=err_args(1)//'1'
                  CALL xerror(errmsg,80,1,30)
                  nsevere=nsevere+1
               END IF
               IF(strngs(2) .EQ. 'ca') THEN
                  avg_ca=.TRUE.
               ELSE IF(strngs(2) .EQ. 'heavy') THEN
                  avg_he=.TRUE.
               ELSE IF(strngs(2) .EQ. 'backbone') THEN
                  avg_bc=.TRUE.
               ELSE
                  errmsg=err_unr(3)//strngs(2)//err_end(1:14)
     &                 //err_end(16:20)
                  CALL xerror(errmsg,80,1,30)
                  nsevere = nsevere + 1
               END IF

            ELSE IF(strngs(1).EQ. 'eul_angles') THEN
               eul_angles=.TRUE.
               CALL EUL_Create(nores)
               strngs(2)='EUL_DATA'
               INQUIRE(FILE=strngs(2),EXIST=exist)
               IF(exist) THEN
                  CALL openf(EUL_kdump,strngs(2),'FORMATTED',
     &                 'OLD',0)
               ELSE
                  CALL openf(EUL_kdump,strngs(2),'FORMATTED'
     &                 ,'NEW',0)
               END IF

            ELSE IF(strngs(1) .EQ. 'eul_domain') THEN
               eul_angles=.TRUE.
               m=0
               CALL parse_numbers(err_unr,strngs,nword,EUL_index_l,m
     &              ,nsevere)
               EUL_index_l(1)=m
               CALL EUL_Add

c==== Subcommand print rms ==========================================

            ELSE IF(strngs(1).EQ. 'rms_matrix') THEN
               rms_matrix=.TRUE.
               avg_ca=.TRUE.
               strngs(2)='RMS_MATRIX'
               INQUIRE(FILE=strngs(2),EXIST=exist)
               IF(exist) THEN
                  CALL openf(krms_matrix,strngs(2),'FORMATTED',
     &                 'OLD',0)
               ELSE
                  CALL openf(krms_matrix,strngs(2),'FORMATTED'
     &                 ,'NEW',0)
               END IF

            ELSE IF(strngs(1).EQ. 'rms_matrix_plot') THEN
               rms_matrix=.TRUE.
               rms_matrix_plot=.TRUE.
               strngs(2)='RMS_MATRIX.pdb'
               INQUIRE(FILE=strngs(2),EXIST=exist)
               IF(exist) THEN
                  CALL openf(krms_matrix_plot,strngs(2),'FORMATTED',
     &                 'OLD',0)
               ELSE
                  CALL openf(krms_matrix_plot,strngs(2),'FORMATTED'
     &                 ,'NEW',0)
               END IF
               
            ELSE IF(strngs(1).EQ. 'rms_matrix_avg') THEN
               rms_matrix=.TRUE.
               rms_matrix_avg=.TRUE.

            ELSE IF(strngs(1).EQ. 'rms_matrix_freq') THEN
               rms_matrix=.TRUE.
               CALL fndfmt(2,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) dummy
               Write_Freq=INT(dummy)
               RMS_Navg=INT(dummy)

c==== Subcommand print rms ==========================================

            ELSE IF(strngs(1).EQ. 'print' .AND. strngs(2) .EQ. 'rms')
     &              THEN
               IF(nword .LT. 5) THEN
                  errmsg=err_args(1)//'4'
                  CALL xerror(errmsg,80,1,30)
                  nsevere=nsevere+1
               END IF
               IF(strngs(4) .EQ. 'OPEN') THEN
                  avg_rms=.TRUE.
                  CALL fndfmt(2,strngs(3),fmt)
                  READ(strngs(3),fmt) frms
                  CALL uscrpl(strngs(5),80)
                  INQUIRE(FILE=strngs(5),EXIST=exist)
                  IF(exist) THEN
                     CALL openf(krms,strngs(5),'FORMATTED',
     &                    'OLD',0)
                  ELSE
                     CALL openf(krms,strngs(5),'FORMATTED'
     &                    ,'NEW',0)
                  END IF
               ELSE
                  errmsg=err_open
                  CALL xerror(errmsg,80,1,30)
                  nsevere = nsevere + 1
               END IF



            ELSE IF(strngs(1) .EQ. ' ') THEN
               GOTO 1100
               
            ELSE IF(strngs(1).EQ. 'END' ) THEN
               GOTO 100

            ELSE
               errmsg=err_unr(2)//strngs(2)//err_end(1:14)/
     &              /err_end(16:20)
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF
            GOTO 1100

c==== Command TOPOLOGY================================================

      ELSE IF(strngs(1).EQ. 'TOPOLOGY' ) THEN
         not_time_corr=.TRUE.
         prttopl=.TRUE.
900      READ(knlist,'(a78)',END=600) line(1:78)
         CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 900
         CALL parse(line,sep,2,comm,strngs,40,nword,
     x        iret,errmsg)
         IF(strngs(1) .EQ. 'print' ) THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) ftop_print
            IF(strngs(3) .EQ. 'OPEN') THEN
               CALL uscrpl(strngs(4),80)
               INQUIRE(FILE=strngs(4),EXIST=exist)
               IF(exist) THEN
                  CALL openf(ktopol,strngs(4),'FORMATTED','OLD',0)
               ELSE
                  CALL openf(ktopol,strngs(4),'FORMATTED','NEW',0)
               END IF
            ELSE
               errmsg='OPEN keyword not found'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. 'bonds') THEN
            CALL parse_numbers(err_unr,strngs,nword,top_bonds,n_bonds
     &           ,nsevere)
            IF(n_bonds .GT. ntopol) THEN
               errmsg=' Length of the topology-analysis arrays '
     &              / /'insufficient. Increase _ANL_TOPO_.'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. 'bendings') THEN
            CALL parse_numbers(err_unr,strngs,nword,top_bendings
     &           ,n_bendings,nsevere)
            IF(n_bendings .GT. ntopol) THEN
               errmsg=' Length of the topology-analysis arrays '
     &              / /'insufficient. Increase _ANL_TOPO_.'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. 'P-torsions') THEN
            CALL parse_numbers(err_unr,strngs,nword,top_ptors,n_ptors
     &           ,nsevere)
            IF(n_ptors .GT. ntopol) THEN
               errmsg=' Length of the topology-analysis arrays '
     &              / /'insufficient. Increase _ANL_TOPO_.'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. 'I-torsions') THEN
            CALL parse_numbers(err_unr,strngs,nword,top_itors,n_itors
     &           ,nsevere)
            IF(n_itors .GT. ntopol) THEN
               errmsg=' Length of the topology-analysis arrays '
     &              / /'insufficient. Increase _ANL_TOPO_.'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. ' ') THEN
            GOTO 900

         ELSE IF(strngs(1).EQ. 'END' ) THEN
            GOTO 100
         ELSE
            errmsg=err_unr(3)//strngs(1)//err_end(1:14)//err_end(16:20)
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         GOTO 900
         
c==== Command  EQUI =================================================

      ELSE IF(strngs(1).EQ. 'EQUI' ) THEN
          equi=.TRUE.
          IF(strngs(2).EQ.'OPEN') THEN
             IF(nword.ne.3) THEN  
                nsevere = nsevere+1
                errmsg=err_args(1)//'3'
                call xerror(errmsg,80,1,30)
             ELSE
                CALL uscrpl(strngs(3),80)
                
                INQUIRE(FILE=strngs(3),EXIST=exist)
                IF(exist) THEN
                   CALL openf(kequi,strngs(3),'FORMATTED','OLD',0)
                   REWIND kequi
                ELSE
                   CALL openf(kequi,strngs(3),'FORMATTED','NEW',0)
                END IF

             END IF
          END IF
         
c==== Command FORCE_FIELD ============================================

      ELSE IF(strngs(1) .EQ. 'FORCE_FIELD') THEN
         write_ff_pars=.TRUE.
         
c==== Command GYRATION ===============================================

      ELSE IF(strngs(1).EQ. 'GYRATION') THEN
         not_time_corr=.TRUE.
         wrtgyr=.TRUE.
         CALL fndfmt(2,strngs(2),fmt)
         READ(strngs(2),fmt,err=20) fgyr
         IF(strngs(3).EQ. 'OPEN') THEN
            
            CALL uscrpl(strngs(4),80)
            INQUIRE(FILE=strngs(4),EXIST=exist)
            IF(exist) THEN
               CALL openf(kgyr,strngs(4),
     &              'FORMATTED','OLD',0)
            ELSE
               CALL openf(kgyr,strngs(4),
     &              'FORMATTED','NEW',0)
            END IF
         ELSE
            errmsg=err_open
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         
      ELSE IF(strngs(1).EQ. ' ') THEN
         CONTINUE
         
c==== Begininning of next ENVIRONMENT =================================

      ELSE IF(strngs(1)(1:1).EQ. '&'.AND.strngs(1).NE. '&END') THEN
         errmsg=err_unr(1) //strngs(1)(1:8)// err_end
         CALL xerror(errmsg,80,1,30)
         nsevere = nsevere + 1
         GO TO 600

         CONTINUE
         
c==== Command &END ====================================================

      ELSE IF(strngs(1).EQ. '&END') THEN
         GOTO 600
         
      ELSE
         errmsg=err_unr(1) //strngs(1)(1:8)// err_end
         CALL xerror(errmsg,80,1,30)
         nsevere = nsevere + 1
      END IF

      GO TO 100

600   CONTINUE

c=======================================================================
c     Environment parser ends here 
c=======================================================================

c--   syntax errors: abort without verifying input 
      if(nsevere.gt.0.and.nsevere.lt.99) then 
         j=0
         call int_str(nsevere,fmt,j)
         errmsg=fmt(1:j) //' ERRORS WHILE EXECUTING READ_PROPERTIES'
         CALL xerror(errmsg,80,1,2)
         STOP
      ELSE IF(nsevere.gt.99) THEN 
         errmsg= 'MORE THAN 99 ERRORS WHILE EXECUTING READ_PROPERTIES'
         call xerror(errmsg,80,1,2)
         STOP
      END IF
      if(nwarning.gt.0.and.nwarning.lt.99) then 
         iret=0
         j=0
         call int_str(nwarning,fmt,j)
         errmsg= fmt(1:j)//' WARNINGS WHILE EXECUTING READ_PROPERTIES'
         CALL xerror(errmsg,80,1,1)
      ELSE IF(nwarning.gt.99) THEN 
         errmsg= 'MORE THAN 99 WARNINGS WHILE EXECUTING READ_PROPERTIES'
         call xerror(errmsg,80,1,1)
      ENDIF    
      IF(iret.eq.1) THEN
         errmsg='while parsing line: toomany strings'
         CALL xerror(errmsg,80,1,30)
         j=0
         call int_str(nsevere,fmt,j)
         errmsg=fmt(1:j) //' ERRORS WHILE EXECUTING READ_PROPERTIES'
         CALL xerror(errmsg,80,1,2)
         STOP
      END IF
      voronoi_res(1)=n_res
      top_bonds(1)=n_bonds
      top_bendings(1)=n_bendings
      top_ptors(1)=n_ptors
      top_itors(1)=n_itors
      corr_atoms(1)=n_atoms
      IF(hydration) THEN
         st_dummy(1)=n_st_dummy
         sv_dummy(1)=n_sv_dummy
         CALL HYD_Initialize(st_dummy,sv_dummy,khydration)
         DEALLOCATE(st_dummy,sv_dummy)
      END IF
      
      RETURN

c==============================================================================
c     Errors were found
c==============================================================================


 20   CONTINUE
      iret=1
      errmsg='internal reading error: wrong format?? TAB character??'
      CALL xerror(errmsg,80,1,2)
      RETURN

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      END
