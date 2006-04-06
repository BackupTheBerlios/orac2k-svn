#include "pressure.h"
      SUBROUTINE verify_input(fmaxstp,fscale,fprint,fsave,fplot,fascii
     &     ,fprop,fplot_fragm,fplot_center,frject,fconf,fmaxrun,fupdte
     &     ,fxrms,fvi,favg,favg_xrms,fcavities,ftop_print
     &     ,gofr_fprint,gofr_favg,gofr_fcomp,sofk_fprint,sofk_fcomp
     &     ,fprtvaf,ftotvaf,fnovaf,fdipole,fnative,ffragm_dist,fhbonds
     &     ,fhhisto,frms,fgyr,fabmd,freq_ef,freq_dp,fxslt,finst_fit
     &     ,fcalc_cofm,finst_lda,flda_flu,flda_hyd,fprot_hyd,fprot_lda
     &     ,iret,errmsg)

************************************************************************
*   Time-stamp: <97/07/07 11:44:20 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Nov 19 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*
      
      USE VORONOI_Mod, ONLY: voronoi,nvoronoi
      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8  fmaxstp,fscale,fprint,fsave,fplot,fascii,fprop,fplot_fragm
     &     ,fplot_center,frject,fconf,fmaxrun,fupdte,fxrms,favg
     &     ,favg_xrms,fcavities,ftop_print,gofr_fprint
     &     ,gofr_favg,gofr_fcomp,fprtvaf,ftotvaf,fnovaf,fdipole,fnative
     &     ,fvi,ffragm_dist,fhbonds,fhhisto,frms,fgyr,fabmd,freq_ef
     &     ,freq_dp,fxslt,finst_fit,fcalc_cofm,finst_lda,flda_flu
     &     ,flda_hyd,fprot_hyd,fprot_lda,sofk_fprint,sofk_fcomp
      INTEGER iret,nwarning,nsevere,ierr
      CHARACTER*80 errmsg

*----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'parst.h'
      INCLUDE 'cpropar.h'
      INCLUDE 'unit.h'
      INCLUDE 'pme.h'

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8 n_timestep,aux,cut,tsign
      LOGICAL abort1,abort1a,abort1b,abort1c,abort2,near0,ok
      INTEGER i,idum1,idum2,n
      CHARACTER*11 integers 
      CHARACTER*2 nerr

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      aux=DABS(time)
      tsign=aux/time
      time=aux
      WRITE(kprint,'(5x,a)') 'Checking Input ......              ---->'
      integers=" 1234567890"
      nerr="  "
      nwarning = 0
      nsevere  = 0
      IF(.NOT. replicate) THEN
         aaxis=aaxis/DFLOAT(icl)
         baxis=baxis/DFLOAT(icm)
         caxis=caxis/DFLOAT(icn)
      END IF
         
      IF(slt_exist .AND. slv_create) THEN
         slt_create=.TRUE.
      END IF

      IF(slt_exist .AND. sgroup) THEN
         slt_create=.TRUE.
      END IF

      if ( slt_exist . AND . anprot ) THEN

         if ( annpro . gt . ndan1 ) then
             errmsg=
     &'In DEF_SOLUTE(&SOLUTE): SUBUN(config.h) must be encreased'
                     call xerror(errmsg,80,1,20)
                     nsevere  = nsevere  + 1 
         end if

         if ( annpro . gt . 1 ) then
            do n = 1, annpro-1
               if ( anpoint(2,n) . gt . m1 ) then
                     errmsg=
     &'In DEF_SOLUTE(&SOLUTE) one sup_limit > number of atoms'
                     call xerror(errmsg,80,1,20)
                     nsevere  = nsevere  + 1
               end if
            end do
         end if
         if ( anpoint(2,annpro) . gt . m1 ) then
               errmsg=
     &'In DEF_SOLUTE(&SOLUTE) one sup_limit > number of atoms'
               call xerror(errmsg,80,1,20)
               nsevere  = nsevere  + 1
         end if
      end if

      IF(efield) THEN
         rcutm = rspoff
         rtolm = 0.0d0
      END IF

      IF(debug.and.tpgfil) THEN 
         errmsg=
     &'DEBUG(&RUN) ignored when READ_PFR_BIN(&PARAMETERS) is'
     & / / ' specified'       
         call xerror(errmsg,80,1,21)
         nwarning   = nwarning  + 1 
      END IF  
      IF(grpcut .AND. pme) THEN
         errmsg=
     &'PME and GROUP_CUTOFF(&POTENTIAL) are mutually exclusive'
         call xerror(errmsg,80,1,20)
         nsevere  = nsevere  + 1 
      END IF

      IF(grpcut .AND. clewld) THEN
         errmsg
     &='EWALD and GROUP_CUTOFF(&POTENTIAL) are mutually exclusive'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(solven.and.slvpdb.and.linser) THEN
         errmsg=
     &'READ_PDB(&SOLVEN) AND INSERT(&SETUP) have been both specified'
         call xerror(errmsg,80,1,21)
         nwarning = nwarning + 1
      END IF
      
      IF(pme .AND. md_respa) THEN
         IF(shell_pme .EQ. '0') THEN
            errmsg=
     &'From MTS_RESPA(&INTEGRATOR) no shell defined for PME'
            call xerror(errmsg,80,1,21)
            nwarning = nwarning + 1
         END IF
         IF(nfft1 .EQ. 0) THEN
            errmsg=
     & 'In EWALD(&POTENTIAL) no grid speicified for pme.'  
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
      END IF

      IF(clewld .AND. md_respa .AND. (.NOT. pme) ) THEN
         IF(shell_pme .EQ. '0') THEN
            errmsg=
     & 'From MTS_RESPA(&INTEGRATOR) no shell defined'
            call xerror(errmsg,80,1,21)
            nwarning = nwarning + 1
         END IF
      END IF
      
      IF(pme) THEN
         IF(nfft1+1 .GT. mfft1 .OR. nfft2+1 .GT. mfft2 .OR. nfft3+1 .GT.
     &        mfft3) THEN
            errmsg=
     &'pme grid nfft1-2-3 in EWALD(&POTENTIAL) exceed physical'
     &/ /' dimensions'
            call xerror(errmsg,80,1,21)
            nwarning=nwarning+1
c$$$            call xerror(errmsg,80,1,20)
c$$$            nsevere = nsevere + 1 
         END IF
      END IF
      IF(resize_cell .AND. change_cell) THEN
         errmsg=
     & 'In &SETUP. Cannot run with RESIZE_CELL and CHANGE_CELL'
     & / /' options in the same input.'
            call xerror(errmsg,80,1,20)
         nsevere=nsevere+1
      END IF

#ifndef PRESSURE
      IF(cpress .OR. pressure .OR. isostress) THEN
         errmsg=
     & 'In &SIMULATION(STRESS) should recompile with pressure option.'
     & / /' Change pressure.h'
            call xerror(errmsg,80,1,20)
         nsevere=nsevere+1
      END IF
#endif

*=======================================================================
*--- Validate steps                                                  ---
*=======================================================================

      IF(minimize) THEN
         time=1.0D0
         nrespa=1
         n0respa=1
         n1respa=1
         mrespa=1
         lrespa=1
      END IF

      IF( .NOT. analys) THEN
         maxstp=IDINT(fmaxstp/time)
         nprop=IDINT(fprop/time)
         nrject=IDINT(frject/time)
         nsave=IDINT(fsave/time)
         nupdte=IDNINT(fupdte/time)
         nheating=NINT(fscale/time)
         n_timestep=time/DFLOAT(mrespa*lrespa)
         maxrun=IDINT(fmaxrun/n_timestep)
         nprint=IDNINT(fprint/n_timestep)
         nconf=IDNINT(fconf/n_timestep)
         nplot_fragm=NINT(fplot_fragm/n_timestep)
         nplot_center=NINT(fplot_center/n_timestep)
         nascii=IDNINT(fascii/n_timestep)
         nplot=IDNINT(fplot/n_timestep)
         nxrms=IDNINT(fxrms/n_timestep)
         nxslt=IDNINT(fxslt/n_timestep)
         nvi=IDNINT(fvi/n_timestep)
         ndipole=IDNINT(fdipole/n_timestep)
         nnative=IDNINT(fnative/n_timestep)
         nfragm_dist=IDNINT(ffragm_dist/n_timestep)
         nhbonds=IDNINT(fhbonds/n_timestep)
         nhhisto=IDNINT(fhhisto/n_timestep)
         nrms=IDNINT(frms/n_timestep)
         ngyr=IDNINT(fgyr/n_timestep)
         navg_str=IDNINT(favg/n_timestep)
         navg_str_xrms=IDNINT(favg_xrms/n_timestep)
         ncavities=IDNINT(fcavities/n_timestep)
         ntop_print=IDNINT(ftop_print/n_timestep)
         gofr_nprint=NINT(gofr_fprint/n_timestep)
         gofr_navg=NINT(gofr_favg/n_timestep)
         gofr_ncomp=IDINT(gofr_fcomp/n_timestep)
         nabmd=IDINT(fabmd/n_timestep)
         sofk_nprint=NINT(sofk_fprint/n_timestep)
         sofk_ncomp=IDINT(sofk_fcomp/n_timestep)
      ELSE
         IF(dmprnd_o) THEN
            maxstp=IDINT(fmaxstp/fconf)
            maxrun=IDINT(fmaxrun/fconf)
            nconf=1
         END IF

         nplot_fragm=IDNINT(fplot_fragm)
         nplot_center=IDNINT(fplot_center)
         nascii=IDNINT(fascii)
         nplot=IDNINT(fplot)
         nxrms=IDNINT(fxrms)
         nxslt=IDNINT(fxslt)
         nplot_fragm=IDNINT(fplot_fragm)
         nfragm_dist=IDNINT(ffragm_dist)
         ndipole=IDNINT(fdipole)
         nnative=IDNINT(fnative)
         navg_str=IDNINT(favg)
         navg_str_xrms=IDNINT(favg_xrms)
         ncavities=IDNINT(fcavities)
         ntop_print=IDNINT(ftop_print)
         gofr_nprint=IDNINT(gofr_fprint)
         gofr_navg=IDNINT(gofr_favg)
         gofr_ncomp=IDNINT(gofr_fcomp)
         nhbonds=IDNINT(fhbonds)
         nhhisto=IDNINT(fhhisto)
         nrms=IDNINT(frms)
         ngyr=IDNINT(fgyr)
         nfreq_ef=IDNINT(freq_ef)
         nfreq_dp=IDNINT(freq_dp)
         ninst_fit=IDNINT(finst_fit)
         ninst_lda=IDNINT(finst_lda)
         ncalc_cofm=IDNINT(fcalc_cofm)
         nlda_flu=IDNINT(flda_flu)
         nlda_hyd=IDNINT(flda_hyd)
         nprot_hyd=IDNINT(fprot_hyd)
         nprot_lda=IDNINT(fprot_lda)
         sofk_nprint=IDNINT(sofk_fprint)
         sofk_ncomp=IDNINT(sofk_fcomp)
         IF(vacf) THEN
            prtvaf=IDINT(fprtvaf/n_timestep)
            prtvaf=prtvaf/novaf
         END IF
      END IF
      IF(.NOT. stoprun) THEN
         IF( nprint .EQ. 0) THEN
            errmsg=
     &'PRINT(&RUN) is 0 or less than TIMESTEP(&VARIABLES);'/ /
     &           ' NPRINT set to 1.'
            nprint = 1
            call xerror(errmsg,80,1,21)
            nwarning=nwarning+1
         END IF
         IF( nprop .EQ. 0) THEN
            errmsg=
     &'PROPERTY(&RUN) is O or less than TIMESTEP(&VARIABLES)'/ /
     &'; NPROP set to 1000'
            call xerror(errmsg,80,1,21)
            nwarning=nwarning+1
            nprop=1000
         END IF
         IF( nupdte .EQ. 0 .AND. (.NOT. analys)) THEN
            errmsg=
     &'UPDATE(&POTENTIAL) is 0 or less than TIMESTEP'/ /
     &'(&VARIABLES); NUPDTE set to 4'
            call xerror(errmsg,80,1,21)
            nupdte=4
            nwarning = nwarning + 1 
         END IF
            
         IF(gofr) THEN
            IF(gofr_ncomp .EQ. 0) THEN
               errmsg=
     &'compute(GOFR&PROPERTIES) is 0 or less'/ / 
     &'than TIMESTEP(&VARIABLES); reset to 1'
               call xerror(errmsg,80,1,21)
               nwarning = nwarning + 1 
               gofr_ncomp = 1
            END IF
            IF(gofr_nprint .EQ. 0) THEN
               errmsg=
     &'print(GOFR&PROPERTIES) is 0 or less than '/ /
     &'TIMESTEP(&VARIABLES); reset to 1'
               call xerror(errmsg,80,1,21)
               nwarning = nwarning + 1
               gofr_nprint = 1
            END IF
         END IF
      END IF
      IF(linked_cell) THEN
         IF(nupdte_index .EQ. 0) THEN
            errmsg=' Frequency of calls to index cell list is zero. '
     &           / /'Check LINKED_CELL(&POTENTIAL). '
            CALL xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
      END IF
      IF((.NOT. bending) .AND. stretch) THEN
          errmsg='Cannot use stretching and bending at the'
     x  / /' same time. Abort.'
          CALL xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(rspoff .EQ. 0.0D0 .AND. (.NOT. md_respa) .AND. (.NOT. analys)
     &     .AND. (.NOT. stoprun))THEN
         errmsg=
     &'Should specify a cutoff in &POTENTIAL (commands CUTOFF'
     &/ /' or GROUP_CUTOFF)'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(gofr .AND. gofr_neighbor) THEN
         IF(.NOT. md_respa) THEN
            cut=rspoff
            IF(rspcut+cut .LT. gofr_cut) THEN
               gofr_neighbor=.FALSE.
            END IF
            IF(gofr_cut .EQ. 0.0D0) gofr_cut=cut
         ELSE
            cut=rcuth+rtolh
            IF(rneih+cut .LT. gofr_cut) THEN
               gofr_neighbor=.FALSE.
            END IF
            IF(gofr_cut .EQ. 0.0D0) gofr_cut=rcuth+rtolh
         END IF
      END IF

      IF(gofr) THEN
         IF(NINT(gofr_cut/delrg) .GT. maxint) THEN
            errmsg=
     &'**ERROR: Dimensions of pair correlation are insufficient.'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
      END IF
         
*=======================================================================
*----- Consistency checks ---------------------------------------------
*=======================================================================

      IF(minimize) THEN
#ifdef PARALLEL
         errmsg=
     &'The minimize routines have not been parallelized'/ /
     &' yet.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
#endif
         IF(near0(rspoff)) THEN
            errmsg=
     &'Use CUTOFF(&POTENTIAL) to assign a direct space'/ /
     &' cutoff.'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
         rcuth=rspoff-rspon
         rtolh=rspon
         rtoll=rcuth
         rcutl=rspoff-rspon
         rtoll=0.0
         rneih=rspcut+cut
         rneil=0.0
         IF(md_respa) THEN
            errmsg=
     &'Warning MTS_STEP(&INTEGRATOR) option is'/ /
     &' ignored by MINIMIZE'
            call xerror(errmsg,80,1,21)
            nwarning = nwarning + 1
         END IF
         IF(md_respa) THEN
            errmsg=
     &'Warning TIMESTEP(&INTEGRATOR) option is'/ /
     &' ignored by MINIMIZE.Timestep set to 1 fs'
            call xerror(errmsg,80,1,21)
            nwarning = nwarning + 1
         END IF
         md_respa=.FALSE.
         IF(thermos) THEN
            errmsg=
     &'Warning THERMOS(&SIMULATION) option is'/ /
     &' ignored by MINIMIZE'
            call xerror(errmsg,80,1,21)
            nwarning = nwarning + 1
         END IF
         IF(cpress .AND. (.NOT. l_bfgs_b)) THEN
            errmsg=
     &'The (ISO)STRESS(&SIMULATION) option '/ /
     &'can be run only by MINIMIZE(BFGS).'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1
         END IF
         IF(maxstp .EQ. 0) THEN
            errmsg=
     &'The number of minimization steps is zero. '/ /
     &' Use the TIME(&RUN) option to change steps.'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
         IF(.NOT. stretch .OR. stretch_heavy) THEN
            errmsg=
     &'Minimization cannot be run with constraints.'/ /
     &' Change STRETCH options.'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
      END IF
      IF(frequencies .AND. mdsim) THEN
         errmsg=
     &'Calculation of harmonic frequencies is '/ /
     &'incompatible with MDSIM(&SIMULATION).'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(frequencies .AND. md_respa) THEN
         errmsg=
     &'Calculation of harmonic frequencies is '/ /
     &'incompatible with &INTEGRATOR.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(anxrms ) THEN
         IF( nxrms .EQ. 0 .AND. nxslt. EQ. 0 ) THEN
            errmsg='inst_xrms(COMPARE(&PROPERTIES)): Cannot compute'
     & / /       ' without an output file.'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
      END IF
      IF(prot_lda .AND. rspoff.le.0.0d0) THEN
         errmsg=
     &        'LDA(&PROPERTIES): Cannot calculate ' / /
     &        'solvation without '/ /
     &        '*CUTOFF* defined.'
         call xerror(errmsg,80,1,20)
           nsevere = nsevere + 1
        END IF

      IF(prot_lda .AND. update_anl.le.0) THEN
         errmsg=
     &        'LDA(&PROPERTIES): Cannot calculate ' / /
     &        'solvation without '/ /
     &        'UPDATE(&ANALYS) defined.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1
      END IF

      IF(prot_hyd .AND. prot_lda) THEN
         errmsg='&PROPERTIES:Cannot calculate ' / /
     &        'solvation(LDA) and ' / /
     &        'hydration(HYDRATION) simultaneously'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1
      END IF

      IF(prot_hyd .AND. rspoff.le.0.0d0) THEN
         errmsg=
     &        'HYDRATION(&PROPERTIES): Cannot calculate ' / /
     &        'hydration without '/ /
     &        '*CUTOFF* defined.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1
      END IF
      IF(prot_hyd .AND. update_anl.le.0) THEN
         errmsg=
     &        'HYDRATION(&PROPERTIES): Cannot calculate ' / /
     &        'hydration without '/ /
     &        'UPDATE(&ANALYS) defined.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1
      END IF

      IF(efield .AND. rspoff.le.0.0d0) THEN
         errmsg=
     &        'ELECTRIC(&PROPERTIES): Cannot calculate ' / /
     &        'electric field without '/ /
     &        '*CUTOFF* defined.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1
      END IF
      IF(efield .AND. update_anl.le.0) THEN
         errmsg=
     &        'ELECTRIC(&PROPERTIES): Cannot calculate ' / /
     &        'electric field without '/ /
     &        'UPDATE(&ANALYS) defined.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1
      END IF
      IF(efield .AND. .NOT. lmol_ef .AND. .NOT. polar) THEN
         errmsg=
     &        'ELECTRIC(&PROPERTIES): Cannot calculate' / /
     &        'electric field without *sel_mol* defined.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1
      END IF

      IF(polar .AND. .NOT. efield  .AND. analys) THEN
         errmsg=
     &        'POLARIZATION(&PROPERTIES): Cannot calculate' / /
     &        'polarization without ELECTRIC defined.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1
      END IF

      IF(polar .AND. fudge .eq. 1.d0 ) THEN
         write(kprint,*) 
         write(kprint,*) 'WARNING'
         write(kprint,*) 'QQ-FUDGE = 1.0d0'
         write(kprint,*) 'Could set another value in &POLARIZATION'
         write(kprint,*) 'WARNING'
         write(kprint,*) 
      END IF
      
      IF(polar .AND. polar_scale .eq. 1.d0 ) THEN
         write(kprint,*) 
         write(kprint,*) 'WARNING'
         write(kprint,*) 'Are you using enhanced charges?'
         write(kprint,*) 'SCALE (charges) = 1.0d0'
         write(kprint,*) 'Could set another value in &POLARIZATION'
         write(kprint,*) 'WARNING'
         write(kprint,*) 
      END IF
      

      IF((avg_lda .OR. lda_flu .OR. lda_hyd)
     &    .AND. .NOT. sel_lda) THEN
         errmsg=
     &        'STRUCTURES(&PROPERTIES): Cannot calculate' / /
     &        'lda properties without *sel_lda* defined.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1
      END IF

      IF(anxrms .AND. (.NOT. template)) THEN
         errmsg=
     &        'STRUCTURES(&PROPERTIES): Cannot compute X_rms without'
     &   / /     ' template file.'/ /' Provide TEMPLATE(&SETUP) cmd.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1
      END IF
      IF(avg_str .AND. (.NOT. template)) THEN
         errmsg=
     &        'STRUCTURES(&PROPERTIES): Running without' / /
     &        ' template file.'/ /' Provide TEMPLATE(&SETUP) cmd.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1
      END IF
      IF(avg_ca .OR. avg_he .OR. avg_bc) THEN
         IF( .NOT. template ) THEN
         errmsg=
     &        'STRUCTURES(&PROPERTIES): Running without' / /
     &        ' template file.'/ /' Provide TEMPLATE(&SETUP) cmd.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1
         END IF
      END IF
      IF(avg_ca .AND. avg_he) THEN
         errmsg=
     &        'STRUCTURES(&PROPERTIES): Cannot fit' / /
     &        ' structures to heavy and ca atoms at the same time.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1
      END IF
      IF(avg_ca .AND. avg_bc) THEN
         errmsg=
     &        'STRUCTURES(&PROPERTIES): Cannot fit' / /
     &        ' structures to backbone and ca atoms at the same time.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1
      END IF
      IF(avg_he .AND. avg_bc) THEN
         errmsg=
     &        'STRUCTURES(&PROPERTIES): Cannot fit' / /
     &        ' structures to backbone and heavy' / /
     &        ' atoms at the same time.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1
      END IF

      IF(avg_str .OR. avg_rms .OR. inst_fit .OR. calc_cofm .OR.
     &   avg_lda ) THEN
         IF(.NOT. avg_ca .AND. .NOT. avg_he .AND. .NOT. avg_bc) THEN
            errmsg=
     &           'STRUCTURES(&PROPERTIES): Unknown atoms to which' / /
     &           ' fit structure. Use *averaged*.'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1
         END IF
      END IF

      IF(linser .AND. nflag(1) .NE. 0) THEN
         errmsg=
     &'CONTROL(&RUN) must be 0 when INSERT(&SET_UP) is'/ /
     &' specified' 
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1
      END IF

*=======================================================================
*---  Abort if control parameters are wrong (abort also atomic system)--
*=======================================================================
      
      IF(md_respa) THEN

         abort2=sction .OR. noneq
         IF(abort2) THEN
            errmsg=
     &'NON_EQUILIBRIUM(&SOLUTE) and MTS_RESPA(&INTEGRATOR)'/ /
     &'are mutually exclusive'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
         aux=rcutl-rcuth
         abort1=aux.GT.1d-15
         IF(abort1) THEN
            errmsg=
     &'In MTS_RESPA(&INTEGRATOR) radius_l > radius_h.'/ /
     &' Choose radius_l less or equal to radius_h'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
         aux=rcutm-rcutl
         abort1=aux.GT.1d-15
         IF(abort1) THEN
            errmsg=
     &'In MTS_RESPA(&INTEGRATOR) radius_m > radius_l.'/ /
     &'Choose radius_m less or equal to radius_l'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
         aux=km-kl
         abort1=aux.GT.1d-15
         IF(abort1) THEN
            errmsg=
     &'In MTS_RESPA(&INTEGRATOR) wrong k-ewald directive.'/ /
     &' Choose km less or equal to kl' 
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
         aux=kl-rkcut
         abort1=aux.GT.1d-15
         IF(abort1) THEN
            WRITE(kprint,10230)
10230       FORMAT(/ /  
     &'**ERROR: in MTS_RESPA(&INTEGRATOR) wrong k-ewald directive.' /
     &'- choose kl less or equal to rkcut in EWALD(&POTENTIAL)')
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
      END IF
      IF(.NOT. md_respa) THEN
         IF(landersen) THEN
            errmsg=
     &'Directive ANDERSEN(&SIMULATION) can be used only'/ /
     &' when MTS_RESPA(&INTEGRATOR) is specified'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
      END IF
      IF(nflag(1) .NE. 0) THEN
         adjust_cnstr=.FALSE.
      END IF
      IF(readjust_cnstr) THEN
         adjust_cnstr=.TRUE.
      END IF

      IF(nflag(1) .NE. 0 .AND. slv_create) THEN
         errmsg=
     &   'Cannot create solvent and read a restart files'
     & / /' at the same time.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(sgroup .AND. nflag(1) .NE. 0) THEN
         errmsg=
     &'Cannot apply symmetry and read a restart file'
     &/ /' at the same time.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF
      IF(sgroup .AND. tpgfil) THEN
         errmsg=
     &'Cannot read a binary TPGPRM file and apply symmetry'
     &/ /' at the same time.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(nplot_fragm.gt.0.and.nfragm.eq.0) THEN 
         errmsg=
     &   'No fragment defined in &READ_SOLUTE:' / /
     &   ' no fragment will be printed'
         call xerror(errmsg,80,1,21)
         nwarning= nwarning+1
      END IF   

      IF(slv_create .AND. (.NOT. tpgprm_read) ) THEN
         errmsg=
     &   'Cannot create solvent without '
     & / /'READ_PRM_ASCII and READ_TPG_ASCII in &PARAMETERS.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(slv_create .AND. nbun_slv .EQ. 0) THEN
         errmsg=
     &   'Cannot create solvent without '
     & / /'defining residues with JOIN SOLVENT in &PARAMETERS.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(slv_add .AND. nbun_slv .EQ. 0) THEN
         errmsg=
     &   'Cannot add solvent without '
     & / /'defining residues with JOIN SOLVENT in &PARAMETERS.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(pdb_read .AND. slv_create) THEN
         errmsg=
     &   'Cannot create solvent and read '
     & / /'from the PDB file defined in READ_PDB.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(pdb_read .AND. slt_create) THEN
         errmsg=
     &   'Cannot generate solvent and read '
     & / /'solute from the PDB file defined in READ_PDB.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(.NOT. analys) THEN
         IF((.NOT. pdb_read .AND. .NOT. slt_create .AND. .NOT.
     &        slv_create).AND. nflag(1) .EQ. 0 .AND. .NOT. stoprun) THEN
            errmsg=
     &   'Cannot run without initial coordinates '
     & / /'from the PDB file defined in READ_PDB.'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
      END IF

      IF(kcoord_slt .EQ. 0 .AND. slt_create) THEN
         errmsg=
     &   'When creating solvent, solute coords. are read'
     & / /' by command COORDINATES in &SOLUTE.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(stoprun .AND. slv_create) THEN
         errmsg=
     &   'Cannot use STOP_RUN (&PARAMETERS) and'
     & / /' create solvent at the same time.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(slv_add .AND. sgroup) THEN
         errmsg=
     &   ' Command SGROUP (&SOLUTE) and ADD_UNITS (&SOLVENT) are'
     & / /' incompatible.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(slv_add .AND. slv_create) THEN
         errmsg=
     &   ' Cannot create solvent and use ADD_UNITS (&SOLVENT) '
     & / /'at the same time.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF


      IF(abmd .AND. DABS(spring) .LT. 1.0D-3) THEN
         errmsg=
     &   'In &POTENTIAL: The ABMD force constant is zero.'
     & / /' Hope this is ok..'
         call xerror(errmsg,80,1,21)
         nwarning= nwarning+1
      END IF

      IF(abmd .AND. diss_atoms) THEN
         errmsg=
     &   'In &POTENTIAL: In this implementation Cannot use'
     & / /'ABMD with subcommand (atoms).'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1
      END IF

      IF(abmd .AND. (dissociate .AND. associate)) THEN
         errmsg=
     &   'In ABMD(&POTENTIAL): Cannot use'
     & / /' dissociate and associate at the same time.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1
      END IF
      IF(((.NOT. restart_read) .AND. (.NOT. restart_write)) .AND.
     &     nflag(1) .NE. 0) THEN
         errmsg='Need restart file to run. Use '
     &   / /'&INOUT(RESTART)) do define file.'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1
      END IF
      IF((.NOT. restart_read) .AND. restart_write .AND. nflag(1) .NE. 0)
     &     THEN
         restart_in=restart_out
      END IF
      IF(restart_read .AND. (.NOT. restart_write) .AND. nflag(1) .NE. 0)
     &     THEN
         errmsg='You are reading from a restart file, but writing'
     &        / /' to none. Is this is correct?'
         call xerror(errmsg,80,1,21)
         nwarning=nwarning+1
      END IF
      IF(restart_read .AND. nflag(1) .EQ. 0)
     &     THEN
         errmsg='The simulation starting from scratch, '
     &        / /'&INOUT(RESTART(read)) has no effect.'
         call xerror(errmsg,80,1,21)
         nwarning=nwarning+1
      END IF
      IF(abmd) THEN
         ok=.TRUE.
         IF(fold .AND. dissociate) ok=.FALSE.
         IF(fold .AND. associate) ok=.FALSE.
         IF(fold .AND. abmd_tors) ok=.FALSE.
         IF(dissociate .AND. abmd_tors) ok=.FALSE.
         IF(associate .AND. abmd_tors) ok=.FALSE.
         IF(.NOT. ok) THEN
            errmsg=
     &           'In ABMD(&POTENTIAL): Cannot use'
     & / /' two ABMD options at the same time.'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1
         END IF
      END IF
      IF(coupl_mol .AND. coupl_grp) THEN
         errmsg=
     &        'In COUPLING(&SIMULATION): Can use only one type of '
     &     / /   'coupling.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1
      END IF
         
#ifdef PARALLEL
      IF(.NOT. analys .AND. coupl_mol) THEN
         errmsg=
     &        'In COUPLING(&SIMULATION): You are using group '
     &        / /'coupling in ORAC''s parallel version.'
         call xerror(errmsg,80,1,21)
         nwarning=nwarning+1
      END IF
#endif

*=======================================================================
*---  Error on properties ----------------------------------------------
*=======================================================================

      IF((.NOT. analys) .AND. near0(fconf) .AND. dmprnd_i) THEN
         errmsg=
     &   'Use TRAJECTORY(&INOUT) only with '
     &   / /' the environment &ANALYSIS.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF
      IF(analys .AND. (.NOT. dmprnd_i) .AND. (.NOT. pdb_read)) THEN
         errmsg=
     &   'This run does NOT analyse trajectories.'
     &   / /' Use TRAJECTORY(&INOUT) to analyse trajectories.'
         call xerror(errmsg,80,1,21)
         nwarning   = nwarning  + 1 
      END IF
      IF(check_native .AND. native) THEN
         errmsg=
     &        'Commands CHECK_NATIVE and NATIVE are incompatible.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF
         
      IF((.NOT. analys) .AND. maxrun .EQ. 0 .AND. dmprnd_o ) THEN
         errmsg=
     &   'Cannot dump trajectory without the '
     &   / /'length of the run. Use MAXRUN(&RUN).'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(analys .AND. gofr) THEN
         IF(gofr_neighbor .AND. update_anl .EQ. 0) THEN
            errmsg=
     &   'Cannot compute G of R with neighbor(GOFR(&ANALYSIS))'
     &   / /' without update frequency UPDATE(&ANALYSIS).'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
      END IF
      IF(analys .AND. voronoi  .AND. update_anl .EQ. 0) THEN
         errmsg=
     &   'Cannot compute Voronoi volumes without'
     &   / /' update frequency defined in UPDATE(&ANALYSIS).'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF
      
      IF(cavities .AND. .NOT. voronoi) THEN
         errmsg=
     &   'Cannot compute cavities without Voronoi calculation.'
     &   / /' Define Voronoi parameters.'
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF
      IF(cavities .AND. voronoi) THEN
         IF(nx_cav .EQ. 0 .OR. ny_cav .EQ. 0 .OR. nz_cav .EQ. 0) THEN
            errmsg=
     &   'The cavity grid is incomplete. On one or more directions'
     &   / /' the dimension is zero.'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
      END IF
      IF(cavities .AND. voronoi) THEN
         IF(MOD(ncavities,nvoronoi) .NE. 0) THEN
            IF(ncavities .LT. nvoronoi) THEN
               ncavities=nvoronoi
               errmsg=
     &'print[CAVITIES(&PROPERTIES)] is not a multiple of'/ /
     &'print[VORONOI(&PROPERTIES)]. Reset.'
               call xerror(errmsg,80,1,21)
               nwarning=nwarning+1
            ELSE
               errmsg=
     &'print[CAVITIES(&PROPERTIES)] is not a multiple of'/ /
     &'print[VORONOI(&PROPERTIES)]. Reset.'
               call xerror(errmsg,80,1,21)
               nwarning=nwarning+1
               ncavities=(ncavities/nvoronoi)*nvoronoi
            END IF
         END IF
         IF(cavities .AND. voronoi .AND. ncavities .EQ. 0) THEN
            ncavities=nvoronoi
            errmsg=
     &'print[CAVITIES(&PROPERTIES)] is zero. Reset.'
            call xerror(errmsg,80,1,21)
            nwarning=nwarning+1
         END IF
         IF(kcavities .EQ. 0) THEN
            errmsg=
     &'print[CAVITIES(&PROPERTIES)] no file given.'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
      END IF
      
      IF((hbonds_tot .OR. hbonds_res .OR. hbonds_vor) .AND. nhbonds .EQ.
     &     0) THEN
         errmsg=
     &' Output file for hydrogen bonds computation is undefined. '
         call xerror(errmsg,80,1,20)
         nsevere = nsevere + 1 
      END IF

      IF(nhhisto .NE. 0 .AND. .NOT. hhisto) THEN
         errmsg=
     &' Warning. &PROPERTIES(HBONDS[histogram]) not specified. '
     &/ /' Cannot print histogram.'
         call xerror(errmsg,80,1,21)
         nwarning   = nwarning  + 1 
      END IF

      IF(hbonds_vor) voronoi=.TRUE.

      IF(time_corr) THEN
         IF(not_time_corr) THEN
            errmsg=
     &'Can run TIME_CORRELATION(&PROPERTY) only without other '
     &           / /'properties.'
            call xerror(errmsg,80,1,20)
            nsevere = nsevere + 1 
         END IF
      END IF

*===== Some option in properties are incompatible with md_respa =======

      IF(md_respa) THEN
         IF(vacf) THEN
            errmsg=
     &           ' VACF(&PROPERTIES) has no action while running'
     &           / /' a simulation.'
            call xerror(errmsg,80,1,21)
            nwarning= nwarning+1
         END IF
         IF(diffusion) THEN
            errmsg=
     &           ' DIFFUSION(&PROPERTIES) has no action'
     &  / /' while running a simulation.'
            call xerror(errmsg,80,1,21)
            nwarning= nwarning+1
         END IF
         IF(fragm_dist) THEN
            errmsg=
     &           ' FRAGM_DIST(&PROPERTIES) has no action'
     &  / /' while running a simulation.'
            call xerror(errmsg,80,1,21)
            nwarning= nwarning+1
         END IF
         IF(hbonds_tot .OR. hbonds_res) THEN
            errmsg=
     &           ' HBONDS(&PROPERTIES) has no action'
     &  / /' while running a simulation.'
            call xerror(errmsg,80,1,21)
            nwarning= nwarning+1
         END IF
         IF(avg_rms) THEN
            errmsg=
     &           ' print rms in STRUCTURES(&PROPERTIES) has no action'
     &  / /' while running a simulation.'
            call xerror(errmsg,80,1,21)
            nwarning= nwarning+1
         END IF
            
      END IF
      time=tsign*time

*=======================================================================
*---  Error counts is over; build abort message if errors are found 
*=======================================================================

c--   if errors are severe then 

      if(nsevere.gt.0.and.nsevere.lt.99) then 
         iret=1 
         idum1 = nsevere/10
         idum2 = nsevere - 10*idum1 
         idum1=idum1+1
         idum2=idum2+1
         if(idum2.eq.1) idum2=11
         nerr = integers(idum1:idum1)/ /integers(idum2:idum2)
         errmsg= nerr/ /' ERRORS WHILE EXECUTING VERIFY_INPUT'
         call xerror(errmsg,80,1,2)
      ELSE IF(nsevere.gt.99) THEN 
         iret=1
         errmsg= 'MORE THAN 99 ERRORS WHILE EXECUTING VERIFY_INPUT'
         call xerror(errmsg,80,1,2)
      END IF

c--   if warnings are given then 

      if(nwarning.gt.0.and.nwarning.lt.99) then 
         iret=0 
         idum1 = nwarning/10
         idum2 = nwarning - 10*idum1 
         idum1=idum1+1
         idum2=idum2+1
         if(idum2.eq.1) idum2=11
         nerr = integers(idum1:idum1)/ /integers(idum2:idum2)
         errmsg= nerr/ /' ERRORS WHILE EXECUTING VERIFY_INPUT'
         call xerror(errmsg,80,1,1)
      ELSE IF(nwarning.gt.99) THEN 
         iret=1
         errmsg= 'MORE THAN 99 ERRORS WHILE EXECUTING VERIFY_INPUT'
         call xerror(errmsg,80,1,21)
#ifdef PARALLEL
         CALL MPI_FINALIZE(ierr)
#endif
         STOP
      END IF
      if(nsevere.eq.0.and.nwarning.eq.0) THEN
         WRITE(kprint,'(5x,a)') 'Input OK!!'
      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
