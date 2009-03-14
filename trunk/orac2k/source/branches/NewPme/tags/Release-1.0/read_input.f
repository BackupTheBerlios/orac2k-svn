      SUBROUTINE read_input(sitew,iret,errmsg,myid)

*======================================================================*
*                                                                      *
*                             Last update                              *
*                                                                      *
*     $Id$
*                                                                      *
*     Written by Massimo Marchi CE Saclay                              *
*                                                                      *
*                                                                      *
*     EXTERNALS : Free format input package, MAXA0.                    *
*                                                                      *
*                                                                      *
************************************************************************


*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER sitew,iret,myid
      CHARACTER*80 errmsg

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER nword
      LOGICAL lprint
      CHARACTER*22 err_open
      CHARACTER*37 err_args(3)
      CHARACTER*20 err_end 
      CHARACTER*27 err_unr(4)
      CHARACTER*15 err_fnf
      CHARACTER*80 line,strngs(40)
      CHARACTER*1 sep(2),comm(2)
      REAL*8  fmaxstp,fscale,fprint,fsave,fplot,fascii,fprop,fplot_fragm
     &     ,fplot_center,frject,fconf,fmaxrun,fupdte,fxrms,fvi,fvoronoi
     &     ,fcavities,ftop_print,gofr_fprint,gofr_favg,gofr_fcomp
     &     ,fprtvaf,ftotvaf,fnovaf,fdipole,fnative,favg,favg_xrms
     &     ,ffragm_dist,fhbonds,fhhisto,frms,fgyr,fabmd,freq_ef,freq_dp
     &     ,fxslt,finst_fit,fcalc_cofm,finst_lda,flda_flu,flda_hyd
     &     ,fprot_hyd,fprot_lda,sofk_fprint,sofk_fcomp
      DATA sep/' ',','/comm/'(',')'/

*-------------------- COMMON VARIABLES ---------------------------------

      INCLUDE 'parst.h'
      INCLUDE 'cpropar.h'

      INCLUDE 'unit.h'
      COMMON /card / jrec,nrec,jerr,nerr,nline,istrt
      INTEGER jrec,nrec,jerr,nerr,nline,istrt(80)

*==================== EXECUTABLE STATEMENTS ============================

*--------------  synatx errors messages---------------------------------
*   control error messages in all envs except PARAMETERS and ANALYSIS
  
      err_open = 'OPEN keyword not found'
      err_args(1) = 'Number of arguments must be at least '
      err_args(2) = 'Number of arguments must not exceed  '
      err_args(3) = 'Number of arguments must be only     '
      err_unr(1)='Unrecognized command  ---> '
      err_unr(2)='Unrecognized subcommand -> '
      err_unr(3)='Unrecognized keyword ----> '
      err_unr(4)='UNSUPPORTED  COMMAND ----> '
      err_end = '...or missing &END'  
      err_fnf = ' file not found' 

*------------- Set default values of some input parameters -------------
      CALL blkdta


#ifdef PARALLEL

*=======================================================================
*---- Input for parallel processing with MPI inteface ------------------
*=======================================================================

      IF(myid .NE. 0) THEN
         kprint=6
         open(unit=6,file="/dev/null")
      END IF
#endif

      CALL Open_Input(knlist)

      fmaxstp=0.0D0
      fscale=0.0D0
      fprint=0.0D0
      fsave=0.0D0
      fplot=0.0D0
      fascii=0.0D0
      fprop=0.0D0
      fplot_fragm=0.0D0
      fplot_center=0.0D0
      frject=0.0D0
      fconf=0.0D0
      fmaxrun=0.0D0
      fupdte=0.0D0
      fxrms=0.0D0
      fvi=0.0D0
      fvoronoi=0.0D0
      fcavities=0.0D0
      ftop_print=0.0D0
      gofr_fprint=0.0D0
      gofr_favg=0.0D0
      gofr_fcomp=0.0D0
      fprtvaf=0.0D0
      ftotvaf=0.0D0
      fnovaf=0.0D0
      fdipole=0.0D0
      fnative=0.0D0
      favg=0.0D0
      favg_xrms=0.0D0
      ffragm_dist=0.0D0
      fhbonds=0.0D0
      fhhisto=0.0D0
      frms=0.0D0
      fgyr=0.0D0
      freq_ef=0.0d0
      freq_dp=0.0d0

      lprint=.FALSE.
      WRITE(kprint,1)
      WRITE(kprint,2)
      WRITE(kprint,3)
      WRITE(kprint,2)
      WRITE(kprint,4)
      WRITE(kprint,2)
      WRITE(kprint,5)
      WRITE(kprint,2)
      WRITE(kprint,6)
      WRITE(kprint,2)
      WRITE(kprint,7)
      WRITE(kprint,2)
      WRITE(kprint,2)
      WRITE(kprint,2)
      WRITE(kprint,2)
      WRITE(kprint,10)
      WRITE(kprint,2)
      WRITE(kprint,11)
      WRITE(kprint,2)
      WRITE(kprint,2)
      WRITE(kprint,2)
      WRITE(kprint,12)
      WRITE(kprint,2)
      WRITE(kprint,13)
      WRITE(kprint,2)
      WRITE(kprint,2)
      WRITE(kprint,2)
      WRITE(kprint,2)
      WRITE(kprint,2)
      WRITE(kprint,2)
      WRITE(kprint,1)
      WRITE(kprint,16)

      line(79:80)='  '
      titld0(79:80)='  '

*----- Title card,$ 'blank'

c=======================================================================
c     Input parser starts here for file sys.mddata 
c=======================================================================

100   READ(knlist,'(a78)',END=600) line(1:78)
      CALL wrenc(kprint,line)
      IF(line(1:1) .EQ. '#') GOTO 100 
      CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
      IF(iret.EQ.1) THEN 
         errmsg='while parsing line: toomany strings'
         CALL xerror(errmsg,80,1,30)
         go to 100
      END IF

c=====Environment RUN ==================================================
      IF(strngs(1).EQ. '&RUN' ) THEN
         CALL read_run(fprop,frject,fmaxrun,fmaxstp,fprint,err_args
     &        ,err_unr,err_end)
c=======================================================================

c=====Environment INOUT ================================================
      ELSE IF(strngs(1).EQ. '&INOUT' ) THEN
         CALL read_inout(fconf,fsave,fplot,fplot_fragm,fplot_center
     &        ,fascii,err_open,err_args,err_end,err_unr)
c=======================================================================

c=====Environment SIMULATION============================================
      ELSE IF(strngs(1).EQ. '&SIMULATION' ) THEN
         CALL read_simulation(fscale,err_args,err_unr,err_open,err_end)
c=======================================================================

c=====Environment POTENTIAL=============================================
      ELSE IF(strngs(1).EQ. '&POTENTIAL' ) THEN
         CALL read_potential(fupdte,fabmd,err_args,err_unr,err_end)
c=======================================================================

c=====Environment PARALLEL==============================================
      ELSE IF(strngs(1).EQ. '&PARALLEL' ) THEN
         CALL read_parallel(err_args,err_unr,err_end)
c=======================================================================

c=====Environment SETUP=================================================
      ELSE IF(strngs(1).EQ. '&SETUP' ) THEN
         CALL read_setup(err_open,err_args,err_end,err_unr
     &        ,err_fnf)
c=======================================================================

c=====Environment INTEGRATOR============================================
      ELSE IF(strngs(1).EQ. '&INTEGRATOR' ) THEN
         CALL read_integrator(err_args,err_unr,err_end)
c=======================================================================

c=====Environment PARAMETERS============================================
      ELSE IF(strngs(1).EQ. '&PARAMETERS' ) THEN
         CALL read_parameters(iret,errmsg)
         IF(iret .NE. 0) THEN 
            CALL xerror(errmsg,80,1,2)
         END IF
c=======================================================================

c=====Environment SOLUTE================================================
      ELSE IF(strngs(1).EQ. '&SOLUTE' ) THEN
         CALL read_solute(err_args,err_unr,err_end,err_fnf,err_open)
c=======================================================================

c=====Environment SOLVENT===============================================
      ELSE IF(strngs(1).EQ. '&SOLVENT' ) THEN
         CALL read_solvent(err_open,err_fnf,err_args,err_unr,err_end)
c=======================================================================

c=====Environment ANALYSIS==============================================
      ELSE IF(strngs(1).EQ. '&ANALYSIS' ) THEN
         CALL read_analysis(err_args,err_unr,err_end)

c=======================================================================

c=====Environment PROPERTIES ===========================================
      ELSE IF(strngs(1).EQ. '&PROPERTIES') THEN
         CALL read_properties(fxrms,favg,favg_xrms,fvi,gofr_fprint
     &        ,fcavities,ftop_print,gofr_favg,gofr_fcomp
     &        ,fprtvaf,ftotvaf,fnovaf,fdipole,fnative,ffragm_dist
     &        ,fhbonds,fhhisto,frms,fgyr,freq_ef,freq_dp,fxslt
     &        ,finst_fit,fcalc_cofm,finst_lda,flda_flu,flda_hyd
     &        ,fprot_hyd,fprot_lda,sofk_fprint,sofk_fcomp,err_open
     &        ,err_args,err_unr,err_end)
c=======================================================================

c=====Blank Line =======================================================
      ELSE IF(strngs(1).EQ. ' ' ) THEN
         CONTINUE

c=====Unrecognized Environment==========================================
      ELSE
         errmsg='Unrecognized ENVIRONMENT  ---> '/ / strngs(1)
         call xerror(errmsg,80,1,2) 
      END IF
      
      GOTO 100

*-----------------------------------------------------------------------

600   CONTINUE

      WRITE(kprint,17)
      CALL verify_input(fmaxstp,fscale,fprint,fsave,fplot,fascii
     &     ,fprop,fplot_fragm,fplot_center,frject,fconf,fmaxrun,fupdte
     &     ,fxrms,fvi,favg,favg_xrms,fcavities,ftop_print
     &     ,gofr_fprint,gofr_favg,gofr_fcomp,sofk_fprint,sofk_fcomp
     &     ,fprtvaf,ftotvaf,fnovaf,fdipole,fnative,ffragm_dist,fhbonds
     &     ,fhhisto,frms,fgyr,fabmd,freq_ef,freq_dp,fxslt,finst_fit
     &     ,fcalc_cofm,finst_lda,flda_flu,flda_hyd,fprot_hyd,fprot_lda
     &     ,iret,errmsg)

      IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)

      RETURN
      
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

1     FORMAT(80('='))
2     FORMAT('=',78x,'=')
3     FORMAT
     x    ('=               oooooo       rrrrrrr        aaaaa    ',
     x     '      ccccccc             ='/
     x     '=             oooooooooo     rrrrrrrr     aaaaaaaaa  ',
     x     '     ccccccccc            ='/
     x     '=            ooo      ooo    rr    rrr    aaa   aaa  ',
     x     '    ccc                   ='/
     x     '=            oo        oo    rr    rr    aaa     aaa ',
     x     '    cc                    ='/
     x     '=            oo        oo    rr rrr      aaaaaaaaaaa ',
     x     '    cc                    ='/
     x     '=            ooo      ooo    rr  rrr     aaaaaaaaaaa ',
     x     '    ccc                   ='/
     x     '=             oooooooooo     rr    rr    aa       aa ',
     x     '     ccccccccc            ='/
     x     '=               oooooo       rr    rrr   aa       aa ',
     x     '      ccccccc             =')
4     FORMAT(
     x'=                              (Version 5.0.0)',
     x'                                 =')
5     FORMAT(
     x'=     "A Molecular Dynamics Program to Simulate',
     x' Complex Molecular Systems"     =')
6     FORMAT(
     x'=                     Copyright(C) 1989 - 1999   ',
     x'                              =')
7     FORMAT(
     x'=                        All Right Reserved',
     x'                                    =')
10    FORMAT(
     x'=   ORAC is provided "as is" and without any',
     x' warranty express or implied.      =',/
     x'=   The user assumes all risks of using ORAC.',
     x'                                  =')
11    FORMAT(
     x'=   The user may make copies of ORAC for his',
     x'/her own use, and modify those     =',/
     x'=   copies. The user MAY NOT distribute any copy of the',
     x' original or            =',/
     x'=   modified source code to any users at any sites other',
     x' than his/her own.     =')
12    FORMAT(
     x'=                                                ',
     x'                              =')
13    FORMAT(
     x'=                      Centre  d''Etudes Saclay',
     x'                                 ='/
     x'=                       Gif sur Yvette, FRANCE',
     x'                                 =')
16    FORMAT(
     x     / /'=========================== INPUT TO ',
     x     'THE RUN ==================================='/
     x     '=',78(' '),'=')
17    FORMAT('=',78(' '),'='/80('=')/ /)

      END
