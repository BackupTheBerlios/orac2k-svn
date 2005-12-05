      SUBROUTINE read_potential(fupdte,fabmd,err_args,err_unr,err_end)

************************************************************************
*   Time-stamp: <2005-03-08 16:13:53 marchi>                             *
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

      USE Module_Extra_Forces; USE Module_Thole
      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER iret
      REAL*8  fupdte,fabmd
      CHARACTER*37 err_args(2)
      CHARACTER*20 err_end 
      CHARACTER*27 err_unr(4)

*----------------------- VARIABLES IN COMMON --------------------------*
      
      INCLUDE 'parst.h'
      INCLUDE 'cpropar.h'
      INCLUDE 'unit.h'

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER nword,nsevere,nwarning,j,n_fix_molecules
      CHARACTER*80 errmsg
      CHARACTER*80 line,strngs(40)
      CHARACTER*8 fmt
      CHARACTER*1 sep(2),comm(2)
      LOGICAL  exist
      REAL*8 dummy
      DATA sep/' ',','/comm/'(',')'/

*----------------------- EXECUTABLE STATEMENTS ------------------------*


c=======================================================================
c     Environment parser starts here 
c=======================================================================

      n_fix_molecules=0
      j=0
      nsevere = 0 
      nwarning = 0 
      line(79:80)='  '
100   READ(knlist,'(a78)',END=600) line(1:78)
      CALL wrenc(kprint,line)
      IF(line(1:1) .EQ. '#') GOTO 100 
      CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
      IF(iret.EQ.1) THEN 
         errmsg='while parsing line: toomany strings'
         CALL xerror(errmsg,80,1,2)
         nsevere = nsevere + 1
         go to 100
      END IF

c==== Command GROUP_CUTOFF=============================================

      IF(strngs(1).EQ. 'GROUP_CUTOFF' ) THEN
         if(nword.eq.2.or.nword.eq.4) THEN 
            IF(strngs(2).EQ. 'ON' ) THEN
               grpcut=.TRUE.
               IF(nword .EQ. 4) THEN
                  IF(strngs(3) .EQ. 'CUTSPH') THEN
                     CALL fndfmt(2,strngs(4),fmt)
                     READ(strngs(4),fmt,err=20) rspon
                  END IF
               ELSE
                  nwarning = nwarning  + 1
                  errmsg='Unrecognized or missing keyword '// 
     &                 '; RSPON set to 1.0' 
                  CALL xerror(errmsg,80,1,11)
                  rspon=1.0D0
               END IF
            ELSE IF(strngs(2).EQ. 'OFF' ) THEN
               grpcut=.FALSE.
            ELSE
               nsevere = nsevere + 1
               errmsg=err_unr(3)//strngs(2) 
               CALL xerror(errmsg,80,1,30)
            END IF
         ELSE
            nsevere = nsevere + 1
            errmsg=err_args(1) //'1'
            CALL xerror(errmsg,80,1,30)
         END IF

c==== Command EWALD====================================================

      ELSE IF(strngs(1).EQ. 'EWALD' ) THEN
         IF(nword .EQ. 1) THEN
            nsevere = nsevere + 1
            errmsg=err_args(1) //'1'
            CALL xerror(errmsg,80,1,30)
         ELSE
            IF(strngs(2).EQ. 'ON' ) THEN
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
               
            ELSE IF(strngs(2).EQ. 'PME' ) THEN
               clewld=.true.
               pme=.true.
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
            ELSE IF(strngs(2) .EQ. 'REMOVE_MOMENTUM') THEN
               remove_momentum=.TRUE.
   
            ELSE IF(strngs(2).EQ. 'OFF' ) THEN
               clewld=.FALSE.
               
            ELSE
               nsevere = nsevere + 1
               errmsg=err_unr(3) // strngs(2) 
               CALL xerror(errmsg,80,1,30)
            END IF
         END IF

c==== Command ERFC_SPLINE==============================================

      ELSE IF(strngs(1).EQ. 'ERFC_SPLINE' ) THEN
         erfc_spline=.TRUE.
         IF(nword .NE. 1) THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) erfc_bin
         END IF
         
c==== Command ERF_CORR=================================================

      ELSE IF(strngs(1).EQ. 'ERF_CORR' ) THEN
         erf_corr=.TRUE.
         IF(nword .EQ. 4) THEN
            CALL fndfmt(1,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) nbinew
            CALL fndfmt(2,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) rlew
            CALL fndfmt(2,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) ruew
            delew = (ruew-rlew)/DFLOAT(nbinew) 
          ELSE
            nsevere = nsevere + 1
            errmsg=err_args(1) // '3 after keyword "ERF_CORR"' 
            CALL xerror(errmsg,80,1,30)
         ENDIF

c==== Command EXTRA_FORCE =============================================

      ELSE IF(strngs(1).EQ. 'EXTRA_FORCE_1' ) THEN
         Extra_Force=.TRUE.
         IF(nword .NE. 5) THEN
            nsevere = nsevere + 1
            errmsg=err_args(1) // '4 after keyword "EXTRA_FORCE"' 
            CALL xerror(errmsg,80,1,30)
         ELSE
            inp_1 % type_a=strngs(2)
            inp_1 % type_b=strngs(3)
            CALL fndfmt(2,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) inp_1 % K
            CALL fndfmt(2,strngs(5),fmt)
            READ(strngs(5),fmt,err=20) inp_1 % r0
         END IF

      ELSE IF(strngs(1).EQ. 'EXTRA_FORCE_2' ) THEN
         Extra_Force=.TRUE.
         IF(nword .NE. 5) THEN
            nsevere = nsevere + 1
            errmsg=err_args(1) // '4 after keyword "EXTRA_FORCE"' 
            CALL xerror(errmsg,80,1,30)
         ELSE
            inp_2 % type_a=strngs(2)
            inp_2 % type_b=strngs(3)
            CALL fndfmt(2,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) inp_2 % K
            CALL fndfmt(2,strngs(5),fmt)
            READ(strngs(5),fmt,err=20) inp_2 % r0
         END IF

      ELSE IF(strngs(1).EQ. 'EXTRA_FORCE_3' ) THEN
         Extra_Force=.TRUE.
         IF(nword .NE. 5) THEN
            nsevere = nsevere + 1
            errmsg=err_args(1) // '4 after keyword "EXTRA_FORCE"' 
            CALL xerror(errmsg,80,1,30)
         ELSE
            inp_3 % type_a=strngs(2)
            inp_3 % type_b=strngs(3)
            CALL fndfmt(2,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) inp_3 % K
            CALL fndfmt(2,strngs(5),fmt)
            READ(strngs(5),fmt,err=20) inp_3 % r0
         END IF


c==== Command CUTOFF===================================================

      ELSE IF(strngs(1).EQ. 'CUTOFF' ) THEN
         IF(nword.ne.1) THEN 
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) rspoff
            CALL fndfmt(2,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) rspon
         END IF
         
c==== Command LINKED CELL =============================================

      ELSE IF(strngs(1).EQ. 'LINKED_CELL' ) THEN
         LINKED_CELL = .TRUE.            
         IF(nword.EQ.4) THEN 
            CALL fndfmt(1,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) ncx
            CALL fndfmt(1,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) ncy
            CALL fndfmt(1,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) ncz
         ELSE IF(nword .EQ. 5) THEN
            CALL fndfmt(1,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) ncx
            CALL fndfmt(1,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) ncy
            CALL fndfmt(1,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) ncz
            CALL fndfmt(1,strngs(5),fmt)
            READ(strngs(5),fmt,err=20) nupdte_index
         ELSE
            nsevere = nsevere + 1
            errmsg=err_args(1) //'3'
            CALL xerror(errmsg,80,1,30)
         ENDIF

c==== Command VERLET_LIST==============================================

      ELSE IF(strngs(1).EQ. 'VERLET_LIST' ) THEN
         LINKED_CELL = .FALSE.            

c==== Command UPDATE===================================================

      ELSE IF(strngs(1).EQ. 'UPDATE' ) THEN
         IF(nword .EQ. 3) THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) fupdte
            CALL fndfmt(2,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) rspcut
         ELSE IF(nword .EQ. 5) THEN
            CALL fndfmt(2,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) hrcut
            CALL fndfmt(2,strngs(5),fmt)
            READ(strngs(5),fmt,err=20) hacut
         ELSE
            nsevere = nsevere + 1
            errmsg=err_args(1) //'4'
            CALL xerror(errmsg,80,1,30)
         END IF

c==== Command SAMPLE_ENERGY============================UNSUPPORTED=====

      ELSE IF(strngs(1).EQ. 'SAMPLE_ENERGY' ) THEN
         errmsg = 'UNSUPPORTED COMMAND'
         nwarning = nwarning + 1
         call xerror(errmsg,80,1,11)
         lenerg=.TRUE.
         if(nword.eq.4) THEN 
            CALL fndfmt(1,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) menerg
            IF(strngs(3) .EQ. 'OPEN') THEN
               CALL uscrpl(strngs(4),80)
               INQUIRE(FILE=strngs(4),EXIST=exist)
               IF(exist) THEN
                  CALL openf(kenerg,strngs(4),'FORMATTED','OLD',0)
               ELSE
                  CALL openf(kenerg,strngs(4),'FORMATTED','NEW',0)
               END IF
            ELSE
               nsevere = nsevere + 1
               errmsg='OPEN keyword not found.'
               CALL xerror(errmsg,80,1,30)
            END IF
         ELSE
            nsevere = nsevere + 1
            errmsg=err_args(1) //'3'
            CALL xerror(errmsg,80,1,30)
         END IF

c==== Command U_ENERGY===================================UNSUPPORTED===

      ELSE IF(strngs(1).EQ. 'U_ENERGY' ) THEN
         errmsg = 'UNSUPPORTED COMMAND'
         nwarning = nwarning + 1
         call xerror(errmsg,80,1,11)
         uenerg=.TRUE.
         CALL fndfmt(2,strngs(2),fmt)
         READ(strngs(2),fmt,err=20) uealfa
         CALL fndfmt(2,strngs(3),fmt)
         READ(strngs(3),fmt,err=20) uemin
         CALL fndfmt(2,strngs(4),fmt)
         READ(strngs(4),fmt,err=20) uemax

c==== Command  HBOND =================================================

      ELSEIF(strngs(1).EQ. 'H-BOND' ) THEN
         IF(strngs(2).EQ. 'ON' ) THEN
            hydbnd=.TRUE.
            CALL fndfmt(2,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) hrson
            CALL fndfmt(2,strngs(4),fmt)
            READ(strngs(4),fmt,err=20) hrsoff
            CALL fndfmt(2,strngs(5),fmt)
            READ(strngs(5),fmt,err=20) hanon
            CALL fndfmt(2,strngs(6),fmt)
            READ(strngs(6),fmt,err=20) hanoff
            CALL fndfmt(1,strngs(7),fmt)
            READ(strngs(7),fmt,err=20) nhskip
         ELSE IF(strngs(2).EQ. 'OFF' ) THEN
            hydbnd=.FALSE.
         ELSE
            errmsg=err_unr(3) // strngs(2) 
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF

c==== Command  COMMAND FIX_MOLECULE ==================================

      ELSE IF(strngs(1).EQ. 'HEAVY_MOLECULE' ) THEN
         pfix=.TRUE.
         mass_pfix=1.0D+30
800      READ(knlist,'(a78)',END=900) line(1:78)
         CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 800
         CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
         
c---- Subcommand define  ------------------------------------------------

         IF(strngs(1) .EQ. 'define' ) THEN
            CALL parse_numbers(err_unr,strngs,nword,prot_fix
     &           ,n_fix_molecules,nsevere)

c---- Subcommand mass  --------------------------------------------------

         ELSE IF(strngs(1) .EQ. 'mass' ) THEN
            IF(nword .NE. 2) THEN
               errmsg=err_args(2)//'1'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            ELSE
               CALL fndfmt(2,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) mass_pfix
            END IF

         ELSE IF(strngs(1).EQ. ' ') THEN
            CONTINUE

         ELSE IF(strngs(1) .EQ. 'END') THEN
            GOTO 100
         ELSE
c---        could not fine SUBCOMMAND of END
            errmsg=err_unr(2) // strngs(1)// ' or missing END'
            call xerror(errmsg,80,1,30)
            nsevere = nsevere + 1 
         END IF

         GOTO 800

c==== Command  I-TORSION==============================================

      ELSE IF(strngs(1).EQ. 'I-TORSION' ) THEN
         IF( nword .NE. 2) THEN
            errmsg=err_args(1)//'1'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         ELSE
            IF(strngs(2) .EQ. 'HARMONIC' ) THEN
               itor_ptype=1
            ELSE IF(strngs(2) .EQ. 'COSINE' ) THEN
               itor_ptype=2
            ELSE
               errmsg=err_unr(3) // strngs(2)
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF
         END IF
         
c==== Command  H-MASS=================================================

      ELSE IF(strngs(1).EQ. 'H-MASS' ) THEN
         hmass=.TRUE.
         CALL fndfmt(2,strngs(2),fmt)
         READ(strngs(2),fmt,err=20) hdmass
         
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
         
c==== Command  LJ-FUDGE===============================================

      ELSE IF(strngs(1) .EQ. 'LJ-FUDGE') THEN
         if(nword.eq.2) THEN 
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) lj_fudge
            lj_fudgeb=lj_fudge
         ELSE
            errmsg=err_args(1)//'1'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         
c==== Command  BENDING================================================

      ELSE IF(strngs(1) .EQ. 'BENDING') THEN 
         if(nword.eq.2) THEN 
            IF(strngs(2) .EQ. 'OFF') THEN
               bending=.FALSE.
            ELSE IF(strngs(2) .EQ. 'ON') THEN
               bending=.TRUE.
            ELSE 
               errmsg=err_unr(3)//strngs(2)
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF
         ELSE
            errmsg=err_args(1)//'1'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
            
c==== Command  WIDOM======================================UNSUPPORTED=

      ELSE IF(strngs(1) .EQ. 'WIDOM' ) THEN
         errmsg = err_unr(4)   
         call xerror(errmsg,80,1,11)
         widom=.TRUE.
         CALL fndfmt(1,strngs(2),fmt)
         READ(strngs(2),fmt,err=20) wres
         CALL fndfmt(2,strngs(3),fmt)
         READ(strngs(3),fmt,err=20) wcoef
         
c==== Command  JORGENSEN==================================UNSUPPORTED=

      ELSE IF(strngs(1) .EQ. 'JORGENSEN') THEN
         errmsg = err_unr(4)   
         call xerror(errmsg,80,1,11)
         amphi=.TRUE.
         
c==== Command  AUTO_DIHEDRAL==========================================

      ELSE IF(strngs(1).EQ. 'AUTO_DIHEDRAL' ) THEN
         adihed = .TRUE.
         
c==== Command  SELECT_DIHEDRAL========================================

      ELSE IF(strngs(1).EQ. 'SELECT_DIHEDRAL' ) THEN
         adihed = .FALSE.
         
c==== Command  SELECT_STRETCHING======================================

      ELSE IF(strngs(1).EQ. 'STRETCHING' ) THEN
         stretch = .TRUE.
         IF(nword.gt.1) THEN 
            IF(strngs(2) .EQ. 'HEAVY') THEN
               stretch_heavy= .TRUE.
            ELSE 
               errmsg=err_unr(3)//strngs(2)
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF
         END IF

c==== Command  COMMAND ABMD ==========================================

      ELSE IF(strngs(1) .EQ. 'ABMD' ) THEN
         abmd=.TRUE.
c------- read the line
200      READ(knlist,'(a78)',END=300) line(1:78)
         CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 200
         CALL parse(line,sep,2,comm,strngs,40,nword,
     x        iret,errmsg)
         
c---- Subcommand fold ------------------------------------------------

         IF(strngs(1).EQ. 'fold' ) THEN
            fold=.TRUE.
            rspset=-1.0D0

c---- Subcommand fold ------------------------------------------------

         ELSE IF(strngs(1).EQ. 'unfold' ) THEN
            fold=.TRUE.
            rspset=1.0D0

c---- Subcommand unbiased --------------------------------------------

         ELSE IF(strngs(1).EQ. 'unbiased' ) THEN
            abmd_unbias=.TRUE.
            IF(nword .EQ. 1) THEN
               rspset=0.0D0
            ELSE IF(nword .EQ. 2) THEN
               CALL fndfmt(2,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) rspset
            END IF

c---- Subcommand dissociate ------------------------------------------

         ELSE IF(strngs(1).EQ. 'dissociate' ) THEN
            dissociate=.TRUE.
            IF(strngs(2) .EQ. 'molecule') THEN
               diss_mol=.TRUE.
               CALL fndfmt(1,strngs(3),fmt)
               READ(strngs(3),fmt,err=20) mol_diss(1)
               CALL fndfmt(1,strngs(4),fmt)
               READ(strngs(4),fmt,err=20) mol_diss(2)
            ELSE IF(strngs(2) .EQ. 'atoms') THEN
               diss_atoms=.TRUE.
               CALL fndfmt(1,strngs(3),fmt)
               READ(strngs(3),fmt,err=20) atoms_diss(1,1)
               CALL fndfmt(1,strngs(4),fmt)
               READ(strngs(4),fmt,err=20) atoms_diss(2,1)
               CALL fndfmt(1,strngs(5),fmt)
               READ(strngs(5),fmt,err=20) atoms_diss(1,2)
               CALL fndfmt(1,strngs(6),fmt)
               READ(strngs(6),fmt,err=20) atoms_diss(2,2)
            ELSE
               errmsg=err_unr(3) // strngs(2)//' subcommand dissociate'
               call xerror(errmsg,80,1,30)
               nsevere = nsevere + 1 
            END IF

c---- Subcommand associate -------------------------------------------

         ELSE IF(strngs(1).EQ. 'associate' ) THEN
            associate=.TRUE.
            IF(strngs(2) .EQ. 'molecule') THEN
               diss_mol=.TRUE.
               CALL fndfmt(1,strngs(3),fmt)
               READ(strngs(3),fmt,err=20) mol_diss(1)
               CALL fndfmt(1,strngs(4),fmt)
               READ(strngs(4),fmt,err=20) mol_diss(2)
            ELSE IF(strngs(2) .EQ. 'atoms') THEN
               diss_atoms=.TRUE.
               CALL fndfmt(1,strngs(3),fmt)
               READ(strngs(3),fmt,err=20) atoms_diss(1,1)
               CALL fndfmt(1,strngs(4),fmt)
               READ(strngs(4),fmt,err=20) atoms_diss(2,1)
               CALL fndfmt(1,strngs(5),fmt)
               READ(strngs(5),fmt,err=20) atoms_diss(1,2)
               CALL fndfmt(1,strngs(6),fmt)
               READ(strngs(6),fmt,err=20) atoms_diss(2,2)
            ELSE
               errmsg=err_unr(3) // strngs(2)//' subcommand dissociate'
               call xerror(errmsg,80,1,30)
               nsevere = nsevere + 1 
            END IF

c---- Subcommand fold_native -----------------------------------------

         ELSE IF(strngs(1).EQ. 'native_fold' ) THEN
            abmd_native=.TRUE.
            abmd_nat_dir=1.0D0

c---- Subcommand unfold_native ---------------------------------------

         ELSE IF(strngs(1).EQ. 'native_unfold' ) THEN
            abmd_native=.TRUE.
            abmd_nat_dir=-1.0D0

c---- Subcommand native_theta ----------------------------------------

         ELSE IF(strngs(1).EQ. 'native_para' ) THEN
            abmd_native=.TRUE.
            IF(nword .NE. 3) THEN
               errmsg=err_args(2)//'2'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            ELSE
               CALL fndfmt(2,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) native_theta
               CALL fndfmt(2,strngs(3),fmt)
               READ(strngs(3),fmt,err=20) native_dist
            END IF

c==== Subcommand  native_template ======================================

         ELSE IF(strngs(1).EQ. 'native_template') THEN
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


c---- Subcommand torsion ---------------------------------------------

         ELSE IF(strngs(1).EQ. 'torsion' ) THEN
            IF(nword .NE. 3) THEN
               nsevere = nsevere + 1
               errmsg=err_args(1) //'2'
               CALL xerror(errmsg,80,1,30)
            END IF
            abmd_tors=.TRUE.
            IF(strngs(2) .EQ. 'decrease') THEN
               abmd_tors_dir=-1.0D0
            ELSE IF(strngs(2) .EQ. 'increase') THEN
               abmd_tors_dir=1.0D0
            ELSE
               nsevere = nsevere + 1
               errmsg=err_unr(3)//strngs(2) 
               CALL xerror(errmsg,80,1,30)
            END IF
            CALL fndfmt(1,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) abmd_tors_num


c---- Subcommand print -----------------------------------------------

         ELSE IF(strngs(1).EQ. 'print' ) THEN
            IF(nword.eq.4) THEN 
               CALL fndfmt(2,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) fabmd
               IF(strngs(3) .EQ. 'OPEN') THEN
                  CALL uscrpl(strngs(4),80)
                  INQUIRE(FILE=strngs(4),EXIST=exist)
                  IF(exist) THEN
                     CALL openf(kabmd,strngs(4),'FORMATTED','OLD',0)
                  ELSE
                     CALL openf(kabmd,strngs(4),'FORMATTED','NEW',0)
                  END IF
                  WRITE(kabmd,1000)
               ELSE
                  nsevere = nsevere + 1
                  errmsg='OPEN keyword not found.'
                  CALL xerror(errmsg,80,1,30)
               END IF
            ELSE
               nsevere = nsevere + 1
               errmsg=err_args(1) //'3'
               CALL xerror(errmsg,80,1,30)
            END IF

c---- Subcommand crystalize ------------------------------------------

         ELSE IF(strngs(1).EQ. 'crystalize' ) THEN
            abmd_cryst=.TRUE.
            abmd_cryst_dir=1.0D0
            IF(strngs(2) .EQ. 'K-vector') THEN
               CALL fndfmt(2,strngs(3),fmt)
               READ(strngs(3),fmt,err=20) abmd_kvect
            ELSE
               errmsg=err_unr(3) // strngs(2)//' subcommand crystalize'
               call xerror(errmsg,80,1,30)
               nsevere = nsevere + 1 
            END IF

c---- Subcommand melt ------------------------------------------------

         ELSE IF(strngs(1).EQ. 'melt' ) THEN
            abmd_cryst=.TRUE.
            abmd_cryst_dir=-1.0D0
            IF(strngs(2) .EQ. 'K-vector') THEN
               CALL fndfmt(2,strngs(3),fmt)
               READ(strngs(3),fmt,err=20) abmd_kvect
            ELSE
               errmsg=err_unr(3) // strngs(2)//' subcommand melt'
               call xerror(errmsg,80,1,30)
               nsevere = nsevere + 1 
            END IF

c---- Subcommand spring ----------------------------------------------

         ELSE IF(strngs(1).EQ. 'spring' ) THEN
            IF(nword .NE. 2) THEN
               errmsg=err_args(2)//'1'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            ELSE
               CALL fndfmt(2,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) spring
            END IF

         ELSE IF(strngs(1).EQ. ' ') THEN
            CONTINUE

         ELSE IF(strngs(1) .EQ. 'END') THEN
            GOTO 100
         ELSE
c---        could not fine SUBCOMMAND of END
            errmsg=err_unr(2) // strngs(1)// ' or missing END'
            call xerror(errmsg,80,1,30)
            nsevere = nsevere + 1 
         END IF

         GOTO 200

c==== Command  COMMAND NONBONDED_OFF==================================

      ELSE IF(strngs(1).EQ. 'NONBONDED_OFF' ) THEN
         nonbnd=.FALSE.
         
c==== Command  CONSTRAINT  ===========================================

      ELSE IF(strngs(1).EQ. 'CONSTRAINT') THEN
         IF(nword .GT. 1) THEN 
            IF(strngs(2) .EQ. 'SHAKE') THEN
               mim_lim=0
            ELSE IF(strngs(2) .EQ. 'MIM') THEN
               IF(nword .EQ. 3) THEN
                  CALL fndfmt(1,strngs(3),fmt)
                  READ(strngs(3),fmt,err=20) mim_lim
               ELSE
                  mim_lim=20
               END IF
            ELSE
               nsevere = nsevere + 1
               errmsg=err_unr(3)//strngs(2) 
               CALL xerror(errmsg,80,1,30)
            END IF
         ELSE
            nsevere = nsevere + 1
            errmsg=err_args(1) //'1'
            CALL xerror(errmsg,80,1,30)
         END IF

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

c==== Command  OLD-DIPOLES ===========================================

            ELSE IF(strngs(1).EQ.'OLD-DIPOLES') THEN
               Old_dipoles=.TRUE.

c==== Command  MESOS =================================================

            ELSE IF(strngs(1).EQ.'MESOS') THEN
               mesos=.TRUE.
               if(nword.eq.2) THEN
                  CALL fndfmt(2,strngs(2),fmt)
                  READ(strngs(2),fmt,err=20) dummy
                  mesos_rho=1.0D0/dummy
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
c==== Command  THOLE =================================================

            ELSE IF(strngs(1) .EQ. 'THOLE') THEN
               Thole=.TRUE.

c==== Command  EXTERNAL ==============================================

            ELSE IF(strngs(1) .EQ. 'MODEL') THEN
               IF(nword .EQ. 2) THEN
                  IF(strngs(2) .EQ. 'GAUSS') THEN
                     What_To_do_Pol='Gauss'
                  ELSE IF(strngs(2) .EQ. 'FULL') THEN
                     What_To_do_Pol='Full'
                  ELSE IF(strngs(2) .EQ. 'DIRECT') THEN
                     What_To_do_Pol='Direct'
                  ELSE
                     errmsg=err_unr(3)//strngs(2)
                     CALL xerror(errmsg,80,1,30)
                     nsevere = nsevere + 1
                  END IF
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

c==== subcommand CUTOFF =============================================

            ELSE IF(strngs(1).EQ. 'CUTOFF_EL' ) THEN
               IF(nword.ne.1) THEN
                  CALL fndfmt(2,strngs(2),fmt)
                  READ(strngs(2),fmt,err=20) rcut_EL
               END IF

c==== subcommand CUTOFF =============================================

            ELSE IF(strngs(1).EQ. 'CUTOFF_LJ' ) THEN
               IF(nword.ne.1) THEN
                  CALL fndfmt(2,strngs(2),fmt)
                  READ(strngs(2),fmt,err=20) rcut_LJ
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

c==== Command  KEEP_BONDS =============================================

      ELSE IF(strngs(1).EQ. 'KEEP_BONDS' ) THEN
         readjust_cnstr = .FALSE.

c==== Command  VIRTUAL ================================================

      ELSE IF(strngs(1).EQ. 'VIRTUAL' ) THEN
         Virtual_residue = .TRUE.
         Virtual_types=0
400      READ(knlist,'(a78)',END=500) line(1:78)
         CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 400
         CALL parse(line,sep,2,comm,strngs,40,nword,
     x        iret,errmsg)

         IF(strngs(1) .EQ. 'residue' .AND. nword .EQ. 2) THEN
            Residue_virt(1:8)=strngs(2)(1:8)
            
         ELSE IF(strngs(1).EQ. 'r_min' .AND. nword .EQ. 3) THEN
            Virtual_types=Virtual_types+1
            types_virt(Virtual_types)(1:7)=strngs(2)(1:7)
            CALL fndfmt(2,strngs(3),fmt)
            READ(strngs(3),fmt,err=20) r_min_virt(Virtual_types)

         ELSE IF(strngs(1).EQ. ' ') THEN
            CONTINUE

         ELSE IF(strngs(1) .EQ. 'END') THEN
            GOTO 100
         ELSE
c---        could not fine SUBCOMMAND of END
            errmsg=err_unr(2) // strngs(1)// ' or missing END'
            call xerror(errmsg,80,1,30)
            nsevere = nsevere + 1 
         END IF

         GOTO 400

c==== Command  ADJUST_BONDS ==========================================

      ELSE IF(strngs(1).EQ. 'ADJUST_BONDS' ) THEN
         readjust_cnstr = .TRUE.

c==== Command BLANK LINE===============================================

      ELSE IF(strngs(1).EQ. ' ') THEN

c==== Begininning of next ENVIRONMENT =================================

      ELSE IF(strngs(1)(1:1).EQ. '&'.AND.strngs(1).NE. '&END') THEN
         errmsg= err_unr(1) // strngs(1)(1:8) // err_end
         CALL xerror(errmsg,80,1,30)
         nsevere = nsevere + 1
         GO TO 600

c==== Command &END ====================================================

      ELSE IF(strngs(1).EQ. '&END') THEN
         GOTO 600
         
      ELSE
         errmsg= err_unr(1) // strngs(1)(1:8) // err_end
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
         call int_str(nsevere,fmt,j)
         errmsg=fmt(1:j) //' ERRORS WHILE EXECUTING READ_POTENTIAL'
         call xerror(errmsg,80,1,2)
         STOP
      ELSE IF(nsevere.gt.99) THEN 
         errmsg= 'MORE THAN 99 ERRORS WHILE EXECUTING READ_POTENTIAL'
         call xerror(errmsg,80,1,2)
         STOP
      END IF
      if(nwarning.gt.0.and.nwarning.lt.99) then 
         j=0
         call int_str(nwarning,fmt,j)
         errmsg= fmt(1:j)//' WARNINGS WHILE EXECUTING READ_POTENTIAL'
         CALL xerror(errmsg,80,1,1)
      ELSE IF(nwarning.gt.99) THEN 
         errmsg= 'MORE THAN 99 WARNINGS WHILE EXECUTING READ_POTENTIAL'
         call xerror(errmsg,80,1,1)
      ENDIF    

c==============================================================================
c     Verification Part
c==============================================================================

c--   Check if pme_orded is gt 3 
      if(pme.and.(pme_order.lt.3)) THEN 
         errmsg='"Ewald": PME order must be at least 3'
         CALL xerror(errmsg,80,1,30)
         nsevere = nsevere + 1
      END IF
      if(linked_cell) THEN 
         if(ncx.le.0.or.ncy.le.0.or.ncz.le.0) THEN 
            errmsg='Bad Linked Cell Neighbor list parameters'  
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
      END IF

      if(nsevere.gt.0.and.nsevere.lt.99) then 
         j=0
         call int_str(nsevere,fmt,j)
         errmsg=fmt(1:j) //' ERRORS WHILE EXECUTING READ_POTENTIAL'
         CALL xerror(errmsg,80,1,2)
         STOP 
      ELSE IF(nsevere.gt.99) THEN
         errmsg='MORE THAN 99 ERRORS WHILE EXECUTING READ_POTENTIAL'
         call xerror(errmsg,80,1,2)
         STOP
      END IF
1000  FORMAT('#       T          ABMD_a        ABMD_b        V_ab  '
     &     ,'  NonBond ')

      prot_fix(1)=n_fix_molecules
      nprot_fix=n_fix_molecules

      RETURN

c==============================================================================
c     EOF found while reading ABMD subcommands
c==============================================================================

300   CONTINUE
      iret=1
      errmsg='EOF found while reading ABMD subcommands'
      CALL xerror(errmsg,80,1,2)
      STOP

500   CONTINUE
      iret=1
      errmsg='EOF found while reading VIRTUAL subcommands'
      CALL xerror(errmsg,80,1,2)
      STOP

900   CONTINUE
      iret=1
      errmsg='EOF found while reading HEAVY_MOLECULE subcommands'
      CALL xerror(errmsg,80,1,2)
      STOP

c==============================================================================
c     Errors were found
c==============================================================================


 20   CONTINUE
      iret=1
      errmsg='internal reading error: wrong format?? TAB character??'
      CALL xerror(errmsg,80,1,2)
      STOP

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      END
