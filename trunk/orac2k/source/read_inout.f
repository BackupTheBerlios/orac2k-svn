      SUBROUTINE read_inout(fconf,fsave,fplot,fplot_fragm,fplot_center
     &     ,fascii,err_open,err_args,err_end,err_unr)

************************************************************************
*   Time-stamp: <2006-04-25 12:26:19 marchi>                             *
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

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER iret
      CHARACTER*22 err_open
      CHARACTER*37 err_args(3)
      CHARACTER*20 err_end 
      CHARACTER*27 err_unr(4)
      REAL*8 fplot_center,fconf,fsave,fplot,fascii,fplot_fragm

*----------------------- VARIABLES IN COMMON --------------------------*
      
      INCLUDE 'parst.h'
      INCLUDE 'cpropar.h'
      INCLUDE 'unit.h'

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER nword,j,nsevere,nwarning
      CHARACTER*80 errmsg
      CHARACTER*80 line,strngs(40)
      CHARACTER*8 fmt
      CHARACTER*1 sep(2),comm(2)
      LOGICAL  exist
      DATA sep/' ',','/comm/'(',')'/

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      nsevere=0
      nwarning=0
      line(79:80)='  '

c=======================================================================
c     Environment parser starts here 
c=======================================================================

100   READ(knlist,'(a78)',END=600) line(1:78)
      CALL wrenc(kprint,line)
      IF(line(1:1) .EQ. '#') GOTO 100 
      CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
      IF(iret.EQ.1) THEN 
         errmsg='while parsing line: toomany strings'
         CALL xerror(errmsg,80,1,30)
         nsevere = nsevere + 1
         go to 100
      END IF

c==== Complex command  DUMP ===========================================

      IF(strngs(1).EQ. 'DUMP' ) THEN
         dmprnd_o=.TRUE.
200      READ(knlist,'(a78)',END=600) line(1:78)
         CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 200
         CALL parse(line,sep,2,comm,strngs,40,nword,
     x        iret,errmsg)

c------- Subcommand occupy

         IF(strngs(1) .EQ. 'occupy') THEN
            occupy_space=.TRUE.

c------- Subcommand divide_records

         ELSE IF(strngs(1) .EQ. 'atom_record') THEN
            IF(nword .NE. 2) THEN
               errmsg=err_args(3)//'2'
               CALL xerror(errmsg,80,1,30)
               nsevere=nsevere+1
            END IF
            CALL fndfmt(1,strngs(2),fmt) 
            READ(strngs(2),fmt,err=20) atom_record

c------- Subcommand write

         ELSE IF(strngs(1) .EQ. 'write' .OR. strngs(1) .EQ. 'read') THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) fconf
            IF(strngs(3) .EQ. 'OPEN') THEN
               CALL uscrpl(strngs(4),80)
               dmpfil_o=strngs(4)
               INQUIRE(FILE=dmpfil_o,EXIST=exist)
               IF(.NOT. exist) THEN
                  errmsg=
     &                 'Auxiliary File does not exist. Abort.'
                  CALL xerror(errmsg,80,1,30)
                  nsevere=nsevere+1
               END IF
            ELSE
               errmsg=err_open
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF

         ELSE IF(strngs(1) .EQ. ' ') THEN
            GOTO 200

         ELSE IF(strngs(1).EQ. 'END' ) THEN
            GOTO 100
         ELSE
            errmsg=err_unr(3)//strngs(2)//err_end(1:14)//err_end(16:20)
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         GOTO 200
         
c==== Command  TRAJECTORY ===============================================

      ELSE IF(strngs(1).EQ. 'TRAJECTORY' ) THEN
         dmprnd_i=.TRUE.
         CALL uscrpl(strngs(2),80)
         dmpfil_i=strngs(2)
         INQUIRE(FILE=dmpfil_i,EXIST=exist)
         IF(.NOT. exist) THEN
            errmsg=
     &           'Auxiliary File does not exist. Abort.'
            CALL xerror(errmsg,80,1,30)
            nsevere=nsevere+1
         END IF
         
c==== Command  RESTART  ==================================================

      ELSE IF(strngs(1) .EQ. 'RESTART') THEN
300      READ(knlist,'(a78)',END=600) line(1:78)
         CALL wrenc(kprint,line)
         IF(line(1:1) .EQ. '#') GOTO 300
         CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
         
c-------- Subcommand continue
         
         IF(strngs(1) .EQ. 'old_restart') THEN
            restart_old=.TRUE.
            GOTO 300
         
c-------- Subcommand read
         
         ELSE IF(strngs(1) .EQ. 'read') THEN
            restart_read=.TRUE.
            IF(nword .EQ. 2) THEN
               CALL uscrpl(strngs(2),80)
               INQUIRE(FILE=strngs(2),EXIST=exist)
               restart_in=strngs(2)
               IF(.NOT. exist) THEN
                  errmsg='Restart file was not found.'
                  CALL xerror(errmsg,80,1,30)
                  nsevere=nsevere+1
               ELSE
                  CALL openf(kdump_in,strngs(2),'FORMATTED','OLD',0)
                  CLOSE(kdump_in)
               END IF
            ELSE
               errmsg=err_args(3)//'1'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF
            GOTO 300

c-------- Subcommand write
         
         ELSE IF(strngs(1) .EQ. 'write') THEN
            restart_write=.TRUE.
            IF(nword .EQ. 4) THEN
               CALL fndfmt(2,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) fsave
               IF(strngs(3) .EQ. 'OPEN') THEN
                  CALL uscrpl(strngs(4),80)
                  restart_out=strngs(4)
                  INQUIRE(FILE=strngs(4),EXIST=exist)
                  IF(exist) THEN
                     CALL openf(kdump_out,strngs(4),'FORMATTED','OLD',0)
                  ELSE
                     CALL openf(kdump_out,strngs(4),'FORMATTED','NEW',0)
                  END IF
                  CLOSE(kdump_out)

               ELSE
                  errmsg=err_open
                  CALL xerror(errmsg,80,1,30)
                  nsevere = nsevere + 1
               END IF
            ELSE
               errmsg=err_args(3)//'3'
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF
            GOTO 300

         ELSE IF(strngs(1) .EQ. ' ') THEN
            GOTO 300
            
         ELSE IF(strngs(1).EQ. 'END' ) THEN
            GOTO 100
         ELSE
            errmsg=err_unr(3)//strngs(2)//err_end(1:14)//err_end(16:20)
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         GOTO 300
         
c==== Command  PLOT     ==================================================

      ELSE IF(strngs(1).EQ. 'PLOT' ) THEN
         IF(nword.eq.4) THEN 
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) fplot
            IF(strngs(3) .EQ. 'OPEN') THEN
               CALL uscrpl(strngs(4),80)
               INQUIRE(FILE=strngs(4),EXIST=exist)
               IF(exist) THEN
                  CALL openf(kplot,strngs(4),'FORMATTED','OLD',0)
               ELSE
                  CALL openf(kplot,strngs(4),'FORMATTED','NEW',0)
               END IF
            ELSE
               errmsg=err_open
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF
         ELSE IF(nword .EQ. 5) THEN
            IF(strngs(2) .EQ. 'FRAGMENT') THEN
               CALL fndfmt(2,strngs(3),fmt)
               READ(strngs(3),fmt,err=20) fplot_fragm
               IF(strngs(4) .EQ. 'OPEN') THEN
                  CALL uscrpl(strngs(5),80)
                  INQUIRE(FILE=strngs(5),EXIST=exist)
                  IF(exist) THEN
                     CALL openf(kplot_fragm,strngs(5),'FORMATTED','OLD'
     &                    ,0)
                  ELSE
                     CALL openf(kplot_fragm,strngs(5),'FORMATTED','NEW'
     &                    ,0)
                  END IF
               ELSE
                  errmsg=err_open
                  CALL xerror(errmsg,80,1,30)
                  nsevere = nsevere + 1
               END IF
            ELSE IF(strngs(2) .EQ. 'CENTER') THEN
               CALL fndfmt(2,strngs(3),fmt)
               READ(strngs(3),fmt,err=20) fplot_center
               IF(strngs(4) .EQ. 'OPEN') THEN
                  CALL uscrpl(strngs(5),80)
                  INQUIRE(FILE=strngs(5),EXIST=exist)
                  IF(exist) THEN
                     CALL openf(kplot_center,strngs(5),'FORMATTED','OLD'
     &                    ,0)
                  ELSE
                     CALL openf(kplot_center,strngs(5),'FORMATTED','NEW'
     &                    ,0)
                  END IF
               ELSE
                  errmsg=err_open
                  CALL xerror(errmsg,80,1,30)
                  nsevere = nsevere + 1
               END IF
            ELSE
               errmsg= err_unr(3)//strngs(1)(1:8)
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF
         ELSE
            errmsg=err_args(1)//'3'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         ENDIF

c==== Command  ASCII =====================================================

      ELSE IF(strngs(1).EQ. 'ASCII_OUTBOX' ) THEN
         ascii_nocell=.TRUE.
         ascii_wsc=.FALSE.
         CALL fndfmt(2,strngs(2),fmt)
         READ(strngs(2),fmt,err=20) fascii
         IF(strngs(3) .EQ. 'OPEN') THEN
            CALL uscrpl(strngs(4),80)
            INQUIRE(FILE=strngs(4),EXIST=exist)
            IF(exist) THEN
               CALL openf(kout,strngs(4),'FORMATTED','OLD',0)
            ELSE
               CALL openf(kout,strngs(4),'FORMATTED','NEW',0)
            END IF
         ELSE
            errmsg=err_open
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         
c==== Command  ASCII =====================================================

      ELSE IF(strngs(1).EQ. 'ASCII_WSC' ) THEN
         ascii_nocell=.FALSE.
         ascii_wsc=.TRUE.
         CALL fndfmt(2,strngs(2),fmt)
         READ(strngs(2),fmt,err=20) fascii
         IF(strngs(3) .EQ. 'OPEN') THEN
            CALL uscrpl(strngs(4),80)
            INQUIRE(FILE=strngs(4),EXIST=exist)
            IF(exist) THEN
               CALL openf(kout,strngs(4),'FORMATTED','OLD',0)
            ELSE
               CALL openf(kout,strngs(4),'FORMATTED','NEW',0)
            END IF
         ELSE
            errmsg=err_open
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         
c==== Command  ASCII =====================================================

      ELSE IF(strngs(1).EQ. 'ASCII' ) THEN
         ascii_nocell=.FALSE.
         ascii_wsc=.FALSE.
         CALL fndfmt(2,strngs(2),fmt)
         READ(strngs(2),fmt,err=20) fascii
         IF(strngs(3) .EQ. 'OPEN') THEN
            CALL uscrpl(strngs(4),80)
            INQUIRE(FILE=strngs(4),EXIST=exist)
            IF(exist) THEN
               CALL openf(kout,strngs(4),'FORMATTED','OLD',0)
            ELSE
               CALL openf(kout,strngs(4),'FORMATTED','NEW',0)
            END IF
         ELSE
            errmsg=err_open
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF

c==== Command  ASCII =====================================================

      ELSE IF(strngs(1).EQ. 'PLOT_CENTER' ) THEN
         CALL fndfmt(2,strngs(2),fmt)
         READ(strngs(2),fmt,err=20) fascii
         IF(strngs(3) .EQ. 'OPEN') THEN
            CALL uscrpl(strngs(4),80)
            INQUIRE(FILE=strngs(4),EXIST=exist)
            IF(exist) THEN
               CALL openf(kout,strngs(4),'FORMATTED','OLD',0)
            ELSE
               CALL openf(kout,strngs(4),'FORMATTED','NEW',0)
            END IF
         ELSE
            errmsg=err_open
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         

c==== Blank Line =========================================================


      ELSE IF(strngs(1).EQ. ' ') THEN
         CONTINUE

      ELSE IF(strngs(1)(1:1).EQ. '&'.AND.strngs(1).NE. '&END') THEN
         errmsg= err_unr(1)//strngs(1)(1:8)// err_end
         CALL xerror(errmsg,80,1,30)
         nsevere = nsevere + 1
         GO TO 600

c==== &END keyword =======================================================

      ELSE IF(strngs(1).EQ. '&END') THEN
         GOTO 600

c==== No &END found or wrong keyword =====================================
      ELSE
         errmsg= err_unr(1)// strngs(1)(1:8)// err_end
         call xerror(errmsg,80,1,30)
         nsevere = nsevere + 1 
      END IF

      GOTO 100

600   CONTINUE
      
c=======================================================================
c     Environment parser ends here 
c=======================================================================

c--   syntax errors: abort without verifying input 

      if(nsevere.gt.0.and.nsevere.lt.99) then 
         call int_str(nsevere,fmt,j)
         errmsg = fmt(1:j)//' ERRORS WHILE EXECUTING READ_INOUT'
         CALL xerror(errmsg,80,1,2)
      ELSE IF(nsevere.gt.99) THEN 
         errmsg= 'MORE THAN 99 ERRORS WHILE EXECUTING READ_INOUT'
         call xerror(errmsg,80,1,2)
      END IF

      if(nwarning.gt.0.and.nwarning.lt.99) then 
         iret=0
         j=0
         call int_str(nwarning,fmt,j)
         errmsg= fmt(1:j)//'WARNINGS WHILE EXECUTING READ_INOUT'
         CALL xerror(errmsg,80,1,1)
      ELSE IF(nwarning.gt.99) THEN 
         errmsg= 'MORE THAN 99 WARNINGS WHILE EXECUTING READ_INOUT'
         call xerror(errmsg,80,1,1)
      ENDIF    

      RETURN

 20   CONTINUE
      iret=1
      errmsg='internal reading error: wrong format?? Tab Character??' 
      CALL xerror(errmsg,80,1,2)
      STOP
      END

