      SUBROUTINE read_analysis(err_args,err_unr,err_end,err_open)

************************************************************************
*   Time-stamp: <2005-01-31 10:18:37 marchi>                             *
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
      
      CHARACTER*22 err_open
      CHARACTER*37 err_args(2)
      CHARACTER*20 err_end 
      CHARACTER*27 err_unr(4)

*----------------------- VARIABLES IN COMMON --------------------------*
      
      INCLUDE 'parst.h'
      INCLUDE 'cpropar.h'
      INCLUDE 'unit.h'

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER nword,i,j
      CHARACTER*80 line,strngs(40)
      CHARACTER*8 fmt
      CHARACTER*1 sep(2),comm(2)
      LOGICAL  exist
      INTEGER iret,nsevere
      CHARACTER*80 errmsg
      DATA sep/' ',','/comm/'(',')'/

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      analys=.TRUE.
      line(79:80)='  '
100   READ(knlist,'(a78)',END=600) line(1:78)
      CALL wrenc(kprint,line)
      IF(line(1:1) .EQ. '#') GOTO 100 
      CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
      IF(iret.EQ.0) THEN

c==== Command UPDATE===================================================

         IF(strngs(1).EQ. 'UPDATE' ) THEN
            IF(nword .EQ. 3) THEN
               CALL fndfmt(1,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) update_anl
               CALL fndfmt(2,strngs(3),fmt)
               READ(strngs(3),fmt,err=20) rspoff
            ELSE
               nsevere = nsevere + 1
               errmsg=err_args(1) //'2'
               CALL xerror(errmsg,80,1,30)
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

c==== Command START ===================================================

         ELSE IF(strngs(1).EQ. 'START') THEN
            IF(nword .EQ. 2) THEN
               CALL fndfmt(1,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) start_anl
            ELSE
               nsevere = nsevere + 1
               errmsg=err_args(1) //'1'
               CALL xerror(errmsg,80,1,30)
            END IF

c==== Command STOP ====================================================

         ELSE IF(strngs(1).EQ. 'STOP') THEN
            IF(nword .EQ. 2) THEN
               CALL fndfmt(1,strngs(2),fmt)
               READ(strngs(2),fmt,err=20) stop_anl
            ELSE
               nsevere = nsevere + 1
               errmsg=err_args(1) //'1'
               CALL xerror(errmsg,80,1,30)
            END IF

c==== Command SKIP ====================================================

         ELSE IF(strngs(1).EQ. 'SKIP') THEN
            skip_step=.TRUE.
            IF(nword .NE. 1) THEN
               nsevere = nsevere + 1
               errmsg=err_args(1) //'1'
               CALL xerror(errmsg,80,1,30)
            END IF

c==== Command SCAN ====================================================

         ELSE IF(strngs(1).EQ. 'SCAN') THEN
            scan_traj=.TRUE.
            IF(nword .NE. 1) THEN
               nsevere = nsevere + 1
               errmsg=err_args(1) //'1'
               CALL xerror(errmsg,80,1,30)
            END IF

c==== Command CONTINUE ================================================

         ELSE IF(strngs(1).EQ. 'CONTINUE') THEN

c==== Command &END ====================================================
            
         ELSE IF(strngs(1).EQ. '&END') THEN
            GOTO 600
            
         ELSE IF(strngs(1).EQ. ' ') THEN
            CONTINUE
         ELSE
            errmsg=' Error in READ_ANALYSIS: Subdirective not yet'
     &           //' implemented. Abort.'
            GOTO 500
         END IF
      ELSE
         GOTO 500
      END IF
      GOTO 100

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

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

600   CONTINUE
      iret=0
      RETURN

500   CONTINUE
      CALL xerror(errmsg,80,1,2)
      iret=1

      RETURN
c==============================================================================
c     Errors were found
c==============================================================================

 20   CONTINUE
      errmsg='internal reading error: wrong format?? TAB character??'
      CALL xerror(errmsg,80,1,2)
      STOP
      END
