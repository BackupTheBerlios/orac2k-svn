      SUBROUTINE read_run(fprop,frject,fmaxrun,fmaxstp,fprint,err_args
     &     ,err_unr,err_end)

************************************************************************
*   Time-stamp: <01/03/11 11:09:27 marchi>                             *
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
      
      CHARACTER*37 err_args(2)
      CHARACTER*20 err_end 
      CHARACTER*27 err_unr(4)
      REAL*8 fprop,frject,fmaxrun,fmaxstp,fprint

*----------------------- VARIABLES IN COMMON --------------------------*
      
      INCLUDE 'parst.h'
      INCLUDE 'cpropar.h'
      INCLUDE 'unit.h'

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER iret
      CHARACTER*80 errmsg
      INTEGER nword,nsevere,nwarning,j
      CHARACTER*80 line,strngs(40)
      CHARACTER*8 fmt
      CHARACTER*1 sep(2),comm(2)
      DATA sep/' ',','/comm/'(',')'/

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      line(79:80)='  '
      nsevere = 0
      nwarning = 0

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

c==== Command  CONTROL================================================

      IF(strngs(1).EQ. 'CONTROL' ) THEN
         IF(nword.eq.2) THEN 
            CALL fndfmt(1,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) nflag(1)
         ELSE
            errmsg=err_args(1) //'1'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         
c==== Command  OPTION==================================UNSUPPORTED====

      ELSE IF(strngs(1).EQ. 'OPTION' ) THEN
         errmsg = err_unr(4)
         nwarning = nwarning + 1
         call xerror(errmsg,80,1,11)
         CALL fndfmt(1,strngs(2),fmt)
         READ(strngs(2),fmt,err=20) nflag(2)
         
c==== Command  PROPERTY===============================================

      ELSE IF(strngs(1).EQ. 'PROPERTY' ) THEN
         IF(nword.eq.2) THEN 
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) fprop
         ELSE
            errmsg=err_args(1) //'1'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF

c==== Command  REJECT=================================================

      ELSE IF(strngs(1).EQ. 'REJECT' ) THEN
         IF(nword.eq.2) THEN 
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) frject
         ELSE
            errmsg=err_args(1) //'1'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         
c==== Command  MAXRUN=================================================

      ELSE IF(strngs(1).EQ. 'MAXRUN' ) THEN
         IF(nword.eq.2) THEN 
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) fmaxrun
         ELSE
            errmsg=err_args(1) //'1'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         
c==== Command  TIME===================================================

      ELSE IF(strngs(1).EQ. 'TIME' ) THEN
         IF(nword.eq.2) THEN 
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) fmaxstp
         ELSE
            errmsg=err_args(1) //'1'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         
c==== Command  TIME_LIMIT ============================================

      ELSE IF(strngs(1).EQ. 'TIME_LIMIT' ) THEN
         IF(nword.eq.2) THEN 
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) TimeLimit
            TimeLimit=TimeLimit*3600.0D0
         ELSE
            errmsg=err_args(1) //'1'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         
c==== Command  ADD_TIME ==============================================

      ELSE IF(strngs(1).EQ. 'ADD_TIME' ) THEN
         AddTime=.TRUE.
         IF(nword.eq.2) THEN 
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) fmaxstp
         ELSE
            errmsg=err_args(1) //'1'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         
c==== Command  DEBUG==================================================

      ELSE IF(strngs(1).EQ. 'DEBUG' ) THEN
         DEBUG=.TRUE. 
         if(nword.eq.2) THEN
            if (strngs(2).EQ.'ALL'.OR.strngs(2).EQ.'RESIDUE_SEQUENCE')
     &           THEN
               debug_rs=.TRUE.
            else if (strngs(2).EQ.'ALL'.OR.strngs(2).EQ
     &              .'CONNECTION_TABLE')THEN
               debug_ct=.TRUE.
            else if (strngs(2).EQ.'ALL'.OR.strngs(2).EQ.'BOND_TABLE')
     &              THEN
               debug_st=.TRUE.
            else if (strngs(2).EQ.'ALL'.OR.strngs(2).EQ.'BEND_TABLE')
     &              THEN
               debug_bt=.TRUE.
            else if (strngs(2).EQ.'ALL'.OR.strngs(2).EQ.'PTORS_TABLE')
     &              THEN
               debug_pt=.TRUE.
            else if (strngs(2).EQ.'ALL'.OR.strngs(2).EQ.'ITORS_TABLE')
     &              THEN
               debug_it=.TRUE.
               
            ELSE
               errmsg='Unrecognized keyword '//strngs(2)
               CALL xerror(errmsg,80,1,30)
               nsevere = nsevere + 1
            END IF
         END IF

c==== Command  PRINT==================================================

      ELSE IF(strngs(1).EQ. 'PRINT' ) THEN
         IF(nword.eq.2) THEN 
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20)  fprint
         ELSE
            errmsg=err_args(1) //'1'
            CALL xerror(errmsg,80,1,30)
            nsevere = nsevere + 1
         END IF
         
c==== Blank Line======================================================

      ELSE IF(strngs(1).EQ. ' ') THEN
         CONTINUE

      ELSE IF(strngs(1)(1:1).EQ. '&'.AND.strngs(1).NE. '&END') THEN
         errmsg=err_unr(1) //strngs(1)(1:8)// err_end
         CALL xerror(errmsg,80,1,30)
         nsevere = nsevere + 1
         GO TO 600

c==== &END ===========================================================

      ELSE IF(strngs(1).EQ. '&END') THEN
         GOTO 600

      ELSE
         errmsg=err_unr(1) //strngs(1)(1:8)// err_end
         CALL xerror(errmsg,80,1,30)
         nsevere = nsevere + 1
      END IF

      GOTO 100

600   CONTINUE

c=======================================================================
c     Environment parser ends here 
c=======================================================================

c--   syntax errors: abort without verifying input 
      if(nsevere.gt.0.and.nsevere.lt.99) then 
         j=0
         call int_str(nsevere,fmt,j)
         errmsg=fmt(1:j) //' ERRORS WHILE EXECUTING READ_RUN'
         CALL xerror(errmsg,80,1,2)
         STOP
      ELSE IF(nsevere.gt.99) THEN 
         errmsg= 'MORE THAN 99 ERRORS WHILE EXECUTING READ_RUN'
         call xerror(errmsg,80,1,2)
         STOP
      END IF
      if(nwarning.gt.0.and.nwarning.lt.99) then 
         j=0
         call int_str(nwarning,fmt,j)
         errmsg= fmt(1:j)//'WARNINGS WHILE EXECUTING READ_RUN'
         CALL xerror(errmsg,80,1,1)
      ELSE IF(nwarning.gt.99) THEN 
         errmsg= 'MORE THAN 99 ERRORS WHILE EXECUTING READ_RUN'
         call xerror(errmsg,80,1,1)
      ENDIF    

      RETURN 

c==============================================================================
c     Errors were found
c==============================================================================

 20   CONTINUE
      errmsg='internal reading error: wrong format?? TAB character??'
      CALL xerror(errmsg,80,1,2)
      STOP

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      END
