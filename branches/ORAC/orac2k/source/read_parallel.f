      SUBROUTINE read_parallel(err_args,err_unr,err_end)

************************************************************************
*   Time-stamp: <99/03/04 12:15:19 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Mar  4 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      CHARACTER*37 err_args(2)
      CHARACTER*20 err_end 
      CHARACTER*27 err_unr(4)

*----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'parst.h'
      INCLUDE 'cpropar.h'
      INCLUDE 'unit.h'

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER nword,nsevere,nwarning,j,iret
      CHARACTER*80 errmsg
      CHARACTER*80 line,strngs(40)
      CHARACTER*8 fmt
      CHARACTER*1 sep(2),comm(2)
      LOGICAL  exist
#if defined T3E | defined _CRAY_
      INTEGER*8 one
      DATA one/1/
#endif
      DATA sep/' ',','/comm/'(',')'/

*----------------------- EXECUTABLE STATEMENTS ------------------------*


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

c==== Command P_update===================================================

      IF(strngs(1).EQ. 'P_update' ) THEN
         P_dyn_update=.TRUE.
         IF(nword .NE. 1) THEN
            IF(nword .EQ. 2) THEN
               P_dyn_update_shell=strngs(2)(1:1)
            ELSE
               nsevere = nsevere + 1
               errmsg=err_args(1) / /'2'
               CALL xerror(errmsg,80,1,30)
            END IF
         END IF

#if defined T3E | defined _CRAY_
c==== Command P_stream===================================================

      ELSE IF(strngs(1).EQ. 'P_stream_r8' ) THEN
         CALL set_d_stream(one)

      ELSE IF(strngs(1).EQ. 'P_stream_i' ) THEN
         CALL set_i_stream(one)
#endif

c==== Command P_nei_fact===================================================

      ELSE IF(strngs(1).EQ. 'P_nei_fact' ) THEN
         IF(nword .EQ. 2) THEN
            CALL fndfmt(2,strngs(2),fmt)
            READ(strngs(2),fmt,err=20) P_nei_fact
         ELSE
            nsevere = nsevere + 1
            errmsg=err_args(1) / /'2'
            CALL xerror(errmsg,80,1,30)
         END IF

c==== Command &END ====================================================

      ELSE IF(strngs(1).EQ. '&END') THEN
         GOTO 600
         
      ELSE
         errmsg= err_unr(1) / / strngs(1)(1:8) / / err_end
         CALL xerror(errmsg,80,1,30)
         nsevere = nsevere + 1
      END IF

      GO TO 100

600   CONTINUE



*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN

 20   CONTINUE
      iret=1
      errmsg='Internal reading error: wrong format?? TAB character.'
      CALL xerror(errmsg,80,1,2)

      STOP
      END
