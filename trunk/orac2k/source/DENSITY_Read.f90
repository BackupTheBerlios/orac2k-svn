  SUBROUTINE Read_it(knlist,kprint,nsevere,nword,strngs,iret,errmsg&
       &,read_err)

!!$======================== DECLARATIONS ================================*

    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

    INTEGER :: knlist,kprint,nword,nsevere,iret,read_err
    CHARACTER(80) :: strngs(:),errmsg

!!$------------------------- LOCAL VARIABLES ----------------------------*

    CHARACTER(80) :: line
    REAL(8) :: dummy
    LOGICAL :: exist

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

!!$
!!$==== Command DENSITY  ================================================
!!$

    line(79:80)='  '
    read_err=0
    Density_Calc=.TRUE.
    DO
       READ(knlist,'(a78)',END=600) line(1:78)
       CALL wrenc(kprint,line)
       IF(line(1:1) .EQ. '#') CYCLE
       CALL parser(line,strngs,nword)
!!$
!!$==== Subcommand density grid =========================================
!!$
       SELECT CASE(strngs(1))
       CASE DEFAULT
          errmsg=err_unr(2)//strngs(2)//err_end(1:14)//err_end(16:20)
          CALL xerror(errmsg,80,1,30)
          nsevere = nsevere + 1

       CASE('grid')
          IF(nword .LT. 4) THEN
             errmsg=err_args(1)//'3'
             CALL xerror(errmsg,80,1,30)
             nsevere=nsevere+1
          ELSE
             CALL Read_String(strngs(2),ncx)
             CALL Read_String(strngs(3),ncy)
             CALL Read_String(strngs(4),ncz)
          END IF
!!$
!!$==== subcommand target =============================================
!!$
       CASE('target')
          IF(nword .LT. 2) THEN
             errmsg=err_args(1)//'1'
             CALL xerror(errmsg,80,1,30)
             nsevere=nsevere+1
          ELSE
             Target_Res=strngs(2)
          END IF
!!$
!!$==== subcommand write ==============================================
!!$
       CASE('write')
          IF(nword .LT. 2) THEN
             errmsg=err_args(1)//'1'
             CALL xerror(errmsg,80,1,30)
             nsevere=nsevere+1
          ELSE
             CALL Read_String(strngs(2),dummy)
             n_write=IDINT(dummy)
          END IF
!!$
!!$==== subcommand out_format =========================================
!!$
       CASE('out_format')
          IF(nword .LT. 2) THEN
             errmsg=err_args(1)//'1'
             CALL xerror(errmsg,80,1,30)
             nsevere=nsevere+1
          ELSE
             file_format=strngs(2)(1:8)
          END IF
!!$
!!$==== subcommand Density_avg =========================================
!!$
       CASE('density_avg')
          Density_Avg=.TRUE.
!!$
!!$--------------------------------------------------------------------
!!$
       CASE(' ')
          CYCLE

       CASE('END')
          EXIT
       END SELECT
    END DO

    SELECT CASE(file_format)
    CASE DEFAULT
       filename='DENSITY_FILE.cube'
       INQUIRE(FILE=filename,EXIST=exist)
       IF(exist) THEN
          CALL openf(kdensity,filename,'FORMATTED','OLD',0)
       ELSE
          CALL openf(kdensity,filename,'FORMATTED','NEW',0)
       END IF
    CASE('xplor')
       filename='DENSITY_FILE.xplor'
       INQUIRE(FILE=filename,EXIST=exist)
       IF(exist) THEN
          CALL openf(kdensity,filename,'FORMATTED','OLD',0)
       ELSE
          CALL openf(kdensity,filename,'FORMATTED','NEW',0)
       END IF
    END SELECT
    RETURN

600 read_err=1
    RETURN

20  CONTINUE
    iret=1
    errmsg='internal reading error: wrong format?? TAB character??'
    CALL xerror(errmsg,80,1,2)
    
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

  END SUBROUTINE Read_it
