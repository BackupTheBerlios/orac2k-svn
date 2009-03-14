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
!!$==== Command FOURIER  ================================================
!!$

    
    line(79:80)='  '
    read_err=0
    Electrostatics=.TRUE.
    DO
       READ(knlist,'(a78)',END=600) line(1:78)
       CALL wrenc(kprint,line)
       IF(line(1:1) .EQ. '#') CYCLE
       CALL parser(line,strngs,nword)
!!$
!!$==== Subcommand PME =================================================
!!$
       SELECT CASE(strngs(1))
       CASE DEFAULT
          errmsg=err_unr(2)//strngs(2)//err_end(1:14)//err_end(16:20)
          CALL xerror(errmsg,80,1,30)
          nsevere = nsevere + 1

!!$
!!$==== subcommand unit ===============================================
!!$
       CASE('grid')
          IF(nword .LT. 4) THEN
             errmsg=err_args(1)//'3'
             CALL xerror(errmsg,80,1,30)
             nsevere=nsevere+1
          ELSE
             CALL Read_String(strngs(2),nfft1)
             CALL Read_String(strngs(3),nfft2)
             CALL Read_String(strngs(4),nfft3)
          END IF

!!$
!!$==== subcommand order ==============================================
!!$
       CASE('order')
          IF(nword .LT. 2) THEN
             errmsg=err_args(1)//'1'
             CALL xerror(errmsg,80,1,30)
             nsevere=nsevere+1
          ELSE
             CALL Read_String(strngs(2),order)
          END IF

!!$
!!$==== subcommand alpha ==============================================
!!$
       CASE('alpha')
          IF(nword .LT. 2) THEN
             errmsg=err_args(1)//'1'
             CALL xerror(errmsg,80,1,30)
             nsevere=nsevere+1
          ELSE
             CALL Read_String(strngs(2),alpha)
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
!!$--------------------------------------------------------------------
!!$
       CASE(' ')
          CYCLE

       CASE('END')
          EXIT
       END SELECT
    END DO

    filename='FOURIER_POTENTIAL.xplor'
    INQUIRE(FILE=filename,EXIST=exist)
    IF(exist) THEN
       CALL openf(kfourier,filename,'FORMATTED','OLD',0)
    ELSE
       CALL openf(kfourier,filename,'FORMATTED','NEW',0)
    END IF
    RETURN

600 read_err=1
    RETURN

20  CONTINUE
    iret=1
    errmsg='internal reading error: wrong format?? TAB character??'
    CALL xerror(errmsg,80,1,2)
    
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

  END SUBROUTINE Read_it
