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

    ALLOCATE(Res_u(n_Res))
    line(79:80)='  '
    read_err=0
    RMS_Subtract=.TRUE.
    DO
       READ(knlist,'(a78)',END=600) line(1:78)
       CALL wrenc(kprint,line)
       IF(line(1:1) .EQ. '#') CYCLE
       CALL parser(line,strngs,nword)
!!$
!!$==== Subcommand RMS_SUBTRACT =========================================
!!$
       SELECT CASE(strngs(1))
       CASE DEFAULT
          errmsg=err_unr(2)//strngs(2)//err_end(1:14)//err_end(16:20)
          CALL xerror(errmsg,80,1,30)
          nsevere = nsevere + 1

!!$
!!$==== subcommand unit ===============================================
!!$
       CASE('unit')
          units=units+1
          IF(units > n_Res_u) THEN
             errmsg='Number of Units exceed maximum.'
             CALL xerror(errmsg,80,1,30)
             nsevere=nsevere+1
          END IF
          IF(nword .LT. 2) THEN
             errmsg=err_args(1)//'1'
             CALL xerror(errmsg,80,1,30)
             nsevere=nsevere+1
          ELSE
             Res_u(units)%target=strngs(2)
             Res_u(units)%n0=0
             ALLOCATE(Res_u(units)%beta(n_Res))
             DO 
                READ(knlist,'(a78)',END=600) line(1:78)
                CALL wrenc(kprint,line)
                IF(line(1:1) .EQ. '#') CYCLE
                CALL parser(line,strngs,nword)
                IF(Res_u(units)%n0+nword > n_Res) THEN
                   errmsg='Maximum of Atoms in units reached.'
                   CALL xerror(errmsg,80,1,30)
                   nsevere=nsevere+1
                END IF

                IF(strngs(nword) == 'end') THEN
                   Res_u(units)%beta(Res_u(units)%n0+1:Res_u(units)%n0&
                        &+nword-1)=strngs(1:nword-1)(1:5)
                   Res_u(units)%n0=Res_u(units)%n0+nword-1
                   EXIT
                END IF
                Res_u(units)%beta(Res_u(units)%n0+1:Res_u(units)%n0&
                     &+nword)=strngs(1:nword)(1:5)
                Res_u(units)%n0=Res_u(units)%n0+nword
             END DO
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
!!$==== subcommand charges ============================================
!!$
       CASE('rms_type')
          IF(nword .LT. 2) THEN
             errmsg=err_args(1)//'1'
             CALL xerror(errmsg,80,1,30)
             nsevere=nsevere+1
          ELSE
             SELECT CASE(TRIM(strngs(2)))
             CASE('dipole')
                CONTINUE
             CASE('mechanic')
                CONTINUE
             CASE('cm')
                CONTINUE
             CASE DEFAULT
                errmsg='Illegal rms_type. Only dipole, mechanic and'//&
                     &' cm are legal types.'
                CALL xerror(errmsg,80,1,30)
                nsevere=nsevere+1
             END SELECT
             RMS_TYPE=TRIM(strngs(2))
          END IF

!!$
!!$==== subcommand no_inversion  ======================================
!!$
       CASE('no_inversion')
          RMS_Invert=.FALSE.
!!$
!!$==== subcommand test  ==============================================
!!$
       CASE('test')
          RMS_test=.TRUE.
!!$
!!$--------------------------------------------------------------------
!!$
       CASE(' ')
          CYCLE

       CASE('END')
          EXIT
       END SELECT
    END DO

    filename='RMS_FILE.rms'
    INQUIRE(FILE=filename,EXIST=exist)
    IF(exist) THEN
       CALL openf(krms_sub,filename,'FORMATTED','OLD',0)
    ELSE
       CALL openf(krms_sub,filename,'FORMATTED','NEW',0)
    END IF
    filename_dip='RMS_FILE.dip'
    INQUIRE(FILE=filename_dip,EXIST=exist)
    IF(exist) THEN
       CALL openf(krms_dip,filename_dip,'FORMATTED','OLD',0)
    ELSE
       CALL openf(krms_dip,filename_dip,'FORMATTED','NEW',0)
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
