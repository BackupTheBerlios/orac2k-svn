  SUBROUTINE Read_it(knlist,kprint,nsevere,nword,strngs,iret,errmsg&
       &,read_err,mm1)

!!$======================== DECLARATIONS ================================*

    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

    INTEGER :: knlist,kprint,nword,nsevere,iret,read_err,mm1
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

       CASE('cell')
          IF(nword .LT. 4) THEN
             errmsg=err_args(1)//'3'
             CALL xerror(errmsg,80,1,30)
             nsevere=nsevere+1
          ELSE
             CALL Read_String(strngs(2),a)
             CALL Read_String(strngs(3),b)
             CALL Read_String(strngs(4),c)
          END IF

       CASE('histo')
          Dens_histo=.True.
          IF(nword .Eq. 3) THEN
             CALL Read_String(strngs(2),Dens_Histo_rmax)
             CALL Read_String(strngs(3),Dens_Histo_size)
          END IF

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

       CASE('molecule')
          IF(natoms_Tot == 0) THEN
             ALLOCATE(atoms(mm1))
          END IF
          CALL parse_numbers(err_unr,strngs,nword,atoms(ntot),nats,nsevere)
          atoms(ntot)=Nats
          ntot=ntot+Nats+1
          natoms_Tot=natoms_Tot+1
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
       filename_den='DENSITY_FILE.cube'
       INQUIRE(FILE=filename_den,EXIST=exist)
       IF(exist) THEN
          CALL openf(kdensity,filename_den,'FORMATTED','OLD',0)
       ELSE
          CALL openf(kdensity,filename_den,'FORMATTED','NEW',0)
       END IF


       filename_den2='DENS2_FILE.cube'
       INQUIRE(FILE=filename_den2,EXIST=exist)
       IF(exist) THEN
          CALL openf(kdens2,filename_den2,'FORMATTED','OLD',0)
       ELSE
          CALL openf(kdens2,filename_den2,'FORMATTED','NEW',0)
       END IF
    CASE('xplor')
       filename_den='DENSITY_FILE.xplor'
       INQUIRE(FILE=filename_den,EXIST=exist)
       IF(exist) THEN
          CALL openf(kdensity,filename_den,'FORMATTED','OLD',0)
       ELSE
          CALL openf(kdensity,filename_den,'FORMATTED','NEW',0)
       END IF

       filename_den2='DENS2_FILE.xplor'
       INQUIRE(FILE=filename_den2,EXIST=exist)
       IF(exist) THEN
          CALL openf(kdens2,filename_den2,'FORMATTED','OLD',0)
       ELSE
          CALL openf(kdens2,filename_den2,'FORMATTED','NEW',0)
       END IF
    END SELECT
    IF(natoms_Tot /=0) THEN
       filename_pdb='DENSITY_FILE.pdb'
       INQUIRE(FILE=filename_pdb,EXIST=exist)
       IF(exist) THEN
          CALL openf(kpdb,filename_pdb,'FORMATTED','OLD',0)
       ELSE
          CALL openf(kpdb,filename_pdb,'FORMATTED','NEW',0)
       END IF
    END IF
    If(Dens_Histo) Then
       filename_histo='DENS_HISTO.dat'
       filename_histo_cos='DENS_HISTO_COS.dat'
       Inquire(FILE=filename_histo,EXIST=exist)
       If(exist) Then
          Call openf(kdenhisto,filename_histo,'FORMATTED','OLD',0)
       Else
          Call openf(kdenhisto,filename_histo,'FORMATTED','NEW',0)
       End If
       Inquire(FILE=filename_histo_cos,EXIST=exist)
       If(exist) Then
          Call openf(kdenhisto_cos,filename_histo_cos,'FORMATTED','OLD',0)
       Else
          Call openf(kdenhisto_cos,filename_histo_cos,'FORMATTED','NEW',0)
       End If
    End If
    Return

600 read_err=1
    RETURN

20  CONTINUE
    iret=1
    errmsg='internal reading error: wrong format?? TAB character??'
    CALL xerror(errmsg,80,1,2)
    
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

  END SUBROUTINE Read_it
