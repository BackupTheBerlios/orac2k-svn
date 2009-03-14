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
    LOGICAL :: exist,print_ok

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

!!$
!!$==== Command DENSITY  ================================================
!!$

    line(79:80)='  '
    read_err=0
    voronoi=.TRUE.
    print_OK=.FALSE.
    DO
       READ(knlist,'(a78)',END=600) line(1:78)
       CALL wrenc(kprint,line)
       IF(line(1:1) .EQ. '#') CYCLE
       CALL parser(line,strngs,nword)
!!$
!!$==== Subcommand VORONOI =========================================
!!$
       SELECT CASE(TRIM(strngs(1)))
       CASE DEFAULT
          errmsg=err_unr(2)//strngs(2)//err_end(1:14)//err_end(16:20)
          CALL xerror(errmsg,80,1,30)
          nsevere = nsevere + 1

!!$
!!$==== subcommand unit ===============================================
!!$
       CASE('cutoff')
          CALL Read_String(strngs(2),dummy)
          cutoff=dummy

       CASE('heavy_atoms')
          heavy=.TRUE.

       CASE('rewind')
          rewind_vor=.TRUE.

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


       CASE('compute')
          IF(TRIM(strngs(2)) == 'accessibility') THEN
             access=.TRUE.

          ELSE IF(TRIM(strngs(2)) == 'access') THEN
             access=.TRUE.
             
          ELSE IF(TRIM(strngs(2)) == 'volume') THEN
             volume=.TRUE.
             
          ELSE IF(TRIM(strngs(2)) == 'neighbors') THEN
             neighbor=.TRUE.

          ELSE IF(TRIM(strngs(2)) == 'compress') THEN
             compress=.TRUE.
             neighbor=.TRUE.
             volume=.TRUE.
             noprint=.TRUE.
             IF(nword == 3) THEN
                IF(TRIM(strngs(3)) == 'kba') kba=.TRUE.
             END IF

          ELSE IF(TRIM(strngs(2)) == 'fluctuations' .OR.&
               & TRIM(strngs(2)) == 'fluct') THEN
             fluct=.TRUE.
             IF(nword == 3) THEN
                CALL Read_String(strngs(3),dummy)
                nfluct=IDINT(dummy)
             END IF             
          ELSE
             errmsg=err_unr(2)//strngs(2)
             CALL xerror(errmsg,80,1,30)
             nsevere = nsevere + 1
          END IF

       CASE('dynamics')
          dynamics=.TRUE.

       CASE('every')
          CALL Read_String(strngs(2),dummy)
          nvoronoi=IDINT(dummy)

       CASE('orient_bin')
          CALL Read_String(strngs(2),dummy)
          h_bin1=dummy
          CALL Read_String(strngs(3),dummy)
          h_bin2=dummy

       CASE('orien_bin')
          CALL Read_String(strngs(2),dummy)
          h_bin1=dummy
          CALL Read_String(strngs(3),dummy)
          h_bin2=dummy

       CASE('bin_size' )
          CALL Read_String(strngs(2),dummy)
          bin=dummy
          
       CASE('print_press' )
          CALL Read_String(strngs(2),dummy)
          nprint_press=IDINT(dummy)
          
       CASE('print' )
          print_ok=.TRUE.
          CALL Read_String(strngs(2),dummy)
          nvoronoi=IDINT(dummy)
          IF(strngs(3) .EQ. 'OPEN') THEN
             CALL uscrpl(strngs(4),80)
             filename=strngs(4)
             INQUIRE(FILE=filename,EXIST=exist)
             IF(exist) THEN
                CALL openf(kvoronoi,strngs(4),'FORMATTED','OLD',0)
             ELSE
                CALL openf(kvoronoi,strngs(4),'FORMATTED','NEW',0)
             END IF
          ELSE
             errmsg='OPEN keyword not found'
             CALL xerror(errmsg,80,1,30)
             nsevere = nsevere + 1
          END IF
          
!!$
!!$--------------------------------------------------------------------
!!$
       CASE('no_print' )
          noprint=.TRUE.
       CASE(' ')
          CYCLE

       CASE('END')
          EXIT
       END SELECT
    END DO

    IF(Dynamics) THEN
       filename='VOR_DYNAMICS.dat'
       INQUIRE(FILE=filename,EXIST=exist)
       IF(exist) THEN
          CALL openf(kdynamics,filename,'UNFORMATTED','OLD',0)
       ELSE
          CALL openf(kdynamics,filename,'UNFORMATTED','NEW',0)
       END IF
    END IF

    file_cosa='COSA.hist'
    INQUIRE(FILE=file_cosa,EXIST=exist)
    IF(exist) THEN
       CALL openf(kcosa,file_cosa,'UNFORMATTED','OLD',0)
    ELSE
       CALL openf(kcosa,file_cosa,'UNFORMATTED','NEW',0)
    END IF
    CLOSE(kcosa)

    IF(compress) THEN
       filename='VOR_HISTO.dat'
       INQUIRE(FILE=filename,EXIST=exist)
       IF(exist) THEN
          CALL openf(kvoronoi,filename,'FORMATTED','OLD',0)
       ELSE
          CALL openf(kvoronoi,filename,'FORMATTED','NEW',0)
       END IF
       IF(ncx == -1 .OR. ncy == -1 .OR. ncz == -1) THEN
          errmsg='Need command grid to compute density'
          CALL xerror(errmsg,80,1,30)
          nsevere = nsevere + 1
       END IF
    END IF

    IF(.NOT. print_ok .AND. (.NOT. noprint)) THEN
       filename='VORONOI.dat'
       INQUIRE(FILE=filename,EXIST=exist)
       IF(exist) THEN
          CALL openf(kvoronoi,filename,'FORMATTED','OLD',0)
       ELSE
          CALL openf(kvoronoi,filename,'FORMATTED','NEW',0)
       END IF
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
