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
       SELECT CASE(strngs(1))
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

       CASE('compute')
          IF(strngs(2) .EQ. 'accessibility') THEN
             access=.TRUE.
             
          ELSE IF(strngs(2) .EQ. 'volume') THEN
             volume=.TRUE.
             
          ELSE IF(strngs(2) .EQ. 'neighbors') THEN
             neighbor=.TRUE.
             
          ELSE IF(strngs(2) .EQ. 'fluctuations') THEN
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

       CASE(' ')
          CYCLE

       CASE('END')
          EXIT
       END SELECT
    END DO

    IF(.NOT. print_ok) THEN
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
