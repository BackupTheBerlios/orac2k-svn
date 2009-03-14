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
    WSC__=.TRUE.
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

       CASE('molecule')
          IF(natoms_Tot == 0) THEN
             ALLOCATE(atoms(mm1))
          END IF
          CALL parse_numbers(err_unr,strngs,nword,atoms(ntot),nats,nsevere)
          atoms(ntot)=Nats
          ntot=ntot+Nats+1
          natoms_Tot=natoms_Tot+1

       CASE('simple')
          simple=.TRUE. 
!!$
!!$--------------------------------------------------------------------
!!$
       CASE(' ')
          CYCLE

       CASE('END')
          EXIT
       END SELECT
    END DO
    RETURN

600 read_err=1
    RETURN

20  CONTINUE
    iret=1
    errmsg='internal reading error: wrong format?? TAB character??'
    CALL xerror(errmsg,80,1,2)
    
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

  END SUBROUTINE Read_it
