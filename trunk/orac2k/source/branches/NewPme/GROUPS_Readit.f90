SUBROUTINE Read_it(knlist,kprint,nsevere,nword,strngs,iret,errmsg&
     &,read_err)

!!$======================== DECLARATIONS ================================*

  IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

  INTEGER :: knlist,kprint,nword,nsevere,iret,read_err
    CHARACTER(80) :: errmsg
    CHARACTER(80), TARGET :: strngs(:)

!!$------------------------- LOCAL VARIABLES ----------------------------*

    CHARACTER(80) :: line,aux
    CHARACTER(80), DIMENSION(:), POINTER :: st_remains
    INTEGER :: count=0

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

!!$
!!$==== Command GROUPS  ================================================
!!$

    count=0
    ALLOCATE(molecs(n_molecs))
    molecs(:)%n=0
    molecs(:)%mode=' '
    molecs(:)%type=' '

    line(79:80)='  '
    read_err=0
    groups=.TRUE.
    DO
       READ(knlist,'(a78)',END=600) line(1:78)
       CALL wrenc(kprint,line)
       IF(line(1:1) .EQ. '#') CYCLE
       CALL parser(line,strngs,nword)
       
       SELECT CASE(TRIM(strngs(1)))
       CASE DEFAULT
          errmsg=err_unr(2)//strngs(2)//err_end(1:14)//err_end(16:20)
          CALL xerror(errmsg,80,1,30)
          nsevere = nsevere + 1
          
!!$
!!$==== Subcommand atoms  ===========================================
!!$
       CASE('atoms')
          IF(nword .LT. 3) THEN
             errmsg=err_args(1)//'3'
             CALL xerror(errmsg,80,1,30)
             nsevere=nsevere+1
          ELSE
             count=count+1
             IF(count > n_molecs)  THEN
                WRITE(aux,'(i4)') n_molecs
                errmsg='Total number of groups cannot exceed '//TRIM(aux)//&
                     &' Abort.'
                CALL xerror(errmsg,80,1,30)
                nsevere=nsevere+1
                EXIT
             END IF
             molecs(count)%mode='atoms'
             molecs(count)%type=TRIM(strngs(2))
             st_remains=>strngs(3:nword)
             CALL Parse_Numbers(st_remains,molecs(count)%index)
             molecs(count)%n=SIZE(molecs(count)%index)
          END IF
!!$
!!$==== Subcommand residue  =========================================
!!$
       CASE('residue')
          IF(nword .LT. 3) THEN
             errmsg=err_args(1)//'3'
             CALL xerror(errmsg,80,1,30)
             nsevere=nsevere+1
          ELSE
             count=count+1
             IF(count > n_molecs)  THEN
                WRITE(aux,'(i4)') n_molecs
                errmsg='Total number of groups cannot exceed '//TRIM(aux)//&
                     &' Abort.'
                CALL xerror(errmsg,80,1,30)
                nsevere=nsevere+1
                EXIT
             END IF
             molecs(count)%mode='residue'
             molecs(count)%type=TRIM(strngs(2))
             st_remains=>strngs(3:nword)
             CALL Parse_Numbers(st_remains,molecs(count)%index)
             molecs(count)%n=SIZE(molecs(count)%index)
          END IF
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
