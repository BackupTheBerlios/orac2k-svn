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

    line(79:80)='  '
    read_err=0
    PDBs=.TRUE.
    DO 
       READ(knlist,'(a78)',END=600) line(1:78)                    
       CALL wrenc(kprint,line)
       IF(line(1:1) .EQ. '#') CYCLE
       CALL parser(line,strngs,nword)
!!$
!!$==== subcommand range ==============================================  
!!$
       SELECT CASE(strngs(1))
       CASE DEFAULT                                                
          errmsg=err_unr(2)//strngs(2)//err_end(1:14)&
               &//err_end(16:20)
          CALL xerror(errmsg,80,1,30)                             
          nsevere = nsevere + 1
       CASE('range')
          IF(nword < 3) THEN
             errmsg=err_args(1)//'2'
             CALL xerror(errmsg,80,1,30)                          
             nsevere=nsevere+1                               
          ELSE
             CALL Read_String(strngs(2),Residue_Range(1))
             CALL Read_String(strngs(3),Residue_Range(2))
          END IF
!!$
!!$==== subcommands fits  =============================================  
!!$
       CASE('fit_to')
          IF(nword < 2) THEN
             errmsg=err_args(1)//'1'                              
             CALL xerror(errmsg,80,1,30)                          
             nsevere=nsevere+1                                    
          ELSE
             Template%fit=.TRUE.
             Template%type=strngs(2)(1:3)
          END IF          
!!$
!!$==== subcommand write ==============================================  
!!$
       CASE('write')
          IF(nword < 2) THEN
             errmsg=err_args(1)//'1'                              
             CALL xerror(errmsg,80,1,30)                          
             nsevere=nsevere+1                                    
          ELSE IF(nword == 2) THEN
             CALL Read_String(strngs(2),dummy)
             n_write=IDINT(dummy)
             strngs(3)='PDB_FILE.pdb'
             filename=strngs(3)
             INQUIRE(FILE=strngs(3),EXIST=exist)
             IF(exist) THEN
                CALL openf(kpdb,strngs(3),'FORMATTED','OLD',0)
             ELSE
                CALL openf(kpdb,strngs(3),'FORMATTED','NEW',0)
             END IF

          ELSE IF(nword == 3) THEN
             CALL Read_String(strngs(2),dummy)
             n_write=IDINT(dummy)
             filename=strngs(3)
             INQUIRE(FILE=strngs(3),EXIST=exist)
             IF(exist) THEN
                CALL openf(kpdb,strngs(3),'FORMATTED','OLD',0)
             ELSE
                CALL openf(kpdb,strngs(3),'FORMATTED','NEW',0)
             END IF
          END IF
!!$
!!$==== subcommand compute ============================================  
!!$

       CASE('compute')
          IF(nword < 2) THEN
             errmsg=err_args(1)//'1'                              
             CALL xerror(errmsg,80,1,30)                          
             nsevere=nsevere+1                                    
          ELSE
             CALL Read_String(strngs(2),dummy)
             n_compute=IDINT(dummy)
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

    RETURN

600 read_err=1
    RETURN

20  CONTINUE
    iret=1
    errmsg='internal reading error: wrong format?? TAB character??'
    CALL xerror(errmsg,80,1,2)

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

  END SUBROUTINE Read_it
