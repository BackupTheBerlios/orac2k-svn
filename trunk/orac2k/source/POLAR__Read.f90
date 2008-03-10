  SUBROUTINE Read_it(knlist,kprint,nsevere,nword,strngs,iret,errmsg&
       &,read_err)

!!$======================== DECLARATIONS ================================*

    IMPLICIT none
    INCLUDE 'parst.h'
    INCLUDE 'cpropar.h'

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

    Polar__=.TRUE.                                              
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


       CASE('ewald')
          IF(nword .EQ. 6) THEN
             CALL Read_String(strngs(2),alphal)
             CALL Read_String(strngs(3),nfft1)
             CALL Read_String(strngs(4),nfft2)
             CALL Read_String(strngs(5),nfft3)
             CALL Read_String(strngs(6),pme_order)
          ELSE                                                 
             nsevere = nsevere + 1                             
             errmsg=err_args(1) // '5 after keyword'     
             CALL xerror(errmsg,80,1,30)                       
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

!!$c==== Command ERFC_SPLINE==============================================
                                                                       
       CASE('erfc_spline')
          erfc_spline=.TRUE.                                      
          IF(nword .NE. 1) THEN                                   
             CALL Read_String(strngs(2),erfc_bin)
          END IF
                                                                       
!!$c==== subcommand CUTOFF =============================================  
          
       CASE('cutoff')
          IF(nword.ne.1) THEN
             CALL Read_String(strngs(2),rspoff)
             rneim=0.0D0
             rneil=0.0D0
             rneih=1.5D0
             rtolm=0.0D0
             rtoll=0.0D0
             rtolh=0.5D0
             rcuth=rspoff-rtolm
             rcutl=0.0D0
             rcutm=0.0D0
          END IF
          
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

       CASE('grid')
          IF(nword .LT. 4) THEN
             errmsg=err_args(1)//'3'
             CALL xerror(errmsg,80,1,30)
             nsevere=nsevere+1
          ELSE
             CALL Read_String(strngs(2),nccx)
             CALL Read_String(strngs(3),nccy)
             CALL Read_String(strngs(4),nccz)
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

       CASE('P_of_r')
          do_Pr=.TRUE.
          IF(nword /= 1) THEN
             CALL Read_String(strngs(2),dx_Pr)
             CALL Read_String(strngs(3),max_Pr)
          END IF

       CASE('molecule')
          IF(natoms_Tot == 0) THEN
             ALLOCATE(atoms(m1))
          END IF
          CALL parse_numbers(err_unr,strngs,nword,atoms(ntot),nats,nsevere)
          atoms(ntot)=Nats
          ntot=ntot+Nats+1
          natoms_Tot=natoms_Tot+1

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
