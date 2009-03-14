  SUBROUTINE Read_it(knlist,kprint,nsevere,nword,strngs,iret,errmsg&
       &,read_err)
!!$***********************************************************************
!!$   Time-stamp: <2006-02-09 14:20:32 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Feb  9 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program  ----*
    
!!$======================== DECLARATIONS ================================*
    
    IMPLICIT none
    
!!$----------------------------- ARGUMENTS ------------------------------*

    INTEGER :: knlist,kprint,nword,nsevere,iret,read_err
    CHARACTER(80) :: errmsg
    CHARACTER(80) :: strngs(:)

!!$------------------------- LOCAL VARIABLES ----------------------------*

    TYPE NODES
       CHARACTER(22) :: label
       CHARACTER(22), DIMENSION(:), POINTER :: data
       TYPE (NODES), POINTER :: next
    END TYPE NODES
    TYPE(NODES), POINTER, SAVE :: list,nfirst

    CHARACTER(80) :: line,aux,filename
    INTEGER :: count=0,i,m

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

!!$
!!$==== Command GEOMETRY  ================================================
!!$

    count=0
    line(79:80)='  '
    read_err=0
    Geom_groups=.TRUE.

    ALLOCATE(list)
    NULLIFY(list%next)
    nfirst=>list

    filename='none'
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

       CASE('dist_mass')
          IF(nword .LT. 3) THEN
             errmsg=err_args(1)//'3'
             CALL xerror(errmsg,80,1,30)
             nsevere=nsevere+1
          ELSE
             count=count+1
             list%label=TRIM(strngs(1))
             ALLOCATE(list%data(nword-1))
             list%data(:)=strngs(2:nword)

             ALLOCATE(list%next)
             list=>list%next
             NULLIFY(list%next)
          END IF
       CASE('dist')
          IF(nword .LT. 3) THEN
             errmsg=err_args(1)//'3'
             CALL xerror(errmsg,80,1,30)
             nsevere=nsevere+1
          ELSE
             count=count+1
             list%label=TRIM(strngs(1))
             ALLOCATE(list%data(nword-1))
             list%data(:)=strngs(2:nword)

             ALLOCATE(list%next)
             list=>list%next
             NULLIFY(list%next)
          END IF
       CASE('write')
          IF(nword == 1) THEN
             filename='GEOMETRY.dat'
             CALL openf(k_write,filename,'F','U',0)
          ELSE IF(nword == 2) THEN
             CALL Read_String(strngs(2),n_write)
             filename='GEOMETRY.dat'
             CALL openf(k_write,filename,'F','U',0)
          ELSE IF(nword == 3) THEN
             CALL Read_String(strngs(2),n_write)
             filename=TRIM(strngs(3))
             CALL openf(k_write,filename,'F','U',0)
          END IF

       CASE(' ')
          CYCLE
       CASE('END')
          EXIT
       END SELECT
    END DO
    IF(TRIM(filename) == 'none') THEN
       filename='GEOMETRY.dat'
       CALL openf(k_write,filename,'F','U',0)
    END IF

    ALLOCATE(Geoms(count))
    i=0
    list=>nfirst
    DO WHILE (ASSOCIATED(list%next))
       i=i+1
       Geoms(i)%type=list%label
       m=SIZE(list%data)
       ALLOCATE(Geoms(i)%data(m))
       Geoms(i)%data=list%data
       list => list%next
    END DO
    CALL Cleanup
    RETURN

600 read_err=1
    RETURN

20  CONTINUE
    iret=1
    errmsg='internal reading error: wrong format?? TAB character??'
    CALL xerror(errmsg,80,1,2)
    
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

  CONTAINS
  SUBROUTINE Cleanup
    IMPLICIT none
    TYPE(NODES), POINTER :: dummy

    list=>nfirst
    DO WHILE (ASSOCIATED(list%next))
       dummy=>list%next
       DEALLOCATE(list%data)
       DEALLOCATE(list)
       list => dummy
    END DO
    NULLIFY(dummy,nfirst)

  END SUBROUTINE Cleanup
END SUBROUTINE Read_it
