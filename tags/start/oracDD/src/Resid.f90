MODULE Resid

!!$***********************************************************************
!!$   Time-stamp: <2007-01-11 14:44:04 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Jan  4 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORACS ----*


!!$======================== DECLARATIONS ================================*
  
  USE Constants
  USE Node
  USE Keyword
  USE Errors, ONLY: Add_Errors=>Add, error_other, errmsg_f, Print_Errors
  USE MyParse
  USE Strings, ONLY: MY_Fam
  USE STRPAK
  IMPLICIT none
  PRIVATE
  PUBLIC Resid_, Resids, Topology, Paras

  TYPE Resids
     CHARACTER(len=max_char) :: Type,Residue,FField
     CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE :: line
  END TYPE Resids
  CHARACTER(len=max_char), DIMENSION(:), POINTER :: Input_Type
  
  CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE :: data
  CHARACTER(len=max_char), DIMENSION(5), PARAMETER :: Unused=(/'IC  ','PATC','BILD'&
       &,'END ','HBON'/)
  TYPE(Keyword__Type) :: Keys
  TYPE(Resids), DIMENSION(:), ALLOCATABLE, TARGET, SAVE  :: Topology,Paras
CONTAINS
  SUBROUTINE Resid_(linea)
    CHARACTER(len=*) :: linea
    CHARACTER(len=max_char) :: line,line1,Type
    LOGICAL :: end_of_file=.FALSE.,exist
    INTEGER :: io,iopt,i,nword,ndata,count_a,n,ntype,k,m
    LOGICAL :: ok

    CALL CHANNEL(io)
    INQUIRE(FILE=linea,EXIST=exist)
    IF(.NOT. exist) THEN
       errmsg_f=TRIM(error_other % g(1))//' : '//TRIM(linea)//' does not exist'
       CALL Add_Errors(-1,errmsg_f)
       CALL Print_Errors()
       RETURN
    END IF
    OPEN(unit=io,file=linea,form='FORMATTED')
    IF(.NOT. Node_()) STOP
    DO 
       READ(io,'(a)',IOSTAT=iopt) line
       IF(iopt /= 0) EXIT
       CALL Remove(line)

       DO i=1,SIZE(Comms)
          CALL STR_TRIM(Comms(i),line)
       END DO
       IF(LEN_TRIM(line) == 0) CYCLE
       nword=Myparse_(line)
       ok=.TRUE.
       DO i=1,SIZE(Unused)
          IF(My_Fam(Unused(i),strngs(1))) ok=.FALSE.
       END DO
       IF(.NOT. ok) CYCLE
       CALL Node__Push(line)
    END DO
    ndata=Node__Size()

    ALLOCATE(data(ndata))
    data=' '
    count_a=0
    DO WHILE(Node__Pop(line1))
       count_a=count_a+1
       data(count_a)=line1
    END DO

    CALL Keyword_(data, linea)

    Keys = Keyword__Pop() 
    n=0
    Type = 'tops'
    DO WHILE(.NOT. Keys % Finish )
       n=n+1
       IF(n > 1) THEN
          IF(LEN_TRIM(Keys % Type) == 0) THEN
             Type = 'pars'
          END IF
          EXIT
       END IF
       Keys = Keyword__Pop() 
    END DO

    CALL Keyword__PopReset    
    ntype=Keyword__Size()
    IF(TRIM(Type) == 'tops') THEN
       ALLOCATE(Topology(ntype))
       Keys = Keyword__Pop()
       n=0
       DO WHILE(.NOT. Keys % Finish )
          n=n+1
          ALLOCATE(Topology(n) % line (Keys % End - Keys % Begin + 1))
          Topology (n) % line = data(Keys % Begin: Keys % End)
          Topology (n) % Type = Keys % Type
          Topology (n) % Residue = Keys % Residue 
          Topology (n) % FField = Keys % FField
          Keys = Keyword__Pop()          
       END DO
       DO k=1,SIZE(Topology)
          DO m=1,SIZE(Topology (k) % line)
             CALL TRANLC(Topology (k)% line(m))
          END DO
       END DO
!!$       DO k=1,SIZE(Topology)
!!$          WRITE(*,*) TRIM(Topology  (k)% Type),' ',TRIM(Topology  (k)% Residue)&
!!$               &,' ',TRIM(Topology  (k)% FField)
!!$          DO m=1,SIZE(Topology (k) % line)
!!$             WRITE(*,*) TRIM(Topology (k)% line(m))
!!$          END DO
!!$       END DO

    ELSE IF(TRIM(Type) == 'pars') THEN
       ALLOCATE(Paras(ntype))
       Keys = Keyword__Pop()
       n=0
       DO WHILE(.NOT. Keys % Finish )
          n=n+1
          ALLOCATE(Paras (n) % line (Keys % End - Keys % Begin + 1))
          Paras  (n) % line = data(Keys % Begin: Keys % End)
          Paras (n) % Type = Keys % Type
          Paras (n) % Residue = Keys % Residue 
          Paras (n) % FField = Keys % FField

          Keys = Keyword__Pop()          
       END DO
!!$       DO k=1,SIZE(Paras)
!!$          WRITE(*,*) 'x ',TRIM(Paras (k)% Type),' x ',TRIM(Paras  (k)% Residue)&
!!$               &,' x ',TRIM(Paras  (k)% FField),' x '
!!$          DO m=1,SIZE(paras(k) % line)
!!$             WRITE(*,*) TRIM(Paras (k)% line(m))
!!$          END DO
!!$       END DO
       DO k=1,SIZE(Paras)
          DO m=1,SIZE(paras(k) % line)
             CALL TRANLC(Paras (k)% line(m))
          END DO
       END DO
       
    END IF
    
    DEALLOCATE(data)
  CONTAINS
    SUBROUTINE Remove(line)
      CHARACTER(len=*) :: line
      CHARACTER(len=1), DIMENSION(3), SAVE :: ch=(/CHAR(9),CHAR(11),CHAR(27)/)
      INTEGER :: n 
      DO n=1,SIZE(ch)
         CALL CHR_SAR(ch(n),' ',line)
      END DO
    END SUBROUTINE Remove
    
  END SUBROUTINE Resid_
  

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE Resid
