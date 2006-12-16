MODULE TOPPAR_STORE_MOD
  USE PARAMETERS_GLOBALS
  USE CONSTANTS
  USE NODES_Mod
  USE STRPAK
  IMPLICIT NONE 
  PRIVATE
  PUBLIC Store_Toppar,Create_Keywords,ST_GetTop,ST_GetPar, counter, Init
  CHARACTER(len=max_char), DIMENSION(3), PARAMETER :: Unused=(/'IC  ','PATC','BILD'/)
  CHARACTER(max_char), DIMENSION(:), ALLOCATABLE, SAVE :: labels
  CHARACTER(max_char), DIMENSION(:), ALLOCATABLE, SAVE :: Keywords
  INTEGER, SAVE :: count0=0
  CHARACTER(max_char) :: tok
  CHARACTER(len=*), PARAMETER :: lst=' ,()[]'
  CHARACTER(max_char), PRIVATE, SAVE :: FField,Type,Residue
  INTEGER, SAVE :: counter=0,n_labels
  INTERFACE Store_Toppar
     MODULE PROCEDURE ST_Toppar_run1
     MODULE PROCEDURE ST_Toppar_run2
     MODULE PROCEDURE ST_Toppar_run3
  END INTERFACE
CONTAINS
  SUBROUTINE Create_Keywords(line)
    CHARACTER(len=*) :: line
    CHARACTER(len=max_char), DIMENSION(4), SAVE :: Tops=(/'RESIDUE ',&
         &'PRESIDUE','RESI    ','PRES    '/) 
    CHARACTER(len=max_char), DIMENSION(8), SAVE :: Pars=(/'BONDS    ',&
         &'ANGLES   ','DIHEDRALS','IMPROPER ','NONBONDED',&
         &'BOND     ','BENDINGS ','TORSION  '/)

    IF(ALLOCATED(keywords)) DEALLOCATE(Keywords)           
    IF(TRIM(line) == 'TOPOLOGY') THEN
       ALLOCATE(Keywords(SIZE(TOPS)))
       Keywords=Tops
    ELSE IF(TRIM(line) == 'PARAMETERS') THEN
       ALLOCATE(Keywords(SIZE(Pars)))
       Keywords=Pars
    END IF
  END SUBROUTINE Create_Keywords
  SUBROUTINE Init
    count0=0
    FField=' ';Type=' ';Residue=' '
  END SUBROUTINE Init
  SUBROUTINE ST_TopPar_run1(io,end_of_line)
    IMPLICIT NONE 
    INTEGER :: io
    LOGICAL :: end_of_line
    INTEGER :: nxt,iflag0,iopt,n,c
    CHARACTER(len=max_char) :: line

    end_of_line=.FALSE.
    DO
       READ(io,'(a)',IOSTAT=iopt) line
       IF(iopt /= 0) EXIT
       nxt=0
       CALL Remove_Strange_Characters(line)
       CALL Token(0,lst,line,nxt,tok,iflag0)
       IF(iflag0 == 0) THEN
          c=COUNT(TRIM(tok) == Keywords)
          IF(c == 0) CYCLE
          count0=count0+1
          counter=count0
          RETURN
       END IF
    END DO
    count0=count0+1
    counter=count0
    end_of_line=.TRUE.
    REWIND(io)
  END SUBROUTINE ST_TopPar_run1
  SUBROUTINE ST_TopPar_run2(io,labs1,labs2,end_of_line)
    IMPLICIT NONE 
    INTEGER :: io
    CHARACTER(len=*) :: labs1
    CHARACTER(len=*) :: labs2
    LOGICAL :: end_of_line

    INTEGER :: nxt,iflg,iflag0,j1,j2,iopt,i,n,c
    CHARACTER(len=max_char) :: line,linea
    LOGICAL :: ok

    end_of_line=.FALSE.
    CALL Start()
    DO
       READ(io,'(a)',IOSTAT=iopt) line
       IF(iopt /= 0) EXIT
       CALL Remove_Strange_Characters(line)
       DO i=1,SIZE(Comms)
          CALL STR_TRIM(Comms(i),line)
       END DO
       IF(LEN_TRIM(line) == 0) CYCLE
       nxt=0
       CALL Token(0,lst,line,nxt,tok,iflag0)
       ok=.TRUE.
       DO i=1,SIZE(Unused)
          CALL STR_Fam(Unused(i),tok,j1,j2,iflg)
          IF(iflg == 0) ok=.FALSE.
       END DO
       IF(.NOT. ok) CYCLE
       
       IF(iflag0 == 0) THEN
          c=COUNT(TRIM(tok) == Keywords)
          IF(c == 0) THEN
             CALL Add(line)
             CYCLE
          END IF
          labs1=Residue
          Residue=TRIM(tok)
          CALL Extract(labels,n_labels)
          CALL Cleanup()
          CALL Start()
          labs2=Type
          CALL Token(0,lst,line,nxt,tok,iflag0)
          Type=TRIM(tok)
          IF(LEN_TRIM(labs1) == 0) THEN
             labs1='HEADER'
             labs2='MASS'
          END IF
          count0=count0+1
          counter=count0
          RETURN
       END IF
    END DO
    count0=count0+1
    counter=count0
    end_of_line=.TRUE.
    CALL Extract(labels,n_labels)
    CALL Cleanup()
    CALL Token(0,lst,line,nxt,tok,iflag0)
    labs1=TRIM(Residue)
    labs2=TRIM(Type)
  END SUBROUTINE ST_TopPar_run2
  SUBROUTINE ST_TopPar_run3(io,labs1,end_of_line)
    IMPLICIT NONE 
    INTEGER :: io
    CHARACTER(len=*) :: labs1
    LOGICAL :: end_of_line

    INTEGER :: nxt,iflg,iflag0,j1,j2,iopt,i,n,c
    CHARACTER(len=max_char) :: line,linea
    LOGICAL :: ok

    end_of_line=.FALSE.
    CALL Start()
    DO
       READ(io,'(a)',IOSTAT=iopt) line
       IF(iopt /= 0) EXIT
       CALL Remove_Strange_Characters(line)
       DO i=1,SIZE(Comms)
          CALL STR_TRIM(Comms(i),line)
       END DO
       IF(LEN_TRIM(line) == 0) CYCLE
       nxt=0
       CALL Token(0,lst,line,nxt,tok,iflag0)
       ok=.TRUE.
       DO i=1,SIZE(Unused)
          CALL STR_Fam(Unused(i),tok,j1,j2,iflg)
          IF(iflg == 0) ok=.FALSE.
       END DO
       IF(.NOT. ok) CYCLE
       
       IF(iflag0 == 0) THEN
          c=COUNT(TRIM(tok) == Keywords)
          IF(c == 0) THEN
             CALL Add(line)
             CYCLE
          END IF
          labs1=Residue
          Residue=TRIM(tok)
          CALL Extract(labels,n_labels)
          CALL Cleanup()
          CALL Start()
          IF(LEN_TRIM(labs1) == 0) THEN
             labs1='HEADER'
          END IF
          count0=count0+1
          counter=count0
          RETURN
       END IF
    END DO
    count0=count0+1
    counter=count0
    end_of_line=.TRUE.
    CALL Extract(labels,n_labels)
    CALL Cleanup()
    CALL Token(0,lst,line,nxt,tok,iflag0)
    labs1=TRIM(Residue)
  END SUBROUTINE ST_TopPar_run3
  SUBROUTINE ST_GetTop(Res,labs1,labs2)
    IMPLICIT NONE 
    TYPE(Resid), DIMENSION(:) :: Res
    CHARACTER(len=*) :: labs1,labs2
    CHARACTER(len=max_char) :: linea
    INTEGER :: i,n_end,c
    INTEGER, SAVE :: count1=0

    n_end=n_labels
    c=COUNT(labs1 == Keywords(1:2))
    IF(c /= 0) THEN
       n_end=n_labels-1
    END IF
    Res(counter)%FField='ORAC'
    c=COUNT(labs1 == Keywords(3:4))
    IF(c /= 0 .OR. labs1 == 'HEADER') Res(counter)%FField='CHARMM'

    CALL TRANLC(labs2)
    Res(counter)%Type=labs2
    Res(counter)%Residue=labs1
    ALLOCATE(Res(counter)%line(n_end))
    DO i=1,n_end
       linea=labels(i)
       CALL TRANLC(linea)
       Res(counter)%line(i)=linea
    END DO
    DEALLOCATE(labels)
    count1=count1+1
  END SUBROUTINE ST_GetTop
  SUBROUTINE ST_GetPar(Res,labs1)
    IMPLICIT NONE 
    TYPE(Resid), DIMENSION(:) :: Res
    CHARACTER(len=*) :: labs1
    CHARACTER(len=max_char) :: linea
    INTEGER :: i,n_end,c

    n_end=n_labels
    Res(counter)%FField='ORAC'
    c=COUNT(labs1 == Keywords(1:4))
    IF(c /= 0 .OR. labs1 == 'HEADER') Res(counter)%FField='CHARMM'
    Res(counter)%Type=' '
    Res(counter)%Residue=labs1
    ALLOCATE(Res(counter)%line(n_end))
    DO i=1,n_end
       linea=labels(i)
       CALL TRANLC(linea)
       Res(counter)%line(i)=linea
    END DO
    DEALLOCATE(labels)
  END SUBROUTINE ST_GetPar
  SUBROUTINE Remove_Strange_Characters(line)
    CHARACTER(len=*) :: line
    CHARACTER(len=1), DIMENSION(3), SAVE :: ch=(/CHAR(9),CHAR(11),CHAR(27)/)
    CHARACTER(len=60) :: pippa
    INTEGER :: n 
    DO n=1,SIZE(ch)
       CALL CHR_SAR(ch(n),' ',line)
    END DO
  END SUBROUTINE Remove_Strange_Characters
END MODULE TOPPAR_STORE_MOD
