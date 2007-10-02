MODULE PDB

!!$***********************************************************************
!!$   Time-stamp: <2007-01-04 17:42:34 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Wed Dec 20 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************
!!$---- This subroutine is part of the program ORAC ----*
!!$======================== DECLARATIONS ================================*

  USE Constants, ONLY: max_pars,max_data,max_char, Used
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, error_args, errmsg_f
  USE Strings, ONLY: My_Fam
  USE Node
  IMPLICIT none
  PRIVATE
  PUBLIC PDB_,PDB__Coordinates,ResiduePDB, AtomPdb
  TYPE :: AtomPDB
     INTEGER :: Serial
     REAL(8) :: x,y,z
     CHARACTER(len=4) :: AtmName
  END TYPE AtomPDB
  TYPE ResiduePDB
     INTEGER :: No
     CHARACTER(len=3) :: ResName
     TYPE(AtomPdb), DIMENSION(:), ALLOCATABLE :: atm
  END TYPE ResiduePDB
CONTAINS
  SUBROUTINE PDB_(Type,line,PDB_string)
    CHARACTER(len=*) :: Type,line
    CHARACTER(len=max_char), DIMENSION(:), POINTER :: PDB_string
    LOGICAL :: ok,end_of_list
    CHARACTER(len=max_char) :: linea
    CHARACTER(len=max_char), DIMENSION(:), POINTER :: line1=>NULL()
    INTEGER :: n,io, count_out,count_a,iopt
    
    
    CALL Channel(io)
    INQUIRE(file=line,EXIST=ok)
    IF(ok) THEN
       OPEN(unit=io,file=line,form='FORMATTED',status='OLD')
    ELSE
       errmsg_f=TRIM(Type)//' .PDB file '''//TRIM(line)//''' does not exist'
       CALL Add_Errors(-1,errmsg_f)
       RETURN
    END IF

    IF(.NOT. Node_()) STOP
    DO 
       READ(unit=io,fmt='(a80)',IOSTAT=iopt) linea
       IF(iopt /= 0) EXIT
       ok=.FALSE.
       DO n=1,SIZE(Used)
          IF(My_Fam(TRIM(Used(n)),linea)) ok=.TRUE.
       END DO
       IF(.NOT. ok) CYCLE
       CALL Node__Push(linea)
    END DO
    count_out=Node__Size()
    ALLOCATE(PDB_string(count_out))
    
    count_a=0
    DO WHILE(Node__Pop(line1))
       count_A=count_A+1
       PDB_string(count_a)=ADJUSTL(TRIM(line1(1)))
    END DO
    CALL Node__Delete()

    CLOSE(io)
  END SUBROUTINE PDB_
  SUBROUTINE PDB__Coordinates(Type,PDB_string,ResPdb)
!!$******************************************************************************
!!$PDB Format Description Version 2.2 from http://www.rcsb.org
!!$COLUMNS       DATA TYPE       FIELD         DEFINITION
!!$-------------------------------------------------------------------------------
!!$ 1 -  6       Record name     "ATOM  "
!!$ 7 - 11       Integer         serial        Atom serial number.
!!$13 - 16       Atom            name          Atom name.
!!$17            Character       altLoc        Alternate location indicator.
!!$18 - 20       Residue name    resName       Residue name.
!!$22            Character       chainID       Chain identifier.
!!$23 - 26       Integer         resSeq        Residue sequence number.
!!$27            AChar           iCode         Code for insertion of residues.
!!$31 - 38       Real(8.3)       x             Orthogonal coordinates for X in
!!$                                            Angstroms.
!!$39 - 46       Real(8.3)       y             Orthogonal coordinates for Y in
!!$                                            Angstroms.
!!$47 - 54       Real(8.3)       z             Orthogonal coordinates for Z in
!!$                                            Angstroms.
!!$55 - 60       Real(6.2)       occupancy     Occupancy.
!!$61 - 66       Real(6.2)       tempFactor    Temperature factor.
!!$73 - 76       LString(4)      segID         Segment identifier, left-justified.
!!$77 - 78       LString(2)      element       Element symbol, right-justified.
!!$79 - 80       LString(2)      charge        Charge on the atom.
!!$-------------------------------------------------------------------------------
!!$-------------------------------------------------------------------------------
    CHARACTER(len=*) :: Type
    CHARACTER(len=max_char), DIMENSION(:), POINTER :: PDB_string

    TYPE(ResiduePDB), DIMENSION(:), POINTER :: ResPdb

    INTEGER :: n,m, Id_Res0,n_Last
    CHARACTER(len=max_char) :: Res0
    INTEGER, DIMENSION(:), ALLOCATABLE :: Ids
    INTEGER :: Serial,ResSeq,Id_Old,Id_New,count_out,count_a
    REAL(8) :: x,y,z
    CHARACTER(len=5) :: AtmName,ResName
    LOGICAL :: end_of_list


    n=1
    DO WHILE(PDB_string(n)(1:6) /= 'ATOM  ' .AND. PDB_string(n)(1:6) /= 'HETATM')
       IF(n  == SIZE(PDB_string)) THEN
          errmsg_f=TRIM(Type)//' .PDB File does not contains coordinates'
          CALL Add_Errors(-1,errmsg_f)
          RETURN
       END IF
       n=n+1
    END DO
    READ(PDB_string(n),'(17x,a3)' ) Res0
    READ(PDB_string(n),'(22x,i4)' ) Id_Res0

    m=0
    id_old=Id_Res0
    DO n=1,SIZE(PDB_string)
       IF(PDB_string(n)(1:6) == 'ATOM  ' .OR. PDB_string(n)(1:6) == 'HETATM') THEN
          READ(PDB_string(n),'(22x,i4)' ) Id_New
          IF(Id_New /= Id_old) THEN
             m=m+1
             Id_old=Id_new
          END IF
       END IF
    END DO

    ALLOCATE(ResPDB(m+1))

    m=0
    IF(.NOT. Node_()) STOP
    
    count_out=0
    Id_old=Id_Res0
    DO n=1,SIZE(PDB_string)
       IF(PDB_string(n)(1:6) == 'ATOM  ' .OR. PDB_string(n)(1:6) == 'HETATM') THEN
          READ(PDB_string(n),'(22x,i4)' ) Id_New

          IF(Id_New /= Id_old) THEN
             m=m+1
             ALLOCATE(ResPdb(m) % atm(count_out))
             CALL Get_it
             IF(.NOT. Node_()) STOP
             CALL Node__Push(PDB_string(n))
             Id_Old=Id_New
          ELSE
             CALL Node__Push(PDB_string(n))
          END IF
       END IF
    END DO
    m=m+1
    ALLOCATE(ResPdb(m) % atm(count_out))
    CALL Get_it
100 FORMAT(A6,I5,1X,A4,1X,A3,2X,I4,4X,3F8.3,2F6.2)
  CONTAINS
    SUBROUTINE Get_it
      INTEGER count_a
      CHARACTER(len=max_char), DIMENSION(:), POINTER :: line1=>NULL()

      count_A=0
      DO WHILE(Node__Pop(line1))
         count_A=count_A+1
         READ(line1(1),'(6x,i5,1x,a4,1x,a3,2x,i4,4x,3f8.3)' ) &
              & Serial,AtmName,ResName,ResSeq,x,y,z
         ResPdb(m) % atm(count_a) % x=x
         ResPdb(m) % atm(count_a) % y=y
         ResPdb(m) % atm(count_a) % z=z
         ResPdb(m) % atm(count_a) % Serial=Serial
         ResPdb(m) % atm(count_a) % AtmName=AtmName
         ResPdb(m) % No = ResSeq
         ResPdb(m) % ResName = ResName
      END DO
      CALL Node__Delete()
    END SUBROUTINE Get_it
  END SUBROUTINE PDB__Coordinates
  SUBROUTINE PDB__Verify()
  END SUBROUTINE PDB__Verify
END MODULE PDB
