MODULE PDB_Utils

!!$***********************************************************************
!!$   Time-stamp: <2006-12-20 21:54:38 marchi>                           *
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

  USE CONSTANTS, ONLY: max_pars,max_data,max_char, Used
  USE ERROR_Mod, ONLY: Add_Errors=>Add, Print_Errors, error_args, errmsg_f
  USE STRINGS_Mod, ONLY: My_Fam
  USE Linked_Char
  IMPLICIT none
  PRIVATE
  PUBLIC StorePDB,Coordinates,ResiduePDB, AtomPdb
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
  SUBROUTINE StorePDB(Type,line,PDB_)
    CHARACTER(len=*) :: Type,line
    CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE :: PDB_
    LOGICAL :: ok,end_of_list
    CHARACTER(len=max_char) :: linea
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
    CALL Start()
    DO 
       READ(unit=io,fmt='(a80)',IOSTAT=iopt) linea
       IF(iopt /= 0) EXIT
       ok=.FALSE.
       DO n=1,SIZE(Used)
          IF(My_Fam(TRIM(Used(n)),linea)) ok=.TRUE.
       END DO
       IF(.NOT. ok) CYCLE
       CALL Add(linea, count_out)
    END DO
    ALLOCATE(PDB_(count_out))
    count_a=0
    end_of_list=.FALSE.
    CALL Extract(linea, end_of_list, 0)
    DO WHILE(.NOT. end_of_list) 
       count_A=count_A+1
       PDB_(count_a)=ADJUSTL(TRIM(linea))
       CALL Extract(linea, end_of_list)
    END DO
    CALL Cleanup()
    CLOSE(io)
  END SUBROUTINE StorePDB
  SUBROUTINE Coordinates(Type,PDB_,ResPdb)
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
    CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE :: PDB_
    TYPE(ResiduePDB), DIMENSION(:), ALLOCATABLE :: ResPdb

    INTEGER :: n,m, Id_Res0,n_Last
    CHARACTER(len=max_char) :: Res0
    INTEGER, DIMENSION(:), ALLOCATABLE :: Ids
    INTEGER :: Serial,ResSeq,Id_Old,Id_New,count_out,count_a
    REAL(8) :: x,y,z
    CHARACTER(len=5) :: AtmName,ResName
    CHARACTER(len=max_char)  :: linea
    LOGICAL :: end_of_list


    n=1
    DO WHILE(PDB_(n)(1:6) /= 'ATOM  ' .AND. PDB_(n)(1:6) /= 'HETATM')
       IF(n  == SIZE(PDB_)) THEN
          errmsg_f=TRIM(Type)//' .PDB File does not contains coordinates'
          CALL Add_Errors(-1,errmsg_f)
          RETURN
       END IF
       n=n+1
    END DO
    READ(PDB_(n),'(17x,a3)' ) Res0
    READ(PDB_(n),'(22x,i4)' ) Id_Res0

    m=0
    id_old=Id_Res0
    DO n=1,SIZE(PDB_)
       IF(PDB_(n)(1:6) == 'ATOM  ' .OR. PDB_(n)(1:6) == 'HETATM') THEN
          READ(PDB_(n),'(22x,i4)' ) Id_New
          IF(Id_New /= Id_old) THEN
             m=m+1
             Id_old=Id_new
          END IF
       END IF
    END DO

    ALLOCATE(ResPDB(m+1))

    m=0
    CALL Start()
    count_out=0
    Id_old=Id_Res0
    DO n=1,SIZE(PDB_)
       IF(PDB_(n)(1:6) == 'ATOM  ' .OR. PDB_(n)(1:6) == 'HETATM') THEN
          READ(PDB_(n),'(22x,i4)' ) Id_New

          IF(Id_New /= Id_old) THEN
             m=m+1
             ALLOCATE(ResPdb(m) % atm(count_out))
             CALL Get_it
             CALL Start()
             CALL Add(PDB_(n), count_out)
             Id_Old=Id_New
          ELSE
             CALL Add(PDB_(n), count_out)
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
      count_A=0
      end_of_list=.FALSE.
      CALL Extract(linea, end_of_list, 0)
      DO WHILE(.NOT. end_of_list) 
         count_A=count_A+1
         READ(linea,'(6x,i5,1x,a4,1x,a3,2x,i4,4x,3f8.3)' ) &
              & Serial,AtmName,ResName,ResSeq,x,y,z
         ResPdb(m) % atm(count_a) % x=x
         ResPdb(m) % atm(count_a) % y=y
         ResPdb(m) % atm(count_a) % z=z
         ResPdb(m) % atm(count_a) % Serial=Serial
         ResPdb(m) % atm(count_a) % AtmName=AtmName
         ResPdb(m) % No = ResSeq
         ResPdb(m) % ResName = ResName
         CALL Extract(linea, end_of_list)
      END DO
      CALL Cleanup()
    END SUBROUTINE Get_it
  END SUBROUTINE Coordinates
  SUBROUTINE Verify()
  END SUBROUTINE Verify
END MODULE PDB_Utils
