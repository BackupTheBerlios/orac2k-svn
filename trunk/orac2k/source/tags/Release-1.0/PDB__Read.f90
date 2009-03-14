!!$/---------------------------------------------------------------------\
!!$   Copyright  © 2006-2007 Massimo Marchi <Massimo.Marchi at cea.fr>   |
!!$                                                                      |
!!$    This software is a computer program named oracDD whose            |
!!$    purpose is to simulate and model complex molecular systems.       |
!!$    The code is written in fortran 95 compliant with Technical        |
!!$    Report TR 15581, and uses MPI-1 routines for parallel             |
!!$    coding.                                                           |
!!$                                                                      |
!!$    This software is governed by the CeCILL license under             |
!!$    French law and abiding by the rules of distribution of            |
!!$    free software.  You can  use, modify and/ or redistribute         |
!!$    the software under the terms of the CeCILL icense as              |
!!$    circulated by CEA, CNRS and INRIA at the following URL            |
!!$    "http://www.cecill.info".                                         |
!!$                                                                      |
!!$    As a counterpart to the access to the source code and rights      |
!!$    to copy, modify and redistribute granted by the license,          |
!!$    users are provided only with a limited warranty and the           |
!!$    software's author, the holder of the economic rights, and         |
!!$    the successive licensors have only limited liability.             |
!!$                                                                      |
!!$    The fact that you are presently reading this means that you       |
!!$    have had knowledge of the CeCILL license and that you accept      |
!!$    its terms.                                                        |
!!$                                                                      |
!!$    You should have received a copy of the CeCILL license along       |
!!$    with this program; if not, you can collect copies on the URL's    |
!!$    "http://www.cecill.info/licences/Licence_CeCILL_V2-en.html"       |
!!$    "http://www.cecill.info/licences/Licence_CeCILL_V2-fr.html"       |
!!$                                                                      |
!!$----------------------------------------------------------------------/
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

!!$
!!$==== Command DENSITY  ================================================
!!$

    line(79:80)='  '
    read_err=0
    PDB_=.TRUE.
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

       CASE('template')
          IF(nword .LT. 2) THEN
             errmsg=err_args(1)//'1'
             CALL xerror(errmsg,80,1,30)
             nsevere=nsevere+1
          ELSE
             PDB_Filename=strngs(2)
          END IF
       CASE(' ')
          CYCLE

       CASE('END')
          EXIT
       END SELECT
    END DO
    INQUIRE(FILE=PDB_filename,EXIST=exist)
    IF(exist) THEN
       CALL openf(kpdb,PDB_filename,'FORMATTED','OLD',0)
    ELSE
       errmsg='PDB Template file does not exist '
       CALL xerror(errmsg,80,1,30)
       nsevere = nsevere + 1
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
  FUNCTION PDB__Read(Type,PDB_string) RESULT(out)
!!$---- This subroutine is part of the program oracDD ----*
!!$****************************************************************************
!!$PDB Format Description Version 2.2 from http://www.rcsb.org
!!$COLUMNS       DATA TYPE       FIELD         DEFINITION
!!$----------------------------------------------------------------------------
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
!!$73 - 76       LString(4)      segID         Segment identifier,
!!$  left-justified.
!!$77 - 78       LString(2)      element       Element symbol, right-justified.
!!$79 - 80       LString(2)      charge        Charge on the atom.
!!$----------------------------------------------------------------------------
!!$----------------------------------------------------------------------------

    LOGICAL :: out
    CHARACTER(len=*) :: Type
    CHARACTER(len=max_char)  :: PDB_string(:)

    INTEGER :: n,m, Id_Res0,n_Last
    CHARACTER(len=max_char) :: Res0,labs1
    INTEGER, ALLOCATABLE :: Ids(:)
    INTEGER :: Serial,ResSeq,Id_Old,Id_New,count_out,count_a
    REAL(8) :: x,y,z,a0,b0,c0,alpha0,beta0,gamma0
    CHARACTER(len=5) :: AtmName,ResName,chain
    LOGICAL :: end_of_list
    REAL(8), PARAMETER :: eps=5.0D-2
    CHARACTER(len=max_Char) :: labs0

    out=.TRUE.
    n=1
    DO WHILE(PDB_string(n)(1:6) /= 'ATOM  ' .AND. PDB_string(n)(1:6) /= 'HETATM')
       IF(n  == SIZE(PDB_string)) THEN
          errmsg_f=TRIM(Type)//' .PDB File does not contain coordinates'
          CALL Add_Errors(-1,errmsg_f)
          out=.FALSE.
          RETURN
       END IF
       n=n+1
    END DO

    Res0=PDB_string(n)(18:20)
    CALL TRANLC(Res0)

    labs1=PDB_string(n)(23:26)
    CALL MyRead(labs1,Id_Res0)
    m=0
    id_old=Id_Res0
    DO n=1,SIZE(PDB_string)
       IF(PDB_string(n)(1:5) == 'ATOM ' .OR. PDB_string(n)(1:5) == 'HETAT') THEN

          labs1=PDB_string(n)(23:26)
          CALL MyRead(labs1,Id_New)
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
       IF(PDB_string(n)(1:5) == 'ATOM ' .OR. PDB_string(n)(1:5) == 'HETAT') THEN

          IF(Its_a_Number(PDB_string(n)(22:22))) THEN
             labs1=PDB_string(n)(22:26)
          ELSE
             labs1=PDB_string(n)(23:26)
          END IF

          CALL MyRead(labs1,Id_New)
          IF(Id_New /= Id_old) THEN
             m=m+1
             count_out=Node__Size()
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
    count_out=Node__Size()
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

         IF(Its_a_Number(line1(1)(6:6))) THEN
            labs1=line1(1)(6:11)
         ELSE
            labs1=line1(1)(7:11)
         END IF
         CALL Myread(labs1,Serial)

         AtmName=line1(1)(13:16)
         ResName=line1(1)(18:20)
         Chain=line1(1)(22:22)

         CALL TRANLC(AtmName)
         CALL TRANLC(ResName)

         IF(Its_a_Number(line1(1)(22:22))) THEN 
            labs1=Line1(1)(22:26)
         ELSE
            labs1=Line1(1)(23:26)
         END IF
         CALL Myread(labs1,ResSeq)
         labs1=ADJUSTL(Line1(1)(31:38))
         CALL Myread(labs1,x)

         labs1=ADJUSTL(Line1(1)(39:46))
         CALL Myread(labs1,y)

         labs1=ADJUSTL(Line1(1)(47:54))
         CALL Myread(labs1,z)
         ResPdb(m) % atm(count_a) % x=x
         ResPdb(m) % atm(count_a) % y=y
         ResPdb(m) % atm(count_a) % z=z
         ResPdb(m) % atm(count_a) % Serial=Serial
         ResPdb(m) % atm(count_a) % AtmName=ADJUSTL(AtmName)
         ResPdb(m) % No = ResSeq
         ResPdb(m) % Chain = ADJUSTL(Chain)
         ResPdb(m) % ResName = ADJUSTL(ResName)
      END DO
      CALL Node__Delete()
    END SUBROUTINE Get_it
    FUNCTION Its_a_Number(char) RESULT(out)
      CHARACTER(len=1) :: char
      LOGICAL :: out
      out=.FALSE.
      IF(ICHAR(char) >= 49 .AND. ICHAR(char) <= 57) out=.TRUE.
    END FUNCTION Its_a_Number
      
  END FUNCTION PDB__Read
