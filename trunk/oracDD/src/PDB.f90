MODULE PDB

!!$***********************************************************************
!!$   Time-stamp: <2007-01-12 20:39:41 marchi>                           *
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

  USE IndSequence
  USE SystemTpg
  USE Constants, ONLY: max_pars,max_data,max_char, Used
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_w, errmsg_f
  USE Strings, ONLY: My_Fam, MyRead, MyPutNum
  USE Node
  USE MyParse
  IMPLICIT none
  PRIVATE
  PUBLIC PDB_, ResiduePDB, AtomPdb, PDB__Validate
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

  INTEGER, SAVE, POINTER :: Res_Atm(:,:)=>NULL()
  INTEGER, SAVE, POINTER :: Grp_Atm(:,:)=>NULL()
  INTEGER, SAVE, POINTER :: SltSlv(:,:)=>NULL()
  CHARACTER(len=max_char), PARAMETER :: water(8)=(/'wat  ','hoh  ','tip3 '&
       &,'tip3p','tip3f','spc  ','hoh  ','spce '/)
  CHARACTER(len=max_char), PARAMETER :: histidine(5)=(/'his','hsd','hse','hsp','hs2'/)
CONTAINS
  FUNCTION PDB_(Type,PDB_string,ResPdb) RESULT(out)
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

    TYPE(ResiduePDB), POINTER :: ResPdb(:)

    INTEGER :: n,m, Id_Res0,n_Last
    CHARACTER(len=max_char) :: Res0,labs1
    INTEGER, ALLOCATABLE :: Ids(:)
    INTEGER :: Serial,ResSeq,Id_Old,Id_New,count_out,count_a
    REAL(8) :: x,y,z
    CHARACTER(len=5) :: AtmName,ResName
    LOGICAL :: end_of_list

    out=.TRUE.
    n=1
    DO WHILE(PDB_string(n)(1:6) /= 'ATOM  ' .AND. PDB_string(n)(1:6) /= 'HETATM')
       IF(n  == SIZE(PDB_string)) THEN
          errmsg_f=TRIM(Type)//' .PDB File does not contains coordinates'
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
       IF(PDB_string(n)(1:6) == 'ATOM  ' .OR. PDB_string(n)(1:6) == 'HETATM') THEN

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
       IF(PDB_string(n)(1:6) == 'ATOM  ' .OR. PDB_string(n)(1:6) == 'HETATM') THEN

          labs1=PDB_string(n)(23:26)
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
         labs1=line1(1)(7:11)
         CALL Myread(labs1,Serial)

         AtmName=line1(1)(13:16)
         ResName=line1(1)(18:20)

         CALL TRANLC(AtmName)
         CALL TRANLC(ResName)


         labs1=Line1(1)(23:26)
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
         ResPdb(m) % ResName = ADJUSTL(ResName)
      END DO
      CALL Node__Delete()
    END SUBROUTINE Get_it
  END FUNCTION PDB_
  FUNCTION PDB__Validate(Type,ResPdb) RESULT(out)
    CHARACTER(len=*) :: Type
    TYPE(ResiduePDB), POINTER :: ResPdb(:)
    LOGICAL :: out
    LOGICAL, ALLOCATABLE :: oks(:)
    INTEGER :: n,Begins,Ends,m,p,ip_Res,oo,nword
    LOGICAL :: ok
    CHARACTER(len=max_char) :: lab_p,lab_m,res_n,lab_m0,res_m
    CHARACTER(len=max_char), POINTER  :: new_name(:)
    
    Res_Atm=>IndSequence__Res()
    SltSlv=>IndSequence__sltslv_Res()
    ip_res=1
    IF(MY_Fam('solvent',type)) ip_res=2

    out=.TRUE.
    IF(.NOT. ALLOCATED(Tpg % atm)) THEN
       errmsg_f='Cannot validate PDB coordinates if system&
            & topology is not defined'
       CALL Add_Errors(-1,errmsg_f)
       out=.FALSE.
       RETURN
    END IF

    IF(SltSlv(2,ip_Res)-SltSlv(1,ip_Res)+1 /= SIZE(ResPdb)) THEN
       errmsg_f='No. of residues in '//TRIM(Type)//' .pdb file&
            & does not match system topology: Expected '&
            &//TRIM(MyPutnum(SIZE(Res_Atm)))//' found '//TRIM(Myputnum(SIZE(ResPdb)))
       CALL Add_Errors(-1,errmsg_f)
       out=.FALSE.
       RETURN       
    END IF

    ALLOCATE(oks(SIZE(Tpg % atm)))
    oks=.FALSE.

    DO n=SltSlv(1,ip_res),SltSlv(2,ip_res)
       Begins=Res_Atm(1,n)
       Ends  =Res_Atm(2,n)
       res_n=ADJUSTL(Tpg % atm(Begins) % a % Res)

       nword=MyParse_(res_n)
       IF(TRIM(strngs(1)) == 'Link') THEN
          res_n=ADJUSTL(strngs(2))
       END IF
       IF(.NOT. My_Fam(ResPdb(n) % ResName, Res_n)) THEN
          IF(.NOT. ResName_Exceptions(ResPdb(n) % ResName,&
               & res_n)) THEN
             errmsg_w='System Topology residue type '&
                  &//TRIM(res_n)//' does not match &
                  &with .PDB file residue type '&
                  &//ResPdb(n) % ResName//' run continues,&
                  & we hope this is stil ok'
             CALL Add_Errors(1,errmsg_w)
          END IF
       END IF

       DO m=Begins,Ends
          IF(Tpg % atm(m) % a % beta(1:1) ==  'h') oks(m)=.TRUE.
       END DO

       DO m=1,SIZE(ResPdb(n) % atm)
          lab_m0=ResPdb(n) % atm(m) % AtmName
          res_m=ADJUSTL(ResPdb(n) % ResName)
          CALL AtomName_Exceptions(lab_m0, res_m, new_name)             
          ok=.FALSE.
          DO oo=1,SIZE(new_name)
             lab_m=new_name(oo)
             DO p=Begins,Ends
                lab_p=Tpg % atm (p) % a % beta
                IF(TRIM(lab_m) == TRIM(Lab_p)) THEN
                   ok=.TRUE.
                   oks(p)=.TRUE.
                   
                END IF
             END DO
          END DO
          IF(.NOT. ok) THEN
             errmsg_f='While reading coordinates for residue No. '&
                  &//TRIM(Myputnum(n))//', .PDB atom Label '//TRIM(lab_m)&
                  &//' was not found in the topoloy file&
                  & list of residue: '//TRIM(res_n)
             CALL Add_Errors(-1,errmsg_f)
             out=.FALSE.
          END IF
       END DO
       DO m=Begins,Ends
          lab_m=Tpg % atm (m) % a % beta
          IF(.NOT. oks(m)) THEN
             errmsg_f='While reading .PDB coordinates for residue No. '&
                  &//TRIM(Myputnum(n))//', cannot found label '//TRIM(lab_m)&
                  &//' of residue '//TRIM(res_n)
             CALL Add_Errors(-1,errmsg_f)
             out=.FALSE.
          END IF
       END DO
       CALL Print_Errors()
    END DO
    DEALLOCATE(oks)
  CONTAINS
    SUBROUTINE AtomName_Exceptions(lab, res, new_lab)
      CHARACTER(len=*) :: lab,res
      CHARACTER(len=max_char), POINTER :: new_lab(:)

      IF(ASSOCIATED(new_lab)) DEALLOCATE(new_lab)
      IF(TRIM(res) == 'ile') THEN
         IF(TRIM(lab) == 'cd' .OR. TRIM(lab) == 'cd1'  ) THEN
            ALLOCATE(new_lab(2))
            new_lab(1)='cd'
            new_lab(2)='cd1'
            RETURN
         END IF
      END IF
      IF(TRIM(res) == 'hoh' .OR. TRIM(res) == 'tip'&
           & .OR. TRIM(res) == 'wat' .OR. TRIM(res) == 'spc') THEN
         IF(TRIM(lab) == 'o' .OR. TRIM(lab) == 'o1' .OR. TRIM(lab) == 'oh2'  ) THEN
            ALLOCATE(new_lab(3))
            new_lab(1)='o'
            new_lab(2)='o1'
            new_lab(2)='oh2'
            RETURN
         END IF
      END IF
      IF(TRIM(lab) == 'oct1' .OR. TRIM(lab) == 'ot1') THEN
         ALLOCATE(new_lab(2))
         new_lab(1)='oct1'
         new_lab(2)='ot1'
         RETURN
      END IF
         
      IF(TRIM(lab) == 'oct2' .OR. TRIM(lab) == 'ot2') THEN
         ALLOCATE(new_lab(2))
         new_lab(1)='oct2'
         new_lab(2)='ot2'
         RETURN
      END IF
      ALLOCATE(new_lab(1)); new_lab(1)=lab
    END SUBROUTINE AtomName_Exceptions
    FUNCTION ResName_Exceptions(lab_r, lab_t) RESULT(out)
      CHARACTER(len=*) :: lab_r,lab_t
      LOGICAL :: out
      INTEGER :: i,c_r,c_t

      out=.FALSE.
      c_r=0
      c_t=0
      DO i=1,SIZE(water)
         IF(TRIM(lab_r) == TRIM(water(i))) c_r=c_r+1
         IF(TRIM(lab_t) == TRIM(water(i))) c_t=c_t+1
      END DO
      IF(c_r /= 0 .AND. c_t /= 0) THEN
         out=.TRUE.
         RETURN
      END IF
      c_r=0
      c_t=0
      DO i=1,SIZE(histidine)
         IF(TRIM(lab_r) == TRIM(histidine(i))) c_r=c_r+1
         IF(TRIM(lab_t) == TRIM(histidine(i))) c_t=c_t+1
      END DO
      IF(c_r /= 0 .AND. c_t /= 0) THEN
         out=.TRUE.
         RETURN
      END IF

    END FUNCTION ResName_Exceptions
  END FUNCTION PDB__Validate
END MODULE PDB
