!!$/---------------------------------------------------------------------\
!!$                                                                      |
!!$  Copyright (C) 2006-2007 Massimo Marchi <Massimo.Marchi@cea.fr>      |
!!$                                                                      |
!!$      This program is free software;  you  can  redistribute  it      |
!!$      and/or modify it under the terms of the GNU General Public      |
!!$      License version 2 as published  by  the  Free  Software         |
!!$      Foundation;                                                     |
!!$                                                                      |
!!$      This program is distributed in the hope that  it  will  be      |
!!$      useful, but WITHOUT ANY WARRANTY; without even the implied      |
!!$      warranty of MERCHANTABILITY or FITNESS  FOR  A  PARTICULAR      |
!!$      PURPOSE.   See  the  GNU  General  Public License for more      |
!!$      details.                                                        |
!!$                                                                      |
!!$      You should have received a copy of the GNU General  Public      |
!!$      License along with this program; if not, write to the Free      |
!!$      Software Foundation, Inc., 59  Temple  Place,  Suite  330,      |
!!$      Boston, MA  02111-1307  USA                                     |
!!$                                                                      |
!!$\---------------------------------------------------------------------/
FUNCTION PDB__Validate(Type, PDB__Coords) RESULT(out)
!!$***********************************************************************
!!$   Time-stamp: <2007-01-14 17:36:13 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Sun Jan 14 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program oracDD ----*

!!$-----------------------------------------------------------------------\
!!$                                                                       |
!!$--- Validate and write!!$/---------------------------------------------------------------------\
!!$                                                                      |
!!$  Copyright (C) 2006-2007 Massimo Marchi <Massimo.Marchi@cea.fr>      |
!!$                                                                      |
!!$      This program is free software;  you  can  redistribute  it      |
!!$      and/or modify it under the terms of the GNU General Public      |
!!$      License version 2 as published  by  the  Free  Software         |
!!$      Foundation;                                                     |
!!$                                                                      |
!!$      This program is distributed in the hope that  it  will  be      |
!!$      useful, but WITHOUT ANY WARRANTY; without even the implied      |
!!$      warranty of MERCHANTABILITY or FITNESS  FOR  A  PARTICULAR      |
!!$      PURPOSE.   See  the  GNU  General  Public License for more      |
!!$      details.                                                        |
!!$                                                                      |
!!$      You should have received a copy of the GNU General  Public      |
!!$      License along with this program; if not, write to the Free      |
!!$      Software Foundation, Inc., 59  Temple  Place,  Suite  330,      |
!!$      Boston, MA  02111-1307  USA                                     |
!!$                                                                      |
!!$\---------------------------------------------------------------------/

 on PDB__Coords                                  |
!!$                                                                       |
!!$-----------------------------------------------------------------------/

    CHARACTER(len=*) :: Type
    TYPE(AtomPDB), POINTER :: PDB__Coords(:)

    LOGICAL :: out
    LOGICAL, ALLOCATABLE :: oks(:)
    INTEGER :: n,Begins,Ends,m,p,ip_Res,oo,nword,offset,n_offset
    LOGICAL :: ok
    CHARACTER(len=max_char) :: lab_p,lab_m,res_n,lab_m0,res_m
    CHARACTER(len=max_char), POINTER  :: new_name(:)

    IF(.NOT. ASSOCIATED(PDB__Coords)) THEN
       errmsg_f='PDB__Init should be called before PDB__Validate'
       CALL Add_Errors(-1,errmsg_f)
       out=.FALSE.
       RETURN
    END IF
    out=.TRUE.

!!$
!!$--- Check if Tpg is initialized
!!$

    IF(.NOT. ALLOCATED(Tpg % atm)) THEN
       errmsg_f='Cannot validate PDB coordinates if system&
            & topology is not defined'
       CALL Add_Errors(-1,errmsg_f)
       out=.FALSE.
       RETURN
    END IF

!!$
!!$--- Check if the No. of residues correspond with those found
!!$--- in the .PDB file
!!$

    IF(Res_Ends-Res_Begins+1 /= SIZE(ResPdb)) THEN
       errmsg_f='No. of residues in '//TRIM(Type)//' .pdb file&
            & does not match system topology: Expected '&
            &//TRIM(MyPutnum(Res_Ends-Res_Begins+1))&
            &//' found '//TRIM(Myputnum(SIZE(ResPdb)))
       CALL Add_Errors(-1,errmsg_f)
       out=.FALSE.
       RETURN
    END IF

!!$
!!$--- oks are .true. only if all Tpg atoms are found on the .pdb with
!!$--- the exception of the hydrogens
!!$


!!$-- Counts total atoms
    ALLOCATE(oks(SIZE(Tpg % atm)))
    oks=.FALSE.

!!$-- Start the check on the secondary sequence
    offset=0
    IF(Res_Begins /= 1) THEN
       offset=Res_Atm(2,Res_Begins-1)-Res_Atm(1,1)+1
    END IF

    n_offset=Res_Begins-1
    DO n=Res_Begins,Res_Ends
       Begins=Res_Atm(1,n)
       Ends  =Res_Atm(2,n)
       res_n=ADJUSTL(Tpg % atm(Begins) % a % Res)

       nword=MyParse_(res_n)
       IF(TRIM(strngs(1)) == 'Link') THEN  !-- If the residue is linked
          res_n=ADJUSTL(strngs(2))          !-- take its real name
       END IF
!!$--- See if the ResNames match
       IF(.NOT. My_Fam(ResPdb(n-n_offset) % ResName, Res_n)) THEN 
          IF(.NOT. ResName_Exceptions(ResPdb(n-n_offset)&
               & % ResName, res_n)) THEN     ! There are exceptions...
               
             errmsg_w='System Topology residue type '&
                  &//TRIM(res_n)//' does not match &
                  &with .PDB file residue type '&
                  &//ResPdb(n-n_offset) % ResName//' run continues,&
                  & we hope this is stil ok'
             CALL Add_Errors(1,errmsg_w)
          END IF
       END IF

       DO m=1,SIZE(ResPdb(n-n_offset) % atm)
          lab_m0=ResPdb(n-n_offset) % atm(m) % AtmName
          res_m=ADJUSTL(ResPdb(n-n_offset) % ResName)
          CALL AtomName_Exceptions(lab_m0, res_m, new_name)             
          ok=.FALSE.
          DO oo=1,SIZE(new_name)
             lab_m=new_name(oo)
             DO p=Begins,Ends
                lab_p=Tpg % atm (p) % a % beta
                IF(TRIM(lab_m) == TRIM(Lab_p)) THEN
                   ok=.TRUE.
                   oks(p)=.TRUE.
                   PDB__Coords(p-offset) % x = ResPdb(n-n_offset) % atm(m) % x
                   PDB__Coords(p-offset) % y = ResPdb(n-n_offset) % atm(m) % y
                   PDB__Coords(p-offset) % z = ResPdb(n-n_offset) % atm(m) % z
                   PDB__Coords(p-offset) % AtmName = ResPdb(n-n_offset) % atm(m) % AtmName
                   PDB__Coords(p-offset) % Serial = p
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
       ok=.TRUE.
       DO m=Begins,Ends
          IF(Tpg % atm(m) % a % beta(1:1) ==  'h' .AND. oks(m)) ok=.FALSE.
       END DO
       DO m=Begins,Ends
          lab_m=Tpg % atm (m) % a % beta
          IF((.NOT. oks(m)) .AND. (Tpg % atm(m) % a % beta(1:1) /=  'h' )) THEN
             errmsg_f='While reading .PDB coordinates for residue No. '&
                  &//TRIM(Myputnum(n))//', cannot found label '//TRIM(lab_m)&
                  &//' of residue '//TRIM(res_n)
             CALL Add_Errors(-1,errmsg_f)
             out=.FALSE.
          END IF
       END DO
       CALL Print_Errors()
    END DO
    DO n=Res_Begins,Res_Ends
       Begins=Res_Atm(1,n)
       Ends  =Res_Atm(2,n)
       DO p=Begins,Ends
          IF(PDB__Coords(p-offset) % Serial == 0) PDB__Coords(p-offset) % Serial=p
       END DO
    END DO
    DEALLOCATE(oks,ResPdb)
  CONTAINS
    SUBROUTINE AtomName_Exceptions(lab, res, new_lab)
      CHARACTER(len=*) :: lab,res
      CHARACTER(len=max_char), POINTER :: new_lab(:)
      INTEGER :: n,m,o,p,len0
      INTEGER, SAVE :: first_Time=0
      IF(First_Time == 0) THEN
         First_Time=1
         ALLOCATE(ex(1) % res(1)); ALLOCATE(ex(1) % lab(2))
         ALLOCATE(ex(2) % res(4)); ALLOCATE(ex(2) % lab(3))
         ALLOCATE(ex(3) % res(1)); ALLOCATE(ex(3) % lab(2))
         ALLOCATE(ex(4) % res(1)); ALLOCATE(ex(4) % lab(2))
         ex=(/Name_Exception((/'ile'/),(/'cd ','cd1'/))&
              &, Name_Exception((/'hoh','tip','wat','spc'/),(/'o','o1','oh2'/))&
              &, Name_Exception((/' '/),(/'oct1','ot1'/))&
              &, Name_Exception((/' '/),(/'oct2','ot2'/)) /)
      END IF

      IF(ASSOCIATED(new_lab)) DEALLOCATE(new_lab)
      DO n=1,SIZE(ex)

!!$
!!$--- Check if the right ewsidue
!!$
         o=0
         DO m=1,SIZE(ex (n) % res)
            IF(TRIM(ex(n) % res(m)) == TRIM(res)) o=o+1
         END DO

!!$
!!$--- For a blank residue do always a check on the atomic labels
!!$
         
         IF(SIZE(ex(n) % res) == 1) THEN
            len0=LEN_TRIM(ex(n) % res(1))
            IF(len0 == 0) o=1
         END IF       
!!$
!!$--- Finally check on the atomic labels
!!$
         IF(o /= 0) THEN
            p=0
            DO m=1,SIZE(ex(n) % lab)
               IF(TRIM(ex(n) % lab(m)) == TRIM(lab)) p=p+1
            END DO
            IF(p == 0) CYCLE
            
            m=SIZE(ex(n) % lab)
            ALLOCATE(new_lab(m))
            new_lab=ADJUSTL(ex(n) % lab)
            RETURN
         END IF
      END DO
!!$
!!$--- All tests have failed use the original label only
!!$
      ALLOCATE(new_lab(1)); new_lab(1)=lab
    END SUBROUTINE AtomName_Exceptions
    FUNCTION ResName_Exceptions(lab_r, lab_t) RESULT(out)
      CHARACTER(len=*) :: lab_r,lab_t
      LOGICAL :: out
      INTEGER :: i,c_r,c_t
      CHARACTER(len=max_char), PARAMETER :: water(8)=(/'wat  ','hoh  ','tip3 '&
           &,'tip3p','tip3f','spc  ','hoh  ','spce '/)
      CHARACTER(len=max_char), PARAMETER :: histidine(5)=(/'his','hsd'&
           &,'hse','hsp','hs2'/)

      out=.FALSE.
!!$
!!$--- Try to see if the residue name corresponds to a common
!!$--- name of water
!!$
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
!!$
!!$--- Try to see if the residue name corresponds to a common 
!!$--- name of histidine
!!$
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
