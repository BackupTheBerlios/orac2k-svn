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
MODULE AtomCnt

!!$***********************************************************************
!!$   Time-stamp: <2007-01-12 19:02:35 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Tue Nov 28 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program ORAC ----*


  USE IndPatch
  USE Constants
  USE Indsequence
  USE Errors,ONLY: Add_errors=>Add, Print_Errors, errmsg_f
  USE Myparse
  USE SecondarySeq
  USE Tops
  USE STRPAK
  USE Node
  USE Strings, ONLY: My_Fxm,My_Fam
  USE Parameters_Globals
  IMPLICIT none
  PRIVATE
  PUBLIC :: AtomCnt_, AtomCnt__Type, AtomCnts, AtomCnt__Find, AtomCnt__Update

  TYPE AtomCnt__Type
     CHARACTER(len=max_atm) :: Res,beta, betab
     INTEGER :: Res_No, Grp_No, Id_Res, Id_Type, Id_slv
     REAL(8) :: chg,mass
     INTEGER, ALLOCATABLE :: cnt(:)
  END TYPE AtomCnt__Type
  TYPE(AtomCnt__Type), ALLOCATABLE, TARGET, SAVE :: AtomCnts(:)
  INTEGER, SAVE, POINTER :: Res_Atm(:,:)=>NULL()
  INTEGER, SAVE, POINTER :: Grp_Atm(:,:)=>NULL()
CONTAINS
  SUBROUTINE AtomCnt_
    CHARACTER(len=max_char) :: res_i,line
    INTEGER :: nato, n, m, i, j, Res_No, ii, iflag, Grp_No,i_F,jm,nword,o,i_mass
    LOGICAL :: ok
    LOGICAL, SAVE :: Called=.FALSE.
    TYPE(Indpatch__type), DIMENSION(:), POINTER :: IndPa
    TYPE(Indsequence__type), DIMENSION(:), POINTER :: Inds

    indPa=>IndPatch_()
    IF(.NOT. ASSOCIATED(indPa)) CALL Print_Errors()
    inds=>IndSequence_()
    IF(.NOT. ASSOCIATED(inds)) CALL Print_Errors()
    Grp_Atm=>IndSequence__Grp()
    Res_Atm=>IndSequence__Res()

    nato=0
    DO n=1,SIZE(Secondary)
       IF(.NOT. ALLOCATED(Secondary(n) % line)) CYCLE
       DO m=1,SIZE(Secondary(n) % line)
          i_F=inds(n) % i (m) 
          DO i=1,SIZE(App_Char(i_F) % group)
             nato=nato+SIZE(App_Char(i_F)  % group (i) % g)
          END DO
       END DO
    END DO
    ALLOCATE(AtomCnts(nato))

!!$
!!$--- Count atoms
!!$
    i_mass=-1
    DO n=1,SIZE(App_Char)
       IF(My_Fxm('mass',App_Char(n) % Type)) THEN
          i_mass=n
          EXIT
       END IF
    END DO

    IF(i_mass == -1) THEN
       errmsg_f='Atomic masses not found on Topology file: List of &
            &masses in CHARMM format must be added at the&
            & beginning of any Topology file'
       CALL Add_Errors(-1,errmsg_f)
       CALL Print_Errors()
    END IF

    Res_No=0
    Grp_No=0
    nato=0
    DO n=1,SIZE(Secondary)
       IF(.NOT. ALLOCATED(Secondary(n) % line)) CYCLE
       DO m=1,SIZE(Secondary(n) % line)
          Res_No=Res_No+1
          i_F=inds(n) % i (m)
          DO i=1,SIZE(App_Char(i_F) % group)
             Grp_No=Grp_No+1
             jm=SIZE(App_Char(i_F)  % group (i) % g)
             DO j=1,jm
                nato=nato+1
                AtomCnts(nato) % Res = App_Char(i_F) % Type
                AtomCnts(nato) % Res_No = Res_No
                AtomCnts(nato) % Id_Res = i_F
                AtomCnts(nato) % Grp_No = Grp_No
                AtomCnts(nato) % Id_Slv = n
                nword=Myparse_(App_Char(i_F)  % group (i) % g(j))
                AtomCnts(nato) % beta = TRIM(strngs(1))
                AtomCnts(nato) % betab = TRIM(strngs(2))
                CALL SP_Getnum(strngs(3),AtomCnts (nato) % chg, iflag)
                ok=.FALSE.
                DO o=1,SIZE(App_Char(i_mass) % mass,2)
                   IF(My_Fxm(TRIM(AtomCnts(nato) % betab), App_Char(i_mass) % mass (1,o))) THEN
                      AtomCnts(nato) % Id_Type = o
                      ok=.TRUE.
                      EXIT
                   END IF
                END DO
                IF(.NOT. ok) THEN
                   errmsg_f='Atom type '//TRIM(AtomCnts(nato) % betab)&
                        &//' not found in the parameter list'
                   CALL Add_Errors(-1,errmsg_f)
                   CALL Print_Errors()
                END IF
             END DO
          END DO
       END DO
    END DO

    CALL AtomCnt__GetConnections
  CONTAINS
    INCLUDE 'AtomCnt__GetConnections.f90'
  END SUBROUTINE AtomCnt_
  FUNCTION AtomCnt__Find(No,label) RESULT(out)
    INTEGER :: No, out
    CHARACTER(len=*) :: label
    INTEGER :: n, i_Z
    CHARACTER(len=max_char) :: label0

    IF((.NOT. ASSOCIATED(Res_Atm)) .OR. (.NOT. ALLOCATED(App_Char))) THEN
       errmsg_f='Cannot find atom without initialization'
       CALL Add_Errors(-1,errmsg_f)
       CALL Print_Errors()
    END IF

    i_Z=AtomCnts(Res_Atm(1,No)) % Id_Res
    out=-1
    DO n=Res_Atm(1,No),Res_Atm(2,No)
       IF(TRIM(AtomCnts(n) % beta) == TRIM(label)) THEN
          out=n
          RETURN
       END IF
    END DO
    WRITE(label0,'(1x,i4,1x)') No
    errmsg_f='Inter-residue connection of secondary sequence not found. Atom '&
         &//TRIM(label)//' not on residue No. '//TRIM(label0)&
         &//'[ of type '//TRIM(App_Char(i_Z) % Type)//'] &
         &of the secondary sequence.'
    CALL Add_Errors(-1,errmsg_f)
    CALL Print_Errors()
  END FUNCTION AtomCnt__Find

  SUBROUTINE AtomCnt__Update(New_Units)
    INTEGER :: New_Units
    INTEGER, POINTER :: SltSlv(:,:)=>NULL()
    TYPE(AtomCnt__Type), POINTER :: TempAtoms(:)
    INTEGER :: n,m,l,o,offset,nato,new_dim,old_dim,Res_Begins,Res_Ends


    SltSlv=>IndSequence__SltSlv_Res()
    Res_Begins=sltSlv(1,2)
    Res_Ends=sltSlv(2,2)
    nato=Res_Atm(2,Res_Ends)-Res_Atm(1,Res_Begins)+1
    New_Dim=(New_Units-1)*nato
    New_Dim=New_Dim+SIZE(AtomCnts)

    old_Dim=SIZE(AtomCnts)
    ALLOCATE(TempAtoms(old_Dim))
    DO n=1,SIZE(AtomCnts)
       IF(ALLOCATED(AtomCnts (n) % cnt)) THEN
          ALLOCATE(TempAtoms (n) % cnt(SIZE(AtomCnts(n) % cnt)))
       END IF
    END DO

    TempAtoms=AtomCnts
    DEALLOCATE(AtomCnts)
    ALLOCATE(AtomCnts(New_Dim))
    DO n=1,old_dim-nato
       IF(ALLOCATED(TempAtoms (n) % cnt)) THEN
          o=SIZE(TempAtoms(n) % cnt)
          ALLOCATE(AtomCnts (n) % cnt(o))
       END IF
    END DO

    AtomCnts(1:old_Dim-nato)=TempAtoms(1:old_Dim-nato)
    offset=0

    DO l=1,New_Units
       DO n=Res_Begins,Res_Ends
          DO m=Res_atm(1,n),Res_Atm(2,n)
             IF(ALLOCATED(TempAtoms(m) % cnt)) THEN
                o=SIZE(TempAtoms(m) % cnt)
                ALLOCATE(AtomCnts(m+offset) % cnt(o))
             END IF
             AtomCnts(m+offset) = TempAtoms(m)
          END DO
       END DO
       offset=offset+nato
    END DO

    DEALLOCATE(TempAtoms)
  END SUBROUTINE AtomCnt__Update
  SUBROUTINE AtomCnt__Write
    INTEGER :: n
    WRITE(kbin) SIZE(AtomCnts)
    DO n=1,SIZE(AtomCnts)
       WRITE(kbin) AtomCnts(n) % Res, AtomCnts(n) % beta, AtomCnts(n) %&
            & Betab, AtomCnts(n) % Res_no, AtomCnts(n) % Grp_No,&
            & AtomCnts(n) % Id_Res, AtomCnts(n) % Id_Type,&
            & AtomCnts(n) % Id_Slv, AtomCnts(n) % chg, AtomCnts(n) %&
            & mass, SIZE(AtomCnts(n) % cnt)
       WRITE(kbin) AtomCnts(n) % cnt
    END DO
  END SUBROUTINE AtomCnt__Write
  SUBROUTINE AtomCnt__Read
    INTEGER :: n,o,p
    WRITE(kbin) o
    ALLOCATE(AtomCnts(o))
    DO n=1,SIZE(AtomCnts)
       READ(kbin) AtomCnts(n) % Res, AtomCnts(n) % beta, AtomCnts(n) %&
            & Betab, AtomCnts(n) % Res_no, AtomCnts(n) % Grp_No,&
            & AtomCnts(n) % Id_Res, AtomCnts(n) % Id_Type,&
            & AtomCnts(n) % Id_Slv, AtomCnts(n) % chg, AtomCnts(n) %&
            & mass, p
       ALLOCATE(AtomCnts(n) % cnt (p))
       READ(kbin) AtomCnts(n) % cnt(1:p)
    END DO
  END SUBROUTINE AtomCnt__Write
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE AtomCnt
