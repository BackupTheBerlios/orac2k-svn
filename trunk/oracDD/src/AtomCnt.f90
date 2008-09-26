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


  USE PI_
  USE IndPatch
  USE Constants
  USE Indsequence
  USE Errors,ONLY: Add_errors=>Add, Print_Errors, errmsg_f
  USE Myparse
  USE SecondarySeq
  USE Tops
  USE STRPAK
  USE Node
  USE Strings, ONLY: My_Fxm,My_Fam, MyRead
  USE Parameters
  USE Print_Defs
  IMPLICIT none
  PRIVATE
  PUBLIC :: AtomCnt_, AtomCnt__Type, AtomCnts, AtomCnt__Find, AtomCnt__Update&
       &, AtomCnt__Write, AtomCnt__Read

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
             jm=SIZE(App_Char(i_F)  % group (i) % g)
             IF(jm == 0) CYCLE
             Grp_No=Grp_No+1
             DO j=1,jm
                nato=nato+1
                AtomCnts(nato) % Res = App_Char(i_F) % Type
                nword=MyParse_(App_Char(i_F) % Type)
                IF(My_Fam('link', strngs(1))) &
                     &AtomCnts(nato) % Res = TRIM(ADJUSTL(strngs(2)))

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
                      CALL MyRead(App_Char(i_mass) % mass (2,o),AtomCnts(nato) % mass)
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
    INTEGER :: n,m,l,o,offset,nato,new_dim,old_dim,Res_Begins,Res_Ends&
         &,offset2,No_Res_Unit,offset_g,p_g


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
    offset2=0
    offset_g=0
    p_g=TempAtoms(Res_atm(2,Res_Ends)) % Grp_No-TempAtoms(Res_atm(1,Res_Begins)) % Grp_No+1

    No_res_Unit=Res_Ends-Res_Begins+1
    DO l=1,New_Units
       DO n=Res_Begins,Res_Ends
          DO m=Res_atm(1,n),Res_Atm(2,n)
             IF(ALLOCATED(TempAtoms(m) % cnt)) THEN
                o=SIZE(TempAtoms(m) % cnt)
                ALLOCATE(AtomCnts(m+offset) % cnt(o))
             END IF
             AtomCnts(m+offset) = TempAtoms(m)
             AtomCnts(m+offset) % Grp_no = TempAtoms(m) % Grp_No + offset_g
             AtomCnts(m+offset) % Res_no = TempAtoms(m) % Res_No + offset2
          END DO
          offset2=offset2+No_Res_Unit
       END DO
       offset_g=offset_g+p_g
       offset=offset+nato
    END DO
    DEALLOCATE(TempAtoms)
  END SUBROUTINE AtomCnt__Update
  SUBROUTINE AtomCnt__Write
    INTEGER :: n,m
    INTEGER :: ierr
    IF(PI_node == 0) THEN
       WRITE(kbinary) SIZE(AtomCnts)
       DO n=1,SIZE(AtomCnts)
          WRITE(kbinary) AtomCnts(n) % Res, AtomCnts(n) % beta, AtomCnts(n) %&
               & Betab, AtomCnts(n) % Res_no, AtomCnts(n) % Grp_No,&
               & AtomCnts(n) % Id_Res, AtomCnts(n) % Id_Type,&
               & AtomCnts(n) % Id_Slv, AtomCnts(n) % chg, AtomCnts(n) %&
               & mass
          IF(ALLOCATED(AtomCnts(n) % cnt)) THEN
             m=SIZE(AtomCnts(n) % cnt)
             WRITE(kbinary) m
             WRITE(kbinary) AtomCnts(n) % cnt
          ELSE
             m=0
             WRITE(kbinary) m
          END IF
       END DO
    END IF
#ifdef HAVE_MPI
    CALL MPI_Barrier(PI_Comm, ierr)
#endif
    WRITE(kprint,*) 'Writing Binary topology/parameter file for system ====>'
  END SUBROUTINE AtomCnt__Write
  FUNCTION AtomCnt__Read() RESULT(out)
    LOGICAL :: out
    INTEGER :: n,o,p

    WRITE(kprint,*) 'Read Binary topology/parameter file for system ====>'
    out=.TRUE.
    READ(kbinary,ERR=100,END=200) o
    IF(o == 0) RETURN
    ALLOCATE(AtomCnts(o))
    DO n=1,SIZE(AtomCnts)
       READ(kbinary,ERR=100,END=200) AtomCnts(n) % Res, AtomCnts(n) % beta&
            &, AtomCnts(n) % Betab, AtomCnts(n) % Res_no, AtomCnts(n) % Grp_No,&
            & AtomCnts(n) % Id_Res, AtomCnts(n) % Id_Type,&
            & AtomCnts(n) % Id_Slv, AtomCnts(n) % chg, AtomCnts(n) %&
            & mass
       READ(kbinary,ERR=100,END=200) p
       IF(p /= 0) THEN
          ALLOCATE(AtomCnts(n) % cnt (p))
          READ(kbinary,ERR=100,END=200) AtomCnts(n) % cnt
       END IF
    END DO
    RETURN
100 errmsg_f='Error while reading AtomCnt Parameters'
    CALL Add_Errors(-1,errmsg_f)
    out=.FALSE.
    RETURN
200 errmsg_f='End of file found while reading AtomCnt Parameters'
    CALL Add_Errors(-1,errmsg_f)
    out=.FALSE.
    RETURN    
  END FUNCTION AtomCnt__Read
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE AtomCnt
