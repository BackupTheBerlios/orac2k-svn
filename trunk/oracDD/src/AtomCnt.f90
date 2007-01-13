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
  IMPLICIT none
  PRIVATE
  PUBLIC :: AtomCnt_, AtomCnt__Type, AtomCnts, AtomCnt__Find 

  TYPE AtomCnt__Type
     CHARACTER(len=max_atm) :: Res,beta, betab
     INTEGER :: Res_No, Grp_No, Id_Res, Id_Type, Id_slv
     REAL(8) :: chg,mass
     INTEGER, ALLOCATABLE :: cnt(:)
  END TYPE AtomCnt__Type
  TYPE(AtomCnt__Type), ALLOCATABLE, TARGET, SAVE :: AtomCnts(:)
  INTEGER, SAVE, POINTER :: Res_Atm(:,:)=>NULL()
  INTEGER, SAVE, POINTER :: Grp_Atm(:,:)=>NULL()
  INTEGER, SAVE, POINTER :: SltSlv(:,:)=>NULL()
CONTAINS
  SUBROUTINE AtomCnt_
    CHARACTER(len=max_char) :: res_i,line
    INTEGER :: nato, n, m, i, j, Res_No, ii, iflag, Grp_No,i_F,jm,nword
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

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE AtomCnt
