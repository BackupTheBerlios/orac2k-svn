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
MODULE IndSequence

!!$***********************************************************************
!!$   Time-stamp: <2007-01-12 18:58:16 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Tue Jan  9 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*

  USE Parameters_Globals
  USE Errors, ONLY: Add_Errors=>Add, errmsg_f, Print_errors
  USE Constants
  USE SecondarySeq
  USE IndPatch
  USE Tops
  
  IMPLICIT none
  PRIVATE
  PUBLIC IndSequence_, IndSequence__type, IndSequence__Grp, IndSequence__Res&
       &, IndSequence__Pickres, Indsequence__SltSlv_Res, Indsequence__SltSlv_Grp
  TYPE IndSequence__Type
     INTEGER, DIMENSION(:), ALLOCATABLE :: i
  END TYPE IndSequence__Type
  TYPE(IndSequence__Type), DIMENSION(2), SAVE, TARGET :: Indexa
  INTEGER, ALLOCATABLE, SAVE, TARGET :: Res_Atm(:,:)
  INTEGER, ALLOCATABLE, SAVE, TARGET :: Grp_Atm(:,:)  
  INTEGER, SAVE, TARGET :: SltSlv_Res(2,2),SltSlv_Grp(2,2)

  LOGICAL, ALLOCATABLE, TARGET, SAVE :: ok_Residue(:)
CONTAINS
  FUNCTION IndSequence_() RESULT(out)
    TYPE(IndSequence__Type), DIMENSION(:), POINTER :: out
    INTEGER :: n,m,nato,Res_No,Grp_No,i_f,jm,j,i,i_L
    CHARACTER(len=Max_Char) :: res_i

    out=>NULL()
    DO n=1,SIZE(Secondary)
       ALLOCATE(Indexa(n) % i (SIZE(Secondary(n) % line)))
    END DO

    Res_No=0
    Grp_No=0

    DO n=1,SIZE(Secondary)
       sltslv_res(1,n)=Res_no+1
       sltslv_Grp(1,n)=Grp_no+1
       IF(.NOT. ALLOCATED(Secondary(n) % line)) CYCLE
       DO m=1,SIZE(Secondary(n) % line)
          Res_No=Res_No+1
          res_i=Secondary(n) % line(m)
          i_f=IndSequence__PickRes(res_i)
          Indexa(n) % i (m) = i_f

!!$
!!$--- Count atoms
!!$
          IF(Indexa(n) % i (m) == -1) THEN
             errmsg_f='Residue '//TRIM(res_i)//&
                  &' not found in force field database'
             CALL Add_Errors(-1,errmsg_f)
             RETURN
          END IF
          DO i=1,SIZE(App_Char(i_F) % group)
             Grp_No=Grp_No+1
             nato=nato+SIZE(App_Char(i_F)  % group (i) % g)
          END DO
       END DO
       sltslv_Res(2,n)=Res_No
       sltslv_Grp(2,n)=Grp_No
    END DO
    out=>Indexa
    
    ALLOCATE(Res_Atm(2,Res_No))
    ALLOCATE(Grp_Atm(2,Grp_No))

    nato=0
    Res_No=0
    Grp_No=0
    DO n=1,SIZE(Secondary)
       IF(.NOT. ALLOCATED(Secondary(n) % line)) CYCLE
       DO m=1,SIZE(Secondary(n) % line)

          Res_No=Res_No+1
          
          i_F=Indexa(n) % i (m) 
!!$
!!$--- Count atoms
!!$
          DO i=1,SIZE(App_Char(i_F) % group)
             Grp_No=Grp_No+1
             jm=SIZE(App_Char(i_F)  % group (i) % g)
             DO j=1,jm
                nato=nato+1
                IF(i == 1 .AND. j == 1) Res_Atm(1,Res_No)=nato
                IF(j == 1) Grp_Atm(1,Grp_No)=nato
             END DO
             Grp_Atm(2,Grp_No)=nato
          END DO
          Res_Atm(2,Res_No)=nato
       END DO
    END DO
  END FUNCTION IndSequence_
  FUNCTION IndSequence__Assign() RESULT(out)
    TYPE(IndSequence__Type), DIMENSION(:), POINTER :: out
    out=>NULL()
    IF((.NOT. ALLOCATED(Indexa(1) % i))&
         & .AND. (.NOT. ALLOCATED(Indexa(2) % i))) RETURN
    out=>Indexa
  END FUNCTION IndSequence__Assign
    
  FUNCTION IndSequence__Grp() RESULT (out)
    INTEGER, DIMENSION(:,:), POINTER :: out
    out=>NULL()
    IF(.NOT. ALLOCATED(Grp_Atm)) RETURN
    out=>Grp_Atm
  END FUNCTION IndSequence__Grp
    
  FUNCTION IndSequence__Res() RESULT (out)
    INTEGER, DIMENSION(:,:), POINTER :: out
    out=>NULL()
    IF(.NOT. ALLOCATED(Res_Atm)) RETURN
    out=>Res_Atm
  END FUNCTION IndSequence__Res

  FUNCTION IndSequence__SltSlv_Res() RESULT (out)
    INTEGER, DIMENSION(:,:), POINTER :: out
    out=>SltSlv_res
  END FUNCTION IndSequence__SltSlv_Res
  FUNCTION IndSequence__SltSlv_Grp() RESULT (out)
    INTEGER, DIMENSION(:,:), POINTER :: out
    out=>SltSlv_Grp
  END FUNCTION IndSequence__SltSlv_Grp
    
  FUNCTION IndSequence__PickRes(res) RESULT (out)
    CHARACTER(len=*) :: res
    INTEGER :: i,i_found,out
    i_found=-1
    DO i=1,SIZE(App_Char)
       IF(res  == App_Char(i)%Type) THEN
          i_found=i
          EXIT
       END IF
    END DO
    out=i_found
  END FUNCTION IndSequence__PickRes
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE IndSequence
