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
MODULE Parameters

!!$***********************************************************************
!!$   Time-stamp: <2007-01-09 10:56:45 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Nov 23 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program ORAC ----*

!!$---- DATA Only Modules -----------------------------------------------*

  USE SecondarySeq
  USE Constants, ONLY: max_pars,max_data, max_char
  USE Parameters_globals

!!$---- Modules ---------------------------------------------------------*

  USE Tree
  USE Errors, ONLY: Add_Errors=>Add, error_other, error_unr, error_args&
       &, errmsg_f, errmsg_w
  USE Strings, ONLY: MY_Fxm
  USE Myparse 
  USE STRPAK
  USE Resid

!!$---- DATA Statements -------------------------------------------------*

  IMPLICIT none
  PRIVATE
  PUBLIC Parameters__Scan, kbinary, Called_Binary, Called_Tpg, Called_Prm
  CHARACTER(len=max_char), SAVE :: input
  INTEGER, SAVE :: kbinary=0
  LOGICAL, SAVE :: Called_Join=.FALSE., Called_Patch=.FALSE., Called_Tpg=.FALSE.&
       &, Called_Prm=.FALSE., Called_Binary=.FALSE.
  CHARACTER(len=max_char), SAVE :: ftpg_read,fpar_read,fbinary='TOPPARAM.bin'
CONTAINS

!!$---- EXTECUTABLE Statements ------------------------------------------*

  SUBROUTINE Parameters__Scan
    CHARACTER(len=max_pars) :: line,linea
    INTEGER :: n,nword
    TYPE(Branch), SAVE :: check
    CALL Tree__Check_Tree('&PARAM',check)
    IF(.NOT. ASSOCIATED(check%children)) RETURN

    DO n=1,SIZE(check%children)
       line=TRIM(check%children(n))
       nword=MYParse_(line)
       
       linea=strngs(1)

       IF(MY_Fxm('READ_TPG',linea)) THEN
          Called_Tpg=.TRUE.
          CALL Resid_(strngs(2))
       ELSE IF(MY_Fxm('READ_PRM',linea)) THEN
          Called_Prm=.TRUE.
          CALL Resid_(strngs(2))
       ELSE IF(MY_Fxm('JOIN',linea)) THEN
          Called_Join=.TRUE.
          CALL Join(TRIM(line))
       ELSE IF(MY_Fxm('PATCH',linea)) THEN
          Called_Patch=.TRUE.
          CALL Patch_it(TRIM(line))
       ELSE IF(MY_Fxm('BINAR',linea) ) THEN
          Called_Binary=.TRUE.
          CALL Binary
       ELSE
          errmsg_f='Illegal commmands found:'//TRIM(linea)
          CALL Add_Errors(-1,errmsg_f)
       END IF
    END DO
    CALL Validate
    CONTAINS
      SUBROUTINE Validate        
        IF(Called_Tpg .NEQV. Called_Prm) THEN
           errmsg_f='Both Topology and Parameter files are needed. &
                &It seems that one is read but the other is not'
           CALL Add_Errors(-1,errmsg_f)
        ELSE IF((.NOT. Called_Tpg) .AND. (.NOT. Called_Prm) ) THEN
           IF(.NOT. Called_Binary) THEN
              errmsg_f='System topology and parameters are not avalaible&
                   & for this run: Need a binary Topology and Parameter file'
              CALL Add_Errors(-1,errmsg_f)
           END IF
        ELSE IF(Called_Tpg .AND. Called_Prm) THEN
           IF(.NOT. Called_Binary) THEN
              errmsg_w='You are not saving the system topology and parameters&
                   & on a binary file. Do you really want this?'
              CALL Add_Errors(1,errmsg_w)
           END IF
           IF(.NOT. Called_join) THEN
              errmsg_f='The command ''JOIN'' has not been used: &
                   &Need to read a sequence of residue to build&
                   & system'
              CALL Add_Errors(-1,errmsg_f)
           END IF
        END IF
      END SUBROUTINE Validate
  END SUBROUTINE Parameters__Scan
  SUBROUTINE Binary
    INTEGER ::  io,nword

    nword=SIZE(strngs)

    SELECT CASE(nword)
    CASE(1)
       CALL CHANNEL(io)
       kbinary=io
       OPEN(unit=kbinary,file=fbinary,form='UNFORMATTED',status='UNKNOWN')
    CASE(2) 
       CALL CHANNEL(io)
       kbinary=io
       fbinary=TRIM(strngs(2))
       OPEN(unit=kbinary,file=fbinary,form='UNFORMATTED',status='UNKNOWN')
    CASE DEFAULT
       errmsg_f=error_args % g (2)//' 0 or 1 '
       CALL Add_Errors(-1,errmsg_f)
       RETURN
    END SELECT
  END SUBROUTINE Binary
  SUBROUTINE Join(name)
    CHARACTER(len=*) :: name
    TYPE(Branch), SAVE :: checks
    CHARACTER(len=max_pars) :: line
    INTEGER :: c,n,m,rept,p,iflag,nword
    INTEGER, SAVE :: count_J=0
    CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE  :: Seq
    CHARACTER(len=max_char)  :: Type_Seq

    CALL Tree__Check_Tree(name,checks)
    IF(.NOT. ASSOCIATED(checks%children)) RETURN

    
    line=TRIM(checks%name)
    nword=MYParse_(line)
    IF(MY_Fxm('SOLUTE',strngs(2))) THEN
       count_J=1
    ELSE IF(MY_Fxm('SOLVENT',strngs(2))) THEN
       count_J=2
    ELSE 
       errmsg_f='JOIN must be followed by either SOLUTE&
            & or SOLVENT. Found '//TRIM(strngs(2))//' instead'
       CALL Add_Errors(-1,errmsg_f)
       RETURN
    END IF
    Type_Seq=strngs(2)

    m=SIZE(checks%children)
    ALLOCATE(Seq (m))
    c=0
    DO n=m,1,-1
       Seq(m-n+1)=TRIM(checks%children(n))
       line=TRIM(checks%children(n))
       nword=MYParse_(line) 
       p=0
       DO WHILE(p < nword)
          p=p+1
          IF(strngs(p) /= 'x') THEN
             c=c+1
          ELSE
             CALL SP_Getnum(strngs(p+1), rept,iflag)
             c=c+rept
             p=p+1
          END IF
       END DO
    END DO
    CALL SecondarySeq_(count_J,Type_Seq,Seq)
    DEALLOCATE(Seq)
  END SUBROUTINE Join
  SUBROUTINE Patch_it(name)
    CHARACTER(len=*) :: name
    TYPE(Branch), SAVE :: checks
    CHARACTER(len=max_pars) :: line,linea
    INTEGER :: n,nword,iflag

    CALL Tree__Check_Tree(name,checks)
    IF(.NOT. ASSOCIATED(checks%children)) RETURN
    ALLOCATE(patches(SIZE(checks%children)))
    DO n=1,SIZE(checks%children)
       line=TRIM(checks%children(n))
       nword=MYParse_(line)
       linea=strngs(1)
       IF(MY_Fxm('resi',linea)) THEN
          IF(nword /= 4) THEN
             errmsg_f=error_args % g (4)//' 4'
             CALL Add_Errors(-1,errmsg_f)
             EXIT
          END IF
          patches(n)%Type='resi'
          patches(n)%New_Res=TRIM(strngs(2))
          patches(n)%pres=TRIM(strngs(3))
          patches(n)%res=TRIM(strngs(4))
       ELSE IF(MY_Fxm('link',linea)) THEN
          IF(nword /= 6) THEN
             errmsg_f=error_args % g (4)//' 6'
             CALL Add_Errors(-1,errmsg_f)
             EXIT
          END IF
          patches(n)%Type='link'
          patches(n)%pres=TRIM(strngs(2))
          patches(n)%Res_l(1)=TRIM(strngs(3))
          patches(n)%Res_l(2)=TRIM(strngs(5))
          CALL SP_GETNUM(strngs(4),patches(n)%one,iflag)
          IF(iflag /= 0) THEN
             errmsg_f='Reading No. of Residue Failed'
             CALL Add_Errors(-1,errmsg_f)
          END IF
          CALL SP_GETNUM(strngs(6),patches(n)%two,iflag)
          IF(iflag /= 0) THEN
             errmsg_f='Reading No. of Residue Failed'
             CALL Add_Errors(-1,errmsg_f)
          END IF
       ELSE
          errmsg_f=error_unr % g (2)//TRIM(linea)
          CALL Add_Errors(-1,errmsg_f)
          EXIT
       END IF
    END DO

  END SUBROUTINE Patch_it
END MODULE Parameters
