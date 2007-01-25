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
MODULE Cell
  USE Constants, ONLY: max_pars,max_data,max_char
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, error_args, errmsg_f
  USE STRPAK, ONLY: SP_Getnum
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Cell_, Cell__Volume, co,oc,a,b,c,alpha,beta,gamma, volume
  REAL(8), DIMENSION(:,:), ALLOCATABLE, SAVE :: co,oc
  REAL(8), SAVE :: a, alpha=90.0D0
  REAL(8), SAVE :: b, c, beta,gamma,Volume
CONTAINS
  SUBROUTINE Cell_(strngs)
    USE Numerics
    USE Units
    CHARACTER(len=max_pars), DIMENSION(:) :: strngs
    INTEGER ::nword,iflags
    nword=SIZE(strngs)
    SELECT CASE(nword)
    CASE(2)
       CALL SP_Getnum(strngs(2),a,iflags)
       b=a
       c=a
       beta=alpha
       gamma=alpha
    CASE(4)
       CALL SP_Getnum(strngs(2),a,iflags)
       CALL SP_Getnum(strngs(3),b,iflags)
       CALL SP_Getnum(strngs(4),c,iflags)
       beta=alpha
       gamma=alpha
       IF(iflags /= 0) THEN
          errmsg_f='Internal reading error: Module Cell'
          CALL Add_Errors(-1,errmsg_f)
          RETURN
       END IF
    CASE(7)
       CALL SP_Getnum(strngs(2),a,iflags)
       CALL SP_Getnum(strngs(3),b,iflags)
       CALL SP_Getnum(strngs(4),c,iflags)
       CALL SP_Getnum(strngs(5),alpha,iflags)
       CALL SP_Getnum(strngs(6),beta,iflags)
       CALL SP_Getnum(strngs(7),gamma,iflags)
       IF(iflags /= 0) THEN
          errmsg_f='Internal reading error: Module Cell'
          CALL Add_Errors(-1,errmsg_f)
          RETURN
       END IF
    CASE DEFAULT 
       errmsg_f=error_args % g (4)//' 1, 3 or 6'
       CALL Add_Errors(-1,errmsg_f)
       RETURN
    END SELECT
    IF(iflags /= 0) THEN
       errmsg_f='Internal reading error: Module Setup_Cell'
       CALL Add_Errors(-1,errmsg_f)
       RETURN
    END IF
    ALLOCATE(co(3,3),oc(3,3))
    CALL GetCO
    Volume=MatInv(co,oc)*boxl**3
    
    IF(Determinant == 0.0D0) THEN
       errmsg_f='Cell parameters are probably wrong: CO matrix is singular'
       CALL Add_Errors(-1,errmsg_f)
       CALL Print_Errors()
    END IF
  CONTAINS
    SUBROUTINE GetCO
      REAL(8) ::  degrad,qt,alf,bet,gam,ax,bx,by,cx,cy,cz
      INTEGER :: n
      
      co(2,1)=0.0D0
      co(3,1)=0.0D0
      co(3,2)=0.0D0
      degrad=pi/180.0d0
      ax=a
      alf=DCOS(degrad*alpha)
      bet=DCOS(degrad*beta)
      qt=DSIN(degrad*gamma)
      gam=DCOS(degrad*gamma)
      bx=b*gam
      by=b*qt
      cx=c*bet
      cy=c*(alf-bet*gam)/qt
      cz=dsqrt(c*c-cx*cx-cy*cy)
      co(1,1)=ax
      co(1,2)=bx
      co(1,3)=cx
      co(2,2)=by
      co(2,3)=cy
      co(3,3)=cz
      DO n=1,3
         co(n,1)=co(n,1)/boxl
         co(n,2)=co(n,2)/boxl
         co(n,3)=co(n,3)/boxl
      END DO
      WHERE(ABS(co) < 1.0D-8)
         co=0.0D0
      END WHERE
    END SUBROUTINE GetCO
  END SUBROUTINE Cell_
  FUNCTION Cell__Volume(a0,b0,c0,alf0,bet0,gam0) RESULT(out)
    USE Units
    USE Numerics
    REAL(8) :: a0,b0,c0
    REAL(8) :: alf0,bet0,gam0
    REAL(8) :: out
    REAL(8) ::  degrad,qt,alf,bet,gam,ax,bx,by,cx,cy,cz
    REAL(8) :: co0(3,3),oc0(3,3),volume0
    INTEGER :: n

    co0(2,1)=0.0D0
    co0(3,1)=0.0D0
    co0(3,2)=0.0D0
    degrad=pi/180.0d0
    ax=a0
    alf=DCOS(degrad*alf0)
    bet=DCOS(degrad*bet0)
    qt=DSIN(degrad*gam0)
    gam=DCOS(degrad*gam0)
    bx=b0*gam
    by=b0*qt
    cx=c0*bet
    cy=c0*(alf-bet*gam)/qt
    cz=dsqrt(c0*c0-cx*cx-cy*cy)
    co0(1,1)=ax
    co0(1,2)=bx
    co0(1,3)=cx
    co0(2,2)=by
    co0(2,3)=cy
    co0(3,3)=cz
    DO n=1,3
       co0(n,1)=co0(n,1)/boxl
       co0(n,2)=co0(n,2)/boxl
       co0(n,3)=co0(n,3)/boxl
    END DO
    WHERE(ABS(co0) < 1.0D-8)
       co0=0.0D0
    END WHERE
    Volume0=MatInv(co0,oc0)*boxl**3
    out=Volume0
  END FUNCTION Cell__Volume
END MODULE Cell
MODULE Solute
  USE Constants, ONLY: max_pars,max_data, max_char
  USE ReadStore
  USE Tree
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, error_args, errmsg_f
  USE Strings, ONLY: MY_Fxm
  USE Myparse
  USE STRPAK, ONLY: SP_Getnum
  IMPLICIT NONE 
  PRIVATE
  PUBLIC :: Solute_,PDB_Solute, PDB_Template, Solute__Exclusion
  CHARACTER(len=max_char), ALLOCATABLE, SAVE :: PDB_Solute(:), PDB_Template(:)
  REAL(8), SAVE :: Solute__Exclusion=1.0D0
CONTAINS
  SUBROUTINE Solute_(name)
    CHARACTER(len=*) :: name
    TYPE(Branch), SAVE :: check
    CHARACTER(len=max_pars) :: line,linea
    CHARACTER(len=max_char) :: lab0
    INTEGER :: n,nword,iflags

    CALL Tree__Check_Tree(name,check)
    IF(.NOT. ASSOCIATED(check%children)) RETURN

    DO n=1,SIZE(check%children)
       line=TRIM(check%children(n))
       nword=MYParse_(line)
       linea=strngs(1)
       
       IF(MY_Fxm('coord',linea)) THEN
          IF(nword /= 2) THEN
             WRITE(lab0,'(i2)') nword
             errmsg_f=error_args % g (4) //' 1 whereas it was '&
                  &//TRIM(lab0)//' : '//TRIM(line)
             CALL Add_Errors(-1,errmsg_f)
             RETURN
          END IF

          WRITE(*,*) 'Storing SOLUTE .pdb ====>'
          IF(.NOT. ReadStore_(strngs(2))) CALL Print_Errors()
          ALLOCATE(PDB_Solute(SIZE(RS__string)))
          PDB_Solute=RS__string
          CALL ReadStore__Delete
       ELSE IF(MY_Fxm('temp',linea)) THEN
          IF(nword /= 2) THEN
             WRITE(lab0,'(i2)') nword
             errmsg_f=error_args % g (4) //' 1 whereas it was '&
                  &//TRIM(lab0)//' : '//TRIM(line)
             CALL Add_Errors(-1,errmsg_f)
             RETURN
          END IF
          WRITE(*,*) 'Storing TEMPLATE .pdb ====>'
          IF(.NOT. ReadStore_(strngs(2))) CALL Print_Errors()
          ALLOCATE(PDB_Template(SIZE(RS__string)))
          PDB_Template=RS__string
          CALL ReadStore__Delete
       ELSE IF(MY_Fxm('scale',linea)) THEN
          CONTINUE
       ELSE IF(MY_Fxm('exclu',linea)) THEN
          IF(nword /= 2) THEN
             WRITE(lab0,'(i2)') nword
             errmsg_f=error_args % g (4) //' 1 whereas it was '&
                  &//TRIM(lab0)//' : '//TRIM(line)
             CALL Add_Errors(-1,errmsg_f)
             RETURN
          END IF
          CALL SP_Getnum(strngs(2),Solute__exclusion,iflags)
          IF(iflags /= 0) THEN
             errmsg_f='Internal reading error: Module Solute'
             CALL Add_Errors(-1,errmsg_f)
             RETURN
          END IF
          
       ELSE
          errmsg_f='Illegal commmands found:'//TRIM(linea)
          CALL Add_Errors(-1,errmsg_f)
       END IF
    END DO
    CALL Print_Errors()
  END SUBROUTINE Solute_
END MODULE Solute
MODULE Solvent
  USE Constants, ONLY: max_pars,max_data, max_char
  USE ReadStore
  USE Tree
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, error_args, errmsg_f
  USE Strings, ONLY: MY_Fxm
  USE Myparse
  USE STRPAK
  IMPLICIT NONE 
  PRIVATE
  CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE, SAVE :: PDB_Solvent
  TYPE :: Solvent__Type
     LOGICAL :: Build=.FALSE.
     INTEGER :: Added=0
     INTEGER :: replicate(3)=(/0,0,0/)
     CHARACTER(len=max_Char) :: Cell_Type='SC'
     REAL(8) :: rho=-1.0D0
  END TYPE Solvent__Type
  TYPE :: Solvent__Cell
     SEQUENCE
     CHARACTER(len=3) :: type
     INTEGER :: nt
     REAL(8) :: t(3,4)
  END TYPE Solvent__Cell

  TYPE(Solvent__Type), SAVE :: Solvent__Param
  TYPE(Solvent__Cell), PARAMETER :: Solvent__Box(3)=(/ Solvent__Cell('SC ',1&
       &,RESHAPE((/ &
       & 0.0D0, 0.0D0, 0.0D0, 0.0D0,&
       & 0.0D0, 0.0D0, 0.0D0, 0.0D0,&
       & 0.0D0, 0.0D0, 0.0D0, 0.0D0  /),(/3,4/)) ), Solvent__Cell('BBC',2&
       &,RESHAPE((/ &
       & 0.0D0, 0.5D0, 0.0D0, 0.0D0,&
       & 0.0D0, 0.5D0, 0.0D0, 0.0D0,&
       & 0.0D0, 0.5D0, 0.0D0, 0.0D0 /),(/3,4/)) ), Solvent__Cell('FCC',4&
       &,RESHAPE((/ &
       & 0.0D0, 0.5D0, 0.0D0,-0.5D0,&
       & 0.0D0,-0.5D0, 0.5D0, 0.0D0,&
       & 0.0D0, 0.0D0,-0.5D0, 0.5D0 /),(/3,4/)) )/)

  PUBLIC :: Solvent_,PDB_Solvent, Solvent__Param, Solvent__Type&
       &, Solvent__Cell, Solvent__Box
CONTAINS
  SUBROUTINE Solvent_(name)
    CHARACTER(len=*) :: name
    TYPE(Branch), SAVE :: check
    CHARACTER(len=max_pars) :: line,linea
    CHARACTER(len=max_char) :: lab0
    INTEGER :: n,nword,iflags

    CALL Tree__Check_Tree(name,check)
    IF(.NOT. ASSOCIATED(check%children)) RETURN

    DO n=1,SIZE(check%children)
       line=TRIM(check%children(n))
       nword=MYParse_(line)
       linea=strngs(1)
       
       IF(MY_Fxm('coord',linea)) THEN
          IF(nword /= 2) THEN
             WRITE(lab0,'(i2)') nword
             errmsg_f=error_args % g (4) //' 1 whereas it was '&
                  &//TRIM(lab0)//' : '//TRIM(line)
             CALL Add_Errors(-1,errmsg_f)
             CYCLE
          END IF
          WRITE(*,*) 'Storing SOLVENT .pdb ====>'
          IF(.NOT. ReadStore_(strngs(2))) CALL Print_Errors()
          ALLOCATE(PDB_Solvent(SIZE(RS__string)))
          PDB_Solvent=RS__string
          CALL ReadStore__Delete

       ELSE IF(My_Fxm('repli',linea)) THEN
          Solvent__Param % Build=.TRUE.
          SELECT CASE(nword-1)
          CASE(2)
             Solvent__Param % Cell_Type = TRIM(ADJUSTL(strngs(2)))
             CALL SP_Getnum(strngs(3),Solvent__Param % replicate(1),iflags)
             IF(iflags /=0) THEN
                errmsg_f='Cannot convert to integer: '''&
                     &//TRIM(strngs(3))//''' '
                CALL Add_Errors(-1,errmsg_f)
                CYCLE
             END IF
             Solvent__Param % replicate(2)=Solvent__Param % replicate(1)             
             Solvent__Param % replicate(3)=Solvent__Param % replicate(1)             
          CASE(4)
             Solvent__Param % Cell_Type = TRIM(ADJUSTL(strngs(2)))
             CALL SP_Getnum(strngs(3),Solvent__Param % replicate(1),iflags)
             CALL SP_Getnum(strngs(4),Solvent__Param % replicate(2),iflags)
             CALL SP_Getnum(strngs(5),Solvent__Param % replicate(2),iflags)
             IF(iflags /=0) THEN
                errmsg_f='Cannot convert to integer: '''&
                     &//TRIM(strngs(2))//' '//TRIM(strngs(3))&
                     &//' '//TRIM(strngs(4))//''' '
                CALL Add_Errors(-1,errmsg_f)
                CYCLE
             END IF
          CASE DEFAULT 
             errmsg_f=error_args % g (4)//' 2 or 4'
             CALL Add_Errors(-1,errmsg_f)
             CYCLE
          END SELECT
       ELSE IF(My_Fxm('add',linea)) THEN
          Solvent__Param % Build=.TRUE.
          IF(nword /= 2) THEN
             errmsg_f=error_args % g (4)//' 1 '
             CALL Add_Errors(-1,errmsg_f)
          END IF
          CALL SP_Getnum(strngs(2),Solvent__Param % added,iflags)
          IF(iflags /=0) THEN
             errmsg_f='Cannot convert to integer: '''//TRIM(strngs(3))//''' '
             CALL Add_Errors(-1,errmsg_f)
             CYCLE
          END IF

       ELSE IF(My_Fxm('molvol',linea)) THEN
          Solvent__Param % Build=.TRUE.
          IF(nword /= 2) THEN
             errmsg_f=error_args % g (4)//' 1 '
             CALL Add_Errors(-1,errmsg_f)
             CYCLE
          END IF
          CALL SP_Getnum(strngs(2),Solvent__Param % rho,iflags)
          IF(iflags /=0) THEN
             errmsg_f='Cannot convert to double: '''&
                  &//TRIM(strngs(2))//''' '
             CALL Add_Errors(-1,errmsg_f)
             CYCLE
          END IF
          Solvent__Param % rho=1.0D0/Solvent__Param % rho
       ELSE
          errmsg_f='Illegal commmands found: '//TRIM(linea)
          CALL Add_Errors(-1,errmsg_f)
       END IF
    END DO
    CALL Print_Errors()
  END SUBROUTINE Solvent_
END MODULE Solvent
MODULE Setup

!!$***********************************************************************
!!$   Time-stamp: <2006-12-20 11:50:31 marchi>                           *
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

  USE Constants, ONLY: max_pars,max_data, max_char
  USE Parameters_globals

!!$---- Modules ---------------------------------------------------------*

  USE Tree
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, error_other, error_unr&
       &, error_args, errmsg_f
  USE Strings, ONLY: MY_Fxm
  USE Myparse
  USE Cell
  USE Solute
  USE Solvent
  USE ReadStore

!!$---- DATA Statements -------------------------------------------------*

  IMPLICIT none
  PRIVATE
  PUBLIC Setups__Scan, Setup__PDB
  CHARACTER(len=max_char), ALLOCATABLE, SAVE :: Setup__PDB(:)
CONTAINS

!!$---- EXTECUTABLE Statements ------------------------------------------*

  SUBROUTINE Setups__Scan
    CHARACTER(len=max_pars) :: line,linea
    REAL(8) :: a
    INTEGER :: n,m,iflag,nword
    TYPE(Branch), SAVE :: check

    CALL Tree__Check_Tree('&SETUP',check)
    IF(.NOT. ASSOCIATED(check%children)) RETURN

    DO n=1,SIZE(check%children)
       line=TRIM(check%children(n))
       nword=Myparse_(line)
       linea=strngs(1)

       IF(MY_Fxm('CELL',linea)) THEN
          CALL Cell_(strngs)
       ELSE IF(MY_Fxm('SOLU',linea)) THEN
          CALL Solute_(TRIM(line))
       ELSE IF(MY_Fxm('SOLV',linea)) THEN
          CALL Solvent_(TRIM(line))
       ELSE IF(My_Fxm('PDB',linea)) THEN 
          CALL Store_SystemPDB
!!$       ELSE IF(MY_Fxm('RESET',linea)) THEN
!!$          CALL Reset_CM
       ELSE
          errmsg_f='Illegal commmands found:'//TRIM(linea)
          CALL Add_Errors(-1,errmsg_f)
       END IF
    END DO
  CONTAINS
    SUBROUTINE Store_SystemPDB
      WRITE(*,*) 'Storing System .pdb ====>'
      IF(.NOT. ReadStore_(strngs(2))) CALL Print_Errors()
      ALLOCATE(Setup__PDB(SIZE(RS__string)))
      Setup__PDB=RS__string
      CALL ReadStore__Delete
    END SUBROUTINE Store_SystemPDB
  END SUBROUTINE Setups__Scan
END MODULE Setup
