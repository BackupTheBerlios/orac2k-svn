MODULE Cell
  USE Constants, ONLY: max_pars,max_data,max_char
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, error_args, errmsg_f
  USE STRPAK, ONLY: SP_Getnum
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Cell_, co,oc,a,b,c,alpha,beta,gamma
  REAL(8), DIMENSION(:,:), ALLOCATABLE, SAVE :: co,oc
  REAL(8), SAVE :: a, alpha=90.0D0
  REAL(8), SAVE :: b, c, beta,gamma
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
    CASE(7)
       CALL SP_Getnum(strngs(2),a,iflags)
       CALL SP_Getnum(strngs(3),b,iflags)
       CALL SP_Getnum(strngs(4),c,iflags)
       CALL SP_Getnum(strngs(5),alpha,iflags)
       CALL SP_Getnum(strngs(6),beta,iflags)
       CALL SP_Getnum(strngs(7),gamma,iflags)
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
    CALL MatInv(co,oc)
    
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
END MODULE Cell
MODULE Solute
  USE Constants, ONLY: max_pars,max_data, max_char
  USE PDB
  USE Tree
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, error_args, errmsg_f
  USE Strings, ONLY: MY_Fxm
  USE Myparse
  IMPLICIT NONE 
  PRIVATE
  PUBLIC :: Solute_,PDB_Solute, PDB_Template
  CHARACTER(len=max_char), DIMENSION(:), POINTER, SAVE :: PDB_Solute, PDB_Template
CONTAINS
  SUBROUTINE Solute_(name)
    CHARACTER(len=*) :: name
    TYPE(Branch), SAVE :: check
    CHARACTER(len=max_pars) :: line,linea
    CHARACTER(len=max_char) :: lab0
    INTEGER :: n,nword

    CALL Tree__Check_Tree(name,check)
    IF(.NOT. ASSOCIATED(check%children)) RETURN

    DO n=1,SIZE(check%children)
       line=TRIM(check%children(n))
       nword=MYParse_(line)
       linea=strngs(1)
       
       IF(MY_Fxm('coord',linea)) THEN
          IF(nword /= 2) THEN
             WRITE(lab0,'(i2)') nword
             errmsg_f=error_args % g (4) //' 1 whereas it was '//TRIM(lab0)//' : '//TRIM(line)
             CALL Add_Errors(-1,errmsg_f)
             RETURN
          END IF
          WRITE(*,*) 'Storing SOLUTE .pdb ====>'
          CALL PDB_('Solute',strngs(2),PDB_Solute)
       ELSE IF(MY_Fxm('temp',linea)) THEN
          IF(nword /= 2) THEN
             WRITE(lab0,'(i2)') nword
             errmsg_f=error_args % g (4) //' 1 whereas it was '//TRIM(lab0)//' : '//TRIM(line)
             CALL Add_Errors(-1,errmsg_f)
             RETURN
          END IF
          WRITE(*,*) 'Storing TEMPLATE .pdb ====>'
          CALL PDB_('Template',strngs(2),PDB_Template)
       ELSE IF(MY_Fxm('scale',linea)) THEN
          CONTINUE
       ELSE
          errmsg_f='Illegal commmands found:'//TRIM(linea)
          CALL Add_Errors(-1,errmsg_f)
       END IF
    END DO
    CALL Print_Errors()
  END SUBROUTINE Solute_
END MODULE Solute
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

!!$---- DATA Statements -------------------------------------------------*

  IMPLICIT none
  PRIVATE
  PUBLIC Setups__Scan
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

!!$       ELSE IF(MY_Fxm('RESET',linea)) THEN
!!$          CALL Reset_CM
!!$       ELSE IF(MY_Fxm('SOLV',linea)) THEN
!!$          CALL Solvent(TRIM(line))
!!$       ELSE IF(MY_Fxm('TEMPLATE',linea)) THEN
!!$          CALL Template
       ELSE
          errmsg_f='Illegal commmands found:'//TRIM(linea)
          CALL Add_Errors(-1,errmsg_f)
       END IF
    END DO
  END SUBROUTINE Setups__Scan
END MODULE Setup
