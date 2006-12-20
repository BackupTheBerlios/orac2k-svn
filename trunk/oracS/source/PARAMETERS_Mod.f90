MODULE PARAMETERS_Mod

!!$***********************************************************************
!!$   Time-stamp: <2006-12-20 15:07:23 marchi>                           *
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

  USE CONSTANTS, ONLY: max_pars,max_data, max_char
  USE PARAMETERS_GLOBALS

!!$---- Modules ---------------------------------------------------------*

  USE CLASS_Tree
  USE ERROR_Mod, ONLY: Add_Errors=>Add, error_other, error_unr, error_args
  USE STRINGS_Mod, ONLY: MY_Fxm
  USE MYPARSE_Mod, my_parse=>parse
  USE TOPS_Mod, ONLY: Transform
  USE TOPPAR_STORE_Mod, ONLY: ST_TopPar=>Store_TopPar, ST_GetTop, ST_GetPar, &
       &counter,Create_Keywords, ST_Init=>Init

!!$---- DATA Statements -------------------------------------------------*

  IMPLICIT none
  PRIVATE
  PUBLIC Scan
  CHARACTER(len=max_pars), DIMENSION(:), ALLOCATABLE, PRIVATE ::&
       & strngs
  CHARACTER(len=max_char), SAVE :: input
  CHARACTER(len=max_data) :: errmsg_w,errmsg_f
CONTAINS

!!$---- EXTECUTABLE Statements ------------------------------------------*

  SUBROUTINE Scan
    CHARACTER(len=max_pars) :: line,linea
    REAL(8) :: a
    INTEGER :: n,m,iflag
    TYPE(Branch), SAVE :: check

    CALL Check_Tree('&PARAM',check)
    IF(.NOT. ASSOCIATED(check%children)) RETURN

    DO n=1,SIZE(check%children)
       line=TRIM(check%children(n))
       CALL MY_Parse(line,strngs)
       linea=strngs(1)

       IF(MY_Fxm('READ_TPG',linea)) THEN
          CALL Store_Topology(strngs(2))
          IF(ALLOCATED(Topology)) THEN
             CALL Transform(Topology)
          END IF
       ELSE IF(MY_Fxm('READ_PRM',linea)) THEN
          CALL Store_Parameters(strngs(2))
       ELSE IF(MY_Fxm('JOIN',linea)) THEN
          CALL Join(TRIM(line))
       ELSE IF(MY_Fxm('PATCH',linea)) THEN
          CALL Patch_it(TRIM(line))
       ELSE IF(MY_Fxm('BINAR',linea) .OR. MY_Fxm('WRITE',linea) ) THEN
          CALL Binary
       ELSE
          errmsg_f='Illegal commmands found:'//TRIM(linea)
          CALL Add_Errors(-1,errmsg_f)
       END IF
    END DO
  END SUBROUTINE Scan
  SUBROUTINE Binary
    INTEGER ::nword,io
    nword=SIZE(strngs)
    IF(nword /= 2) THEN
       errmsg_f=error_args % g (2)//' 2'
       CALL Add_Errors(-1,errmsg_f)
       RETURN
    ELSE
       CALL CHANNEL(io)
       kbin=io
       OPEN(unit=io,file=strngs(2),form='UNFORMATTED',status='UNKNOWN')
    END IF
  END SUBROUTINE Binary
  SUBROUTINE Join(name)
    CHARACTER(len=*) :: name
    TYPE(Branch), SAVE :: checks
    CHARACTER(len=max_pars) :: line
    INTEGER :: c,n,o,m,rept,p,iflag
    INTEGER, SAVE :: count_J=0
    CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE  :: Seq

    CALL Check_Tree(name,checks)
    IF(.NOT. ASSOCIATED(checks%children)) RETURN

    
    line=TRIM(checks%name)
    CALL MY_Parse(line,strngs)
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
    Secondary_Seq(count_J) % type=strngs(2)

    m=SIZE(checks%children)
    ALLOCATE(Seq (m))
    c=0
    DO n=m,1,-1
       Seq(m-n+1)=TRIM(checks%children(n))
       line=TRIM(checks%children(n))
       CALL MY_Parse(line,strngs)       
       p=0
       DO WHILE(p < SIZE(strngs))
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
    ALLOCATE(Secondary_Seq (count_j) % line(c))
    c=0
    DO n=1,m
       line=Seq(n)
       CALL MY_Parse(line,strngs)
       p=0
       DO WHILE(p < SIZE(strngs))
          p=p+1
          IF(strngs(p) /= 'x') THEN
             c=c+1
             Secondary_Seq(count_j) % line (c)=strngs(p)
          ELSE
             IF(p == SIZE(strngs)) THEN
                errmsg_f='Sequence ends. A integer is expected after a ''x'' '
                CALL Add_Errors(-1,errmsg_f)
                RETURN
             END IF
             CALL SP_Getnum(strngs(p+1), rept,iflag)
             IF(iflag /= 0) THEN
                errmsg_f='Reading error. A integer is expected after a ''x'' '
                CALL Add_Errors(-1,errmsg_f)
                RETURN
             END IF
             Secondary_Seq(count_j) % line (c+1:c+rept)=strngs(p-1)
             c=c+rept
             p=p+1
          END IF
       END DO
    END DO
    DEALLOCATE(Seq)
  END SUBROUTINE Join
  SUBROUTINE Patch_it(name)
    CHARACTER(len=*) :: name
    TYPE(Branch), SAVE :: checks
    CHARACTER(len=max_pars) :: line,linea
    INTEGER :: n,nword,iflag

    CALL Check_Tree(name,checks)
    IF(.NOT. ASSOCIATED(checks%children)) RETURN
    ALLOCATE(patches(SIZE(checks%children)))
    DO n=1,SIZE(checks%children)
       line=TRIM(checks%children(n))
       CALL MY_Parse(line,strngs)
       linea=strngs(1)
       nword=SIZE(strngs)
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
  INCLUDE 'Store_TpgPar.f90'
END MODULE PARAMETERS_Mod
