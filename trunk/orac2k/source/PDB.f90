MODULE PDB

!!$***********************************************************************
!!$   Time-stamp: <2008-03-07 12:41:40 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Wed Oct 31 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*
  
  USE Constants  
  USE ReadStore
  USE INPUT_Mod, ONLY: Read_String, Parser, err_open,err_end,err_unr&
       &,err_fnf,err_args
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors,errmsg_f
  USE Strings, ONLY: My_Fam, MyRead, MyPutNum
  USE STRPAK
  USE Node
  IMPLICIT none
  PRIVATE
  PUBLIC :: Initialize,PDB_,Read_it,Write_it
  TYPE :: AtomPDB
     INTEGER :: Serial
     REAL(8) :: x,y,z
     CHARACTER(len=4) :: AtmName
  END TYPE AtomPDB
  TYPE ResiduePDB
     INTEGER :: No
     CHARACTER(len=3) :: ResName
     CHARACTER(len=1) :: Chain
     TYPE(AtomPdb), DIMENSION(:), ALLOCATABLE :: atm
  END TYPE ResiduePDB
  TYPE(ResiduePDB), ALLOCATABLE, SAVE :: ResPdb(:)

  LOGICAL, SAVE :: PDB_=.FALSE.
  INTEGER, SAVE :: kpdb
  CHARACTER(len=max_char), SAVE :: PDB_Filename
  CHARACTER(len=max_char), ALLOCATABLE, SAVE  :: PDB_template(:)
CONTAINS
  SUBROUTINE Initialize
    INTEGER :: n
    IF(.NOT. ReadStore_(PDB_Filename)) CALL Print_Errors()
    ALLOCATE(PDB_Template(SIZE(RS__string)))
    PDB_Template=RS__string
    CALL ReadStore__Delete
    IF(.NOT. PDB__Read('Template', PDB_Template)) CALL Print_Errors()
  END SUBROUTINE Initialize
  SUBROUTINE Compute
  END SUBROUTINE Compute
  SUBROUTINE write_it(fstep,kpdb,beta,xp0,yp0,zp0,ntap, res, co)
    INTEGER :: ntap,kpdb
    INTEGER :: res(*)
    REAL(8) ::  xp0(*),yp0(*),zp0(*),fstep
    CHARACTER(7) :: beta(*)
    REAL(8), OPTIONAL :: co(3,3)

    INTEGER :: i,j,k,l,m,ndx
    REAL(8) :: xb,yb,zb,charge,dx,dy,dz,dr,xd,yd,zd,sum
    INTEGER :: ResNumb
    CHARACTER(4) :: bet,bet2
    CHARACTER(3) :: ResName
    CHARACTER(1) :: Chain
    REAL(8) ::  a,b,c,alf,bett,gamm
    
    WRITE(kpdb,2) fstep
    IF(PRESENT(co)) THEN
       CALL rotb(a,b,c,alf,bett,gamm,co)
       WRITE(kpdb,5) a,b,c,alf,bett,gamm
    END IF
    charge=0.0D0
    k=0
    dr=0.0D0

    DO i=1,ntap
       xb=xp0(i)
       yb=yp0(i)
       zb=zp0(i)
       bet=ADJUSTL(beta(i))
       CALL low_up(bet,4)
       bet2=ADJUSTR(bet)
       Chain=ResPdb(res(i)) % Chain
       ResName=ResPdb(res(i)) % ResName
       ResNumb=ResPdb(res(i)) % No
       CALL low_up(ResName,3)
       WRITE(kpdb,1) 'ATOM ',i,bet2,ResName,Chain,ResNumb,xb,yb,zb,dr,DFLOAT(k)
    END DO
    
1   FORMAT(a5,i6,1x,a4,x,a3,1x,a1,i4,4x,3f8.3,f8.4,f4.1)
2   FORMAT('REMARK   1 Configuration at time step ',f11.2,'      ',&
         &'                 ')
5     FORMAT('CRYST1',3f9.3,3f7.2,' P  1          1')
  END SUBROUTINE write_it
  
  INCLUDE 'PDB__Read.f90'
END MODULE PDB
