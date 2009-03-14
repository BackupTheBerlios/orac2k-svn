MODULE POLAR_Mod

!!$***********************************************************************
!!$   Time-stamp: <2007-12-06 09:34:04 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Mon Dec  3 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program orac2k ----*


!!$======================== DECLARATIONS ================================*

  USE Node
  USE INPUT_Mod, ONLY: Read_String, Parser, err_open,err_end,err_unr&
       &,err_fnf,err_args
  IMPLICIT none
  PRIVATE
  PUBLIC :: Read_it,Polar_,Polar__, Init

  CHARACTER(8), SAVE :: Target_Res='tip3    '
  LOGICAL, SAVE :: Polar__, erfc_spline
  REAL(8), SAVE :: erfc_bin=0.1D0, alphal, rspoff=10.0D0, rneim,&
       & rntolm, rcutm, rneil, rntoll, rcutl, rneih, rntolh, rcuth,&
       & dx_Pr=0.1D0, max_Pr=20.0D0
  INTEGER, SAVE :: nccx=0,nccy=0,nccz=0, natom,n_write=1,nats=0&
       &,natoms_Tot=0,ntot=1
  LOGICAL, DIMENSION (:), ALLOCATABLE, SAVE :: mask
  INTEGER, DIMENSION (:), ALLOCATABLE, SAVE :: atoms
  REAL(8), ALLOCATABLE, SAVE :: Pr(:)
  INTEGER, ALLOCATABLE, SAVE :: nPr(:)
  INTEGER, ALLOCATABLE, SAVE :: mapnl(:),worka(:),tag_bndg(:)
  REAL(8), ALLOCATABLE, SAVE :: cpu_h(:)
  CHARACTER(80), SAVE :: filename_pol,filename_pdb,filename_pol2,filename_Pr
  INTEGER, SAVE :: kelectric=0,kelectric2=0,kpdb=0,kpr=0
  INTEGER, SAVE :: counter=0
  CHARACTER(8), SAVE :: file_format='xplor   '
  LOGICAL, SAVE :: do_Pr=.FALSE.
  INTEGER, SAVE :: nstart,nend,nstart_a,nend_a,nlocal_a,noeds=0,nprocs=1&
       &,ncube,nodex,nodey,nodez,nstart_h,nend_h,nlocal_h &
       &,nstart_ah,nend_ah,nlocal_ah,nstart_1,nend_1,nlocal_1&
       &,nstart_ex0,nend_ex0,nlocal_ex0,nstart_ex,nend_ex,nlocal_ex&
       &,nstart_2,nend_2,nlocal_2,fudgec,nbyte=4,ncpu_h
  TYPE :: Molec
     INTEGER :: m
     INTEGER, ALLOCATABLE :: p(:)
  END TYPE Molec
  TYPE Coord
     INTEGER :: m
     REAL(8) :: x,y,z
     REAL(8) :: charge
     REAL(8) :: dip(3)
  END TYPE Coord
  TYPE(Coord), ALLOCATABLE, SAVE :: Pt(:)
  TYPE(Molec), ALLOCATABLE, SAVE :: mol(:)
  REAL(8),SAVE :: a=0.0D0,b=0.0D0,c=0.0D0
  REAL(8),SAVE :: ax=0.0D0,bx=0.0D0,cx=0.0D0,alf=0.0D0,bet=0.0D0,gam=0.0D0
CONTAINS
  SUBROUTINE Init(natom0,protl,nprot,prsymb,res1,res2)
    
    IMPLICIT NONE
    CHARACTER(8) :: prsymb(*)
    INTEGER :: res1(*),res2(*),protl(*),nprot
    INTEGER :: natom0
!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: i,na,count,ii,m,n,nx,ny,nz
    LOGICAL :: ok
    REAL(8) :: kT

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    natom=natom0
    
    ALLOCATE(mask(natom))

    mask=.FALSE.
    ok=.FALSE.
    DO i=1,natom0
       IF(Target_Res == prsymb(res2(i))) THEN
          mask(i)=.TRUE.
          ok=.TRUE.
       END IF
    END DO
    na=0
    count=0
    DO ii=1,nprot
       m=protl(1+na)
       IF(m /=0) THEN
          i=protl(1+na+1)
          IF(Mask(i)) THEN
             count=count+1
          END IF
       END IF
       na=na+m+1
    END DO
    ALLOCATE(mol(count))
    na=0
    count=0
    DO ii=1,nprot
       m=protl(1+na)
       IF(m /=0) THEN
          i=protl(1+na+1)
          IF(Mask(i)) THEN
             count=count+1
             mol(count) % m = m
             ALLOCATE(mol(count) % p(m))
             mol(count) % p(:) = protl(1+na+1:1+na+m)
          END IF
       END IF
       na=na+m+1
    END DO

    ALLOCATE(Pt(0:nccx*nccy*nccz-1))
    Pt(:) % dip (1)  = 0.0D0
    Pt(:) % dip (2)  = 0.0D0
    Pt(:) % dip (3)  = 0.0D0
    Pt(:) % charge = 0.0D0
    Pt(:) % m = 0
    

    DO nz=0,nccz-1
       DO ny=0,nccy-1
          DO nx=0,nccx-1
             m=nccy*nccx*nz+nccx*ny+nx
             Pt(m) % x =DBLE(nx)/DBLE(nccx)
             Pt(m) % y =DBLE(ny)/DBLE(nccy)
             Pt(m) % z =DBLE(nz)/DBLE(nccz)
          END DO
       END DO
    END DO
    IF(natoms_Tot /=0) filename_pdb='POLAR_FILE_p.pdb'
    SELECT CASE(file_format)
    CASE DEFAULT
       filename_pol='POLAR_FILE_e.cube'
       filename_pol2='POLAR2_FILE.cube'
    CASE('xplor')
       filename_pol='POLAR_FILE_e.xplor'
       filename_pol2='POLAR2_FILE.xplor'
    END SELECT
  END SUBROUTINE Init
  SUBROUTINE Polar_(xpa,ypa,zpa,coa,aoc)
    IMPLICIT NONE
    REAL(8) :: coa(3,3),aoc(3,3),xpa(*),ypa(*),zpa(*)
    REAL(8) :: a0,b0,c0,gam0,bet0,alf0
    
    counter=counter+1
    CALL rotb(a0,b0,c0,alf0,bet0,gam0,coa)
    ax=ax+a0
    bx=bx+b0
    cx=cx+c0
    alf=alf+alf0
    bet=bet+bet0
    gam=gam+gam0
    CALL Get_Pol

    IF(MOD(counter,1) == 0) WRITE(77,'(i5,2x,3e14.6,2x,f12.5)') counter,Pt(5433) % dip(1)&
         &/DBLE(counter),Pt(5433) % dip(2)/DBLE(counter),Pt(5433) %&
         & dip(3)/DBLE(counter),DBLE(Pt(5433) % m)/DBLE(counter)
  CONTAINS
    SUBROUTINE Get_Pol
      IMPLICIT NONE 
      INCLUDE 'parst.h'
      INCLUDE 'cpropar.h'
      INTEGER :: n,ii,i,nx,ny,nz,m
      REAL(8) :: xcm,ycm,zcm,dip(3),charge,xc,yc,zc,x1,y1,z1,xcmm&
           &,ycmm,zcmm,xd,yd,zd

      
      xd=xpa(1)
      yd=ypa(1)
      zd=zpa(1)
      x1=xd*0.5D0
      y1=yd*0.5D0
      z1=zd*0.5D0
      nx=INT(DBLE(nccx)*(x1-ANINT(x1)+0.5D0))
      ny=INT(DBLE(nccy)*(y1-ANINT(y1)+0.5D0))
      nz=INT(DBLE(nccz)*(z1-ANINT(z1)+0.5D0))
      m=nccy*nccx*nz+nccx*ny+nx
      DO n=1,SIZE(mol)
         xcm=0.0D0
         ycm=0.0D0
         zcm=0.0D0
         dip(:)=0.0D0
         charge=0.0D0
         DO ii=1,mol(n) % m
            i=mol(n) % p (ii) 
            xcm=xcm+xpa(i)
            ycm=ycm+ypa(i)
            zcm=zcm+zpa(i)
            xc=coa(1,1)*xpa(i)+coa(1,2)*ypa(i)+coa(1,3)*zpa(i)
            yc=coa(2,1)*xpa(i)+coa(2,2)*ypa(i)+coa(2,3)*zpa(i)
            zc=coa(3,1)*xpa(i)+coa(3,2)*ypa(i)+coa(3,3)*zpa(i)
            dip(1)=dip(1)+xc*chrge(i)
            dip(2)=dip(2)+yc*chrge(i)
            dip(3)=dip(3)+zc*chrge(i)
            charge=charge+chrge(i)
         END DO
         xcm=xcm/DBLE(mol(n) % m)
         ycm=ycm/DBLE(mol(n) % m)
         zcm=zcm/DBLE(mol(n) % m)
         xcmm=coa(1,1)*xcm+coa(1,2)*ycm+coa(1,3)*zcm
         ycmm=coa(2,1)*xcm+coa(2,2)*ycm+coa(2,3)*zcm
         zcmm=coa(3,1)*xcm+coa(3,2)*ycm+coa(3,3)*zcm
         x1=xcm*0.5D0
         y1=ycm*0.5D0
         z1=zcm*0.5D0
         nx=INT(DBLE(nccx)*(x1-ANINT(x1)+0.5D0))
         ny=INT(DBLE(nccy)*(y1-ANINT(y1)+0.5D0))
         nz=INT(DBLE(nccz)*(z1-ANINT(z1)+0.5D0))
         
         m=nccy*nccx*nz+nccx*ny+nx
         Pt(m) % dip = Pt(m) % dip + dip
         Pt(m) % charge = Pt(m) % charge + charge
         Pt(m) % m = Pt(m) % m + 1
      END DO
    END SUBROUTINE Get_Pol
  END SUBROUTINE Polar_
  INCLUDE 'POLAR__Read.f90'
END MODULE POLAR_Mod
