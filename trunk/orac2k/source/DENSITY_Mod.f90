MODULE DENSITY_Mod

!!$***********************************************************************
!!$   Time-stamp: <2009-07-21 10:01:48 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Oct 14 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f
  USE PBC_Mod
  USE INPUT_Mod, ONLY: Read_String, Parser, err_open,err_end,err_unr&
       &,err_fnf,err_args
  USE PDB, ONLY: PDB_out=>Write_it, PDB_
  Use Tetrahedra
  Use Neighbors
  IMPLICIT None

  PRIVATE
  PUBLIC Initialize, Initialize_Ar, Initialize_Par, Compute, Write_it&
       &, Read_It, Density_Calc, n_write
  REAL(8), DIMENSION (:,:,:), ALLOCATABLE, SAVE :: Density,Dens2,Loc_Dens
  REAL(8), DIMENSION (:,:,:), ALLOCATABLE, SAVE :: Ndens
  REAL(8), DIMENSION (:), ALLOCATABLE, SAVE :: atmass
  REAL(8), DIMENSION (:), ALLOCATABLE, SAVE :: xc,yc,zc
  REAL(8), SAVE :: a=0.0D0,b=0.0D0,c=0.0D0,alph,bet,gamm,coa(3,3),oca(3&
       &,3),volume,xcm,ycm,zcm,Dvolume,rho=0.0D0,unitfluc,UnitDens,Total_Mass
  INTEGER, SAVE :: natom,n_write=1,ncx=0,ncy=0,ncz=0,ntot=1,nats=0&
       &,natoms_Tot=0,natom_Slt
  INTEGER, SAVE :: kdensity=0,kdens2=0,kpdb=0
  INTEGER, SAVE :: counter=0,node=0
  LOGICAL, SAVE :: Density_Calc=.FALSE.,Density_Avg=.FALSE.
  REAL(8), SAVE :: x_ad=0.0D0,y_ad=0.0D0,z_ad=0.0D0
  LOGICAL, DIMENSION (:), ALLOCATABLE, SAVE :: mask,CellMask
  INTEGER, DIMENSION (:), ALLOCATABLE, SAVE :: atoms
  CHARACTER(80), SAVE :: filename_den,filename_den2,filename_pdb
  CHARACTER(8), SAVE :: Target_Res,file_format='cube'
  CHARACTER(7), ALLOCATABLE, SAVE :: beta(:)
  TYPE PDB_t
     CHARACTER(5) :: beta
     CHARACTER(5) :: symb
     INTEGER :: res
     REAL(8) :: x,y,z
  END TYPE PDB_T
  TYPE(PDB_t), DIMENSION(:), ALLOCATABLE, SAVE :: coords
  LOGICAL, SAVE :: Other_Density=.FALSE.,Dens_Histo=.False.
  Real(8), Save :: ddx,ddy,ddz,dx,dy,dz
  Real(8), Allocatable :: xpa(:),ypa(:),zpa(:),Dens_Histo_svol(:,:)&
       &,Dens_Histo_gofr(:)
  Integer, Allocatable :: Dens_Histo_numb(:),idx(:)
  Character(120), Save :: Filename_histo=' ',Filename_histo_cos=' '
  Integer, Save :: kdenhisto=0,kdenhisto_cos=0,dens_histo_size=400,mx,my,mz&
       &,Dens_histo_Size_cos=50
  Real(8), Save :: dens_histo_rmax=10.0_8,dens_histo_bin,Dens_histo_Bin_cos
  Real(8), Allocatable :: Dens_Histo_gofr_ac(:,:)&
       &,Dens_Histo_gofr_cos(:,:,:)
  Integer, Allocatable :: iac(:)
  Type :: Amino
     Character(len=3), Allocatable :: res(:)
     Character(len=3), Allocatable :: atom(:)
     Integer :: Type
  End type Amino
  Type(Amino), Save :: polar,nonpolar
  Logical, Save :: PerResidue=.True.,NoHydrogen=.FALSE.,fixmolecule=.False.
CONTAINS
!!$---- This subroutine is part of the program ORAC ----*

  SUBROUTINE Initialize(nodea,prsymb,betab,betac,res1,res2,mass,ntap&
       &,nato_slt,co,oc,volumea)

!!$======================== DECLARATIONS ================================*

    IMPLICIT none
    INCLUDE 'unit.h'

!!$----------------------------- ARGUMENTS ------------------------------*

    CHARACTER(8) :: prsymb(*)
    CHARACTER(7) :: betab(*),betac(*)
    INTEGER :: res1(*),res2(*),ntap,nodea,nato_slt
    REAL(8) :: mass(*),co(3,3),oc(3,3),volumea

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: i,na,count0,ii,m,n
    LOGICAL :: ok
    REAL(8) :: kT
    Character(len=3) :: char1,char2,reschar1,char0
             
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    Allocate(nonpolar%res(10))
    nonpolar%res=(/'ala','cys','gly','ile','leu','met','phe','pro','trp','val' /)
    nonpolar%type=1
    Allocate(polar%res(10))
    Allocate(polar%atom(28))
    polar%res=(/'arg','asn','asp','glu','gln','lys','ser','thr','tyr','hid' /)
    polar%atom=(/'h  ','hc ','ht ','hs ','c  ','cc ','cd '&
         &,'cm ','n  ','nr1','nr2','nr3','nh1','nh2','nh3','nc2'&
         &,'nph','ny ','np ','o  ','ob ','oc& 
         & ','os ','ot ','om ','oh1','s  ','cal'/)
    polar%type=2

    node=nodea
    kT=300.0D0*gascon/1000.0D0
    unitfluc=efact/unitp/kT
    unitdens=unitm/unitl**3/1000.0D0
    coa=0.0D0
    oca=0.0D0

    If(a == 0.0D0 .Or. b == 0.0D0 .Or. c == 0.0D0) Then
       coa=co; oca=oc
       volume=volumea
    Else
       coa(1,1)=a/2.0D0
       coa(2,2)=b/2.0D0
       coa(3,3)=c/2.0D0
       
       oca(1,1)=1.0D0/coa(1,1)
       oca(2,2)=1.0D0/coa(2,2)
       oca(3,3)=1.0D0/coa(3,3)
       volume=a*b*c
    End If
    natom=ntap
    natom_slt=nato_slt

    IF(Dens_Histo) Then
       ddx=boxl/Dble(ncx)
       ddy=boxl/Dble(ncy)
       ddz=boxl/Dble(ncz)
       Dens_histo_bin=Dens_histo_rmax/Dble(Dens_histo_size)
       Dens_histo_Bin_cos=2.0D0/Dble(Dens_histo_size_cos)
       
       Allocate(xpa(natom),ypa(natom),zpa(natom),idx(natom))
       Allocate(Dens_Histo_gofr(Dens_histo_size+1))
       Allocate(Dens_Histo_svol(Dens_histo_size+1,0:2))
       Allocate(Dens_Histo_numb(Dens_histo_size+1))
       allocate(iac(natom),Dens_Histo_gofr_ac(Dens_histo_size+1,0:2))

       allocate(Dens_Histo_gofr_cos(Dens_histo_size_cos+1,Dens_histo_size+1,0:2))
       
       If(PerResidue) Then
          iac=polar%type
          Do n=1,natom_slt
             char1=Trim(prsymb(res2(n)))
             char0=Trim(betac(n))
             If(char0(1:1) =='h' .And. NoHydrogen) Then
                iac(n)=0
                Cycle
             End If
             Do m=1,Size(NonPolar%res)
                char2=Trim(NonPolar%res(m))
                If(char1 == char2) Then
                   iac(n)=NonPolar%type
                   Exit
                End If
             End Do
          End Do
       Else
          iac=nonpolar%type
          Do n=1,natom_slt
             reschar1=Trim(prsymb(res2(n)))
             char1=Trim(betac(n))
             If(char1(1:1) =='h' .And. NoHydrogen) Then
                iac(n)=0
                Cycle
             End If
             Do m=1,Size(Polar%atom)
                char2=Trim(Polar%atom(m))
                If(char1 == char2) Then
                   iac(n)=Polar%type
                   Exit
                End If
             End Do
             If(reschar1 == 'pro' .And. Trim(char1) == 'n') iac(n)=NonPolar%type
          End Do          
       End If
       Do n=1,natom_slt
          Write(*,'(i6,2x,a6,4x,i6)') n,betac(n),iac(n)
       End Do

       Write(*,*) 'iacs =',Count(iac(1:natom_slt) == polar%type),Count(iac(1:natom_slt) ==&
            & nonpolar%type),Count(iac(1:natom_slt) == 0)
       idx(natom_slt+1:)=2
       idx(1:natom_slt)=1
       Dens_Histo_gofr_ac=0.0_8
       Dens_Histo_gofr=0.0_8
       Dens_Histo_svol=0.0_8
       Dens_Histo_numb=0

       mx=10;my=10;mz=10;
       dx=boxl/Dble(mx); dy=boxl/Dble(my); dz=boxl/Dble(mz);
       Allocate(beta(natom))
       beta=betac(1:natom)
    End IF
    ALLOCATE(atmass(ntap),mask(ntap))

    mask=.FALSE.

    ok=.FALSE.
    Total_Mass=0.0D0
    DO i=1,ntap
       IF(Target_Res == prsymb(res2(i))) THEN
          mask(i)=.TRUE.
          ok=.TRUE.
       END IF
       atmass(i)=mass(i)
       Total_Mass=Total_Mass+mass(i)
    END DO
    IF(.NOT. ok) THEN
       mask=.TRUE.
       WRITE(*,*) '=============>  Hope this is ok, you are computing bare atomic density '
    END IF

    IF(natoms_Tot /= 0) THEN
       na=0
       count0=0
       DO n=1,natoms_Tot
          m=atoms(na+1)
          count0=count0+m
          na=na+1+m
       END DO
       ALLOCATE(coords(count0),xc(ntap),yc(ntap),zc(ntap))
       xc=0.0D0
       yc=0.0D0
       zc=0.0D0
       na=0
       count0=0
       DO n=1,natoms_Tot
          m=atoms(na+1)
          DO ii=1,m
             i=atoms(na+1+ii)
             count0=count0+1
             coords(count0) % beta = TRIM(betab(i))
             coords(count0) % symb = TRIM(prsymb(res2(i)))
             coords(count0) % res = res1(i)
          END DO
          na=na+1+m
       END DO
       ALLOCATE(beta(SIZE(Coords)))
       beta=ADJUSTL(Coords(:) % beta)
    END IF

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*
  END SUBROUTINE Initialize

  SUBROUTINE Initialize_Ar

!!$======================== DECLARATIONS ================================*

    IMPLICIT none

!!$------------------------- LOCAL VARIABLES ----------------------------*


!!$----------------------- EXECUTABLE STATEMENTS ------------------------*
    
    ALLOCATE(Density(ncx,ncy,ncz))
    ALLOCATE(Dens2(ncx,ncy,ncz),ndens(ncx,ncy,ncz),loc_dens(ncx,ncy,ncz))
    Density=0.0D0; Dens2=0.0D0
    Dvolume=Volume/DBLE(ncx*ncy*ncz)
    Allocate(CellMask(ncx*ncy*ncz))

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*
  END SUBROUTINE Initialize_Ar
  SUBROUTINE Initialize_Par(a0,b0,c0,ncx0,ncy0,ncz0,n_write0&
       &,file_format0,atoms0,nats0,natoms_Tot0,filename_den0&
       &,filename_den20,filename_pdb0,Target_Res0)
    REAL(8) :: a0,b0,c0
    INTEGER :: ncx0,ncy0,ncz0,n_write0,atoms0(:),nats0,natoms_Tot0
    CHARACTER(*) :: filename_den0,filename_den20&
         &,filename_pdb0
    CHARACTER(8) :: file_format0,target_res0
    LOGICAL :: exist


    ALLOCATE(atoms(SIZE(atoms0)))
    atoms=atoms0
    a=a0
    b=b0
    c=c0
    ncx=ncx0
    ncy=ncy0
    ncz=ncz0
    n_write=n_write0
    file_format=file_format0
    nats=nats0
    natoms_tot=natoms_tot0
    filename_den=filename_den0
    filename_den2=filename_den20
    filename_pdb=filename_pdb0
    target_res=target_res0
    SELECT CASE(file_format)
    CASE DEFAULT
       INQUIRE(FILE=filename_den,EXIST=exist)
       IF(exist) THEN
          CALL openf(kdensity,filename_den,'FORMATTED','OLD',0)
       ELSE
          CALL openf(kdensity,filename_den,'FORMATTED','NEW',0)
       END IF
       INQUIRE(FILE=filename_den2,EXIST=exist)
       IF(exist) THEN
          CALL openf(kdens2,filename_den2,'FORMATTED','OLD',0)
       ELSE
          CALL openf(kdens2,filename_den2,'FORMATTED','NEW',0)
       END IF
    CASE('xplor')
       INQUIRE(FILE=filename_den,EXIST=exist)
       IF(exist) THEN
          CALL openf(kdensity,filename_den,'FORMATTED','OLD',0)
       ELSE
          CALL openf(kdensity,filename_den,'FORMATTED','NEW',0)
       END IF
       INQUIRE(FILE=filename_den2,EXIST=exist)
       IF(exist) THEN
          CALL openf(kdens2,filename_den2,'FORMATTED','OLD',0)
       ELSE
          CALL openf(kdens2,filename_den2,'FORMATTED','NEW',0)
       END IF
    END SELECT
    IF(natoms_Tot /=0) THEN
       INQUIRE(FILE=filename_pdb,EXIST=exist)
       IF(exist) THEN
          CALL openf(kpdb,filename_pdb,'FORMATTED','OLD',0)
       ELSE
          CALL openf(kpdb,filename_pdb,'FORMATTED','NEW',0)
       END IF
    END IF
  END SUBROUTINE Initialize_Par
  SUBROUTINE Compute(xp0,yp0,zp0,co,oc,Real_Volume,My_Density)

!!$======================== DECLARATIONS ================================*

    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

    REAL(8) :: xp0(*),yp0(*),zp0(*),Real_Volume,co(3,3),oc(3,3)
    REAL(8), OPTIONAL :: My_Density(*)

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: nx,ny,nz,na,m,ii,i,qq,j,k
    REAL(8) :: x1,y1,z1
    INTEGER :: n,nc2x,nc2y,nc2z
    REAL(8) :: xpb,ypb,zpb,aux1,aux2,fract
    Real(8), Save :: dvol
    Integer :: o
    Integer, Allocatable :: loc_nhist(:)
    Real(8), Allocatable :: loc_hist(:)

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

!!$*=======================================================================
!!$*     Compute chain list for system
!!$*=======================================================================

    counter=counter+1
    ndens=0.0D0; loc_dens=0.0D0
    IF(counter == 1) THEN
       IF(PRESENT(My_Density)) THEN
          Other_Density=.TRUE.
       END IF
    END IF
    xcm=0.0D0 ; ycm=0.0D0 ; zcm=0.0D0 

    DO n=1,natom
       xcm=xcm+xp0(n)
       ycm=ycm+yp0(n)
       zcm=zcm+zp0(n)
    END DO
    xcm=xcm/DBLE(natom)
    ycm=ycm/DBLE(natom)
    zcm=zcm/DBLE(natom)
    xp0(1:natom)=xp0(1:natom)-xcm
    yp0(1:natom)=yp0(1:natom)-ycm
    zp0(1:natom)=zp0(1:natom)-zcm

    DO n=1,natom
       IF(Mask(n)) THEN
          IF(.NOT. Other_Density) THEN
             aux1=atmass(n)
          ELSE
             aux1=My_Density(n)
          END IF
          xpb=oca(1,1)*xp0(n)+oca(1,2)*yp0(n)+oca(1,3)*zp0(n)
          ypb=oca(2,1)*xp0(n)+oca(2,2)*yp0(n)+oca(2,3)*zp0(n)
          zpb=oca(3,1)*xp0(n)+oca(3,2)*yp0(n)+oca(3,3)*zp0(n)
          x1=xpb*0.5D0
          y1=ypb*0.5D0
          z1=zpb*0.5D0
          nx=INT(DBLE(ncx-1)*(x1-ANINT(x1)+0.5D0)+0.5D0)
          ny=INT(DBLE(ncy-1)*(y1-ANINT(y1)+0.5D0)+0.5D0)
          nz=INT(DBLE(ncz-1)*(z1-ANINT(z1)+0.5D0)+0.5D0)
          ndens(nx+1,ny+1,nz+1)=ndens(nx+1,ny+1,nz+1)+1.0D0
          loc_dens(nx+1,ny+1,nz+1)=loc_dens(nx+1,ny+1,nz+1)+aux1
       END IF
    END DO

    IF(.NOT. Other_Density) THEN       
       Density=Density+loc_Dens
       Dens2=Dens2+loc_Dens**2
    ELSE
       DO k=1,nz
          DO j=1,ncy
             DO i=1,ncx
                IF(ndens(i,j,k) /= 0.0D0) THEN
                   aux1=loc_dens(i,j,k)/ndens(i,j,k)
                   Density(i,j,k)=Density(i,j,k)+aux1
                   Dens2(i,j,k)=Dens2(i,j,k)+aux1**2
                END IF
             END DO
          END DO
       END DO
    END IF

    IF(Dens_Histo) Then
       xpa=oc(1,1)*xp0(1:natom)+oc(1,2)*yp0(1:natom)+oc(1,3)*zp0(1:natom)
       ypa=oc(2,1)*xp0(1:natom)+oc(2,2)*yp0(1:natom)+oc(2,3)*zp0(1:natom)
       zpa=oc(3,1)*xp0(1:natom)+oc(3,2)*yp0(1:natom)+oc(3,3)*zp0(1:natom)

       If(.Not. Neighbors_(Dens_histo_rmax,mx,my,mz,co)) Call Print_Errors()
       If(.Not. Neighbors__Particles(xpa,ypa,zpa)) Call Print_Errors()
       Allocate(Loc_Hist(Dens_histo_size),loc_nhist(Dens_histo_size))

!!$
!!$-- Calls in this and only order
!!$

       Call ShellDensity
       Call ShellVolume
    End IF
    
    rho=rho+Total_Mass/Real_Volume
    IF(natoms_Tot /= 0) THEN
       na=0
       DO n=1,natoms_Tot
          m=atoms(na+1)
          DO ii=1,m
             i=atoms(1+na+ii)
             xc(i)=xc(i)+xp0(i)
             yc(i)=yc(i)+yp0(i)
             zc(i)=zc(i)+zp0(i)
          END DO
          na=na+m+1          
       END DO
    END IF

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

  Contains
    Subroutine ShellDensity
      Integer :: p,q,r,iv,jv,kv,l,numcell,numcell1
      Integer :: nhisto,o,nnx,nny,nnz,n,mycount,u
      Real(8) :: xa,ya,za,xc0,yc0,zc0,rsq,rsq_Min,rsp,rsq_max,xpi,ypi,zpi
      Real(8) :: vhb(3,4),xf(3),yf(3),zf(3),MyC(3)
      Integer :: nang,ll
      Real(8) :: Norm,cos1,cos2,cos3,cos4


      
      mycount=0
      loc_hist=0.0D0
      rsq_max=Dens_Histo_rmax**2
      CellMask=.False.
      Do n=natom_slt+1,natom,3
         xpi=xpa(n)
         ypi=ypa(n)
         zpi=zpa(n)
         Do m=n,n+2
            xa=xpa(m)
            ya=ypa(m)
            za=zpa(m)
            xf(m-n+1)=co(1,1)*xa+co(1,2)*ya+co(1,3)*za
            yf(m-n+1)=co(2,1)*xa+co(2,2)*ya+co(2,3)*za
            zf(m-n+1)=co(3,1)*xa+co(3,2)*ya+co(3,3)*za
         End Do
         Call Tetrahedra_(xf(1),yf(1),zf(1),xf(2),yf(2),zf(2)&
              &,xf(3),yf(3),zf(3),vhb)
         p=Chain_xyz(n) % i
         q=Chain_xyz(n) % j
         r=Chain_xyz(n) % k
         numcell1=r+mz*(q+my*p)+1
         CellMask(numcell1)=.True.
         rsq_Min=1.0D10
         Do o=1,SIZE(Ind_xyz)
            iv=Ind_xyz(o) % i
            jv=Ind_xyz(o) % j
            kv=Ind_xyz(o) % k
            nnx=mod(mod(p+iv,mx)+mx,mx)
            nny=mod(mod(q+jv,my)+my,my)
            nnz=mod(mod(r+kv,mz)+mz,mz)
            numcell=nnz+mz*(nny+my*nnx)+1
            l=Head_xyz(numcell)
            Do While(l > 0)
               If(idx(l) == 2) Then
                  l=Chain_xyz(l) % p
                  Cycle
               End If
               xa=xpa(l)-xpi
               xa=xa-2.0_8*PBC(xa)
               ya=ypa(l)-ypi
               ya=ya-2.0_8*PBC(ya)
               za=zpa(l)-zpi
               za=za-2.0_8*PBC(za)
               xc0=co(1,1)*xa+co(1,2)*ya+co(1,3)*za
               yc0=co(2,1)*xa+co(2,2)*ya+co(2,3)*za
               zc0=co(3,1)*xa+co(3,2)*ya+co(3,3)*za
               rsq=xc0*xc0+yc0*yc0+zc0*zc0
               If(rsq > rsq_Max) Then
                  l=Chain_xyz(l) % p
                  Cycle
               End If
               If(rsq < rsq_Min) Then
                  Rsq_Min=rsq
                  u=iac(l)
                  MyC=-(/xc0,yc0,zc0/)
                  ll=l
               End If
               l=Chain_xyz(l) % p
            End Do
         End Do
         Rsp=sqrt(Rsq_Min)
         nhisto=Int(Rsp/Dens_Histo_bin)+1
         If(nhisto <= Dens_histo_Size) Then
            Norm=Sqrt(Sum(MyC**2))
            MyC=MyC/Norm
            cos1=Sum(MyC*vhb(:,1))+1.0D0
            nang=Int(cos1/Dens_histo_Bin_cos)+1
            Dens_Histo_gofr_cos(nang,nhisto,u)=Dens_Histo_gofr_cos(nang,nhisto,u)+1.0D0
            cos2=Sum(MyC*vhb(:,2))+1.0D0
            nang=Int(cos2/Dens_histo_Bin_cos)+1
            Dens_Histo_gofr_cos(nang,nhisto,u)=Dens_Histo_gofr_cos(nang,nhisto,u)+1.0D0
            cos3=Sum(MyC*vhb(:,3))+1.0D0
            nang=Int(cos3/Dens_histo_Bin_cos)+1
            Dens_Histo_gofr_cos(nang,nhisto,u)=Dens_Histo_gofr_cos(nang,nhisto,u)+1.0D0
            cos4=Sum(MyC*vhb(:,4))+1.0D0
            nang=Int(cos4/Dens_histo_Bin_cos)+1
            Dens_Histo_gofr_cos(nang,nhisto,u)=Dens_Histo_gofr_cos(nang,nhisto,u)+1.0D0

            Dens_histo_gofr(nhisto)=Dens_histo_gofr(nhisto)+1.0_8
            Dens_Histo_gofr_ac(nhisto,u)=Dens_Histo_gofr_ac(nhisto,u)+1
         End If
      End Do
    End Subroutine ShellDensity
    Subroutine ShellVolume
      Integer :: p,q,r,iv,jv,kv,l,numcell,numcell1
      Integer :: nhisto,o,nnx,nny,nnz,u
      Real(8) :: xa,ya,za,xc0,yc0,zc0,rsq,rsq_Min,rsp,rsq_max,x2,y2&
           &,z2,dvol
      Real(8) :: ranf,xcc,ycc,zcc


      dvol=Real_Volume/Dble(ncx*ncy*ncz)
      
      
 
      loc_nhist=0
      loc_hist=0.0D0
      rsq_max=Dens_Histo_rmax**2
      xcc=0.0D0; ycc=0.0D0; zcc=0.0D0; 
      If(fixmolecule) Then
         xcc=ddx*(2.0D0*ranf()-1.0D0)
         ycc=ddy*(2.0D0*ranf()-1.0D0)
         zcc=ddz*(2.0D0*ranf()-1.0D0)
      End If
      Do k=1,ncz
         z1=(Dble(k)-0.5_8)*ddz+xcc
         Do j=1,ncy
            y1=(Dble(j)-0.5_8)*ddy+ycc
            Do i=1,ncx
               x1=(Dble(i)-0.5_8)*ddx+zcc
               
               x2=x1/dx
               y2=y1/dy
               z2=z1/dz
               p=INT(x2)+(SIGN(1.D0,x2-INT(x2))-1.)/2
               q=INT(y2)+(SIGN(1.D0,y2-INT(y2))-1.)/2
               r=INT(z2)+(sign(1.d0,z2-int(z2))-1.)/2
               p=MOD(MOD(p,mx)+mx,mx)
               q=MOD(MOD(q,my)+my,my)
               r=MOD(MOD(r,mz)+mz,mz)
               numcell1=r+mz*(q+my*p)+1
               If(.Not. CellMask(numcell1)) Cycle
               rsq_Min=1.0D10
               DO o=1,SIZE(Ind_xyz)
                  iv=Ind_xyz(o) % i
                  jv=Ind_xyz(o) % j
                  kv=Ind_xyz(o) % k
                  nnx=mod(mod(p+iv,mx)+mx,mx)
                  nny=mod(mod(q+jv,my)+my,my)
                  nnz=mod(mod(r+kv,mz)+mz,mz)
                  numcell=nnz+mz*(nny+my*nnx)+1
                  l=Head_xyz(numcell)
                  DO WHILE(l > 0)
                     If(idx(l) == 2) Then
                        l=Chain_xyz(l) % p
                        Cycle
                     End If
                     xa=xpa(l)-x1
                     xa=xa-2.0_8*PBC(xa)
                     ya=ypa(l)-y1
                     ya=ya-2.0_8*PBC(ya)
                     za=zpa(l)-z1
                     za=za-2.0_8*PBC(za)
                     xc0=co(1,1)*xa+co(1,2)*ya+co(1,3)*za
                     yc0=co(2,1)*xa+co(2,2)*ya+co(2,3)*za
                     zc0=co(3,1)*xa+co(3,2)*ya+co(3,3)*za
                     rsq=xc0*xc0+yc0*yc0+zc0*zc0
                     If(rsq > rsq_Max) Then
                        l=Chain_xyz(l) % p
                        Cycle
                     End If
                     If(rsq < rsq_Min) Then
                        Rsq_Min=rsq
                        u=iac(l)
                     End If
                     l=Chain_xyz(l) % p
                  End Do
               End DO
               Rsp=sqrt(Rsq_Min)
               nhisto=Int(Rsp/Dens_Histo_bin)+1
               If(nhisto <= Dens_histo_Size) Then
                  Dens_Histo_svol(nhisto,u)=Dens_Histo_svol(nhisto,u)+dvol 
               End If
            End Do
         End Do
      End Do
    End Subroutine ShellVolume
  END SUBROUTINE Compute
  SUBROUTINE Write_it(fstep)

!!$======================== DECLARATIONS ================================*

    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

    REAL(8) :: fstep

!!$------------------------- LOCAL VARIABLES ----------------------------*

    REAL(8) :: Dvolume,avogad=6.0225d23,unitcm=1.0D-8,fact,cob(3,3)&
         &,dummy=0.0D0,a0=0.529177249D0,a,b,c,alf,bet,gamm,avg,sigma&
         &,fact2,Mass_Avg,r,y,z,z0

    INTEGER :: one=1,i,j,k,count,ncx2,ncx_s,ncx_e,ncy2,ncy_s,ncy_e&
         &,ncz2,ncz_s,ncz_e,na,n,m,ii,ntap
    REAL(8), DIMENSION (:,:,:), ALLOCATABLE, SAVE :: Fluct
    Integer :: u
    Real(8) :: sum,rc

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    IF(.NOT. ALLOCATED(Fluct)) ALLOCATE(fluct(ncx,ncy,ncz))
    cob=coa
    cob(1,1)=(2.0D0/a0)*cob(1,1)/DBLE(ncx)
    cob(1,2)=(2.0D0/a0)*cob(1,2)/DBLE(ncy)
    cob(2,2)=(2.0D0/a0)*cob(2,2)/DBLE(ncy)
    cob(1,3)=(2.0D0/a0)*cob(1,3)/DBLE(ncz)
    cob(2,3)=(2.0D0/a0)*cob(2,3)/DBLE(ncz)
    cob(3,3)=(2.0D0/a0)*cob(3,3)/DBLE(ncz)
    xcm=-(coa(1,1)+coa(1,2)+coa(1,3))/a0
    ycm=-(coa(2,1)+coa(2,2)+coa(2,3))/a0
    zcm=-(coa(3,1)+coa(3,2)+coa(3,3))/a0

    Dvolume=Volume/DBLE(ncx*ncy*ncz)

    fact=1.0D0/avogad/unitcm**3/DBLE(counter)
    IF(Other_Density) fact=1.0D0/DBLE(counter)

    CLOSE(kdensity)
    OPEN(unit=kdensity,file=filename_den,form='FORMATTED',status='OLD')
    CLOSE(kdens2)
    OPEN(unit=kdens2,file=filename_den2,form='FORMATTED',status='OLD')

    fluct=0.0D0
    DO i=1,ncx
       DO j=1,ncy
          DO k=1,ncz
             IF(Density(i,j,k) /= 0.0D0) THEN
                IF(Other_Density) THEN
                   Mass_avg=Density(i,j,k)/DBLE(counter)
                   fluct(i,j,k)=(Dens2(i,j,k)/DBLE(counter)-Mass_Avg**2)
                ELSE
                   Mass_avg=Density(i,j,k)/DBLE(counter)
                   fluct(i,j,k)=(Dens2(i,j,k)/DBLE(counter)-Mass_Avg**2)/Mass_avg
                END IF
             END IF
          END DO
       END DO
    END DO
    fact2=unitfluc*DBLE(counter)/rho
    IF(Other_Density) fact2=1.0D0

    WRITE(*,*) '--------------- Averaged Density of system ='&
         &,unitdens*Rho/DBLE(Counter)

    
    SELECT CASE(file_format)
    CASE DEFAULT
!!$========================================================================
!!$--- Write volume Cube file
!!$========================================================================

       WRITE(kdensity,'(''ORAC CUBE FILE Volumes ='',2f15.5)') Volume&
            &,dvolume
       WRITE(kdensity,'(''ORAC CUBE FILE'')')
       WRITE(kdensity,'(i5,3f12.6)') one,xcm,ycm,zcm
       
       WRITE(kdensity,'(i5,3f12.6)') ncx,cob(1,1),cob(2,1),cob(3,1)
       WRITE(kdensity,'(i5,3f12.6)') ncy,cob(1,2),cob(2,2),cob(3,2)
       WRITE(kdensity,'(i5,3f12.6)') ncz,cob(1,3),cob(2,3),cob(3,3)
       
       WRITE(kdensity,'(i5,4f12.6)') one,dummy,dummy,dummy,dummy
       
       DO i=1,ncx
          DO j=1,ncy
             WRITE(kdensity,'(6e13.5)') fact*Density(i,j,1:ncz)
          END DO
       END DO

!!$========================================================================
!!$--- Write volume Cube file
!!$========================================================================

       WRITE(kdens2,'(''ORAC CUBE FILE Volumes ='',2f15.5)') Volume&
            &,dvolume
       WRITE(kdens2,'(''ORAC CUBE FILE'')')
       WRITE(kdens2,'(i5,3f12.6)') one,xcm,ycm,zcm
       
       WRITE(kdens2,'(i5,3f12.6)') ncx,cob(1,1),cob(2,1),cob(3,1)
       WRITE(kdens2,'(i5,3f12.6)') ncy,cob(1,2),cob(2,2),cob(3,2)
       WRITE(kdens2,'(i5,3f12.6)') ncz,cob(1,3),cob(2,3),cob(3,3)
       
       WRITE(kdens2,'(i5,4f12.6)') one,dummy,dummy,dummy,dummy
              
       DO i=1,ncx
          DO j=1,ncy
             WRITE(kdens2,'(6e13.5)') fact2*fluct(i,j,1:ncz)
          END DO
       END DO
    CASE('xplor')
       ncx2=ncx/2
       IF(ncx2*2 /= ncx) ncx2=ncx2+1
       ncx_s=-ncx2
       ncx_e=ncx-1-ncx2

       ncy2=ncy/2
       IF(ncy2*2 /= ncy) ncy2=ncy2+1
       ncy_s=-ncy2
       ncy_e=ncy-1-ncy2

       ncz2=ncz/2
       IF(ncz2*2 /= ncz) ncz2=ncz2+1
       ncz_s=-ncz2
       ncz_e=ncz-1-ncz2

       CALL rotb(a,b,c,alf,bet,gamm,coa)
       CALL Statistics(fact,Density,avg,sigma)
       WRITE(kdensity,*) 
       WRITE(kdensity,'(i8,1x,''!NTITLE'')') 2
       WRITE(kdensity,'('' REMARKS ORAC XPLOR FILE Volumes =''&
            &,2f15.5)') Volume,dvolume
       WRITE(kdensity,'('' REMARKS ORAC XPLOR FILE'')')


       WRITE(kdensity,'(9i8)') ncx,ncx_s,ncx_e,ncy,ncy_s,ncy_e,ncz,ncz_s,ncz_e
       WRITE(kdensity,'(6e12.5)') a,b,c,alf,bet,gamm
       WRITE(kdensity,'(''ZYX'')')
       DO k=1,ncz
          WRITE(kdensity,'(i8)') k-1
          WRITE(kdensity,'(6e12.5)') ((fact*Density(i,j,k),i=1,ncx),j=1,ncy)
       END DO
       WRITE(kdensity,'(i8)') -9999
       WRITE(kdensity,'(2(e12.4,1x))') avg,sigma


       CALL Statistics(fact2,Fluct,avg,sigma)
       WRITE(kdens2,*) 
       WRITE(kdens2,'(i8,1x,''!NTITLE'')') 2
       WRITE(kdens2,'('' REMARKS ORAC XPLOR FILE Volumes =''&
            &,2f15.5)') Volume,dvolume
       WRITE(kdens2,'('' REMARKS ORAC XPLOR FILE'')')


       WRITE(kdens2,'(9i8)') ncx,ncx_s,ncx_e,ncy,ncy_s,ncy_e,ncz,ncz_s,ncz_e
       WRITE(kdens2,'(6e12.5)') a,b,c,alf,bet,gamm
       WRITE(kdens2,'(''ZYX'')')
       DO k=1,ncz
          WRITE(kdens2,'(i8)') k-1
          WRITE(kdens2,'(6e12.5)') ((fact2*Fluct(i,j,k),i=1&
               &,ncx),j=1,ncy)
       END DO
       WRITE(kdens2,'(i8)') -9999
       WRITE(kdens2,'(2(e12.4,1x))') avg,sigma
    END SELECT
    If(dens_histo) Then
       Rewind(kdenhisto_cos)
       Do u=1,2
          Write(kdenhisto_cos,'(''# gofr histogram at step = '',i8,2x,f12.3)')&
               & u,fstep
          Do n=1,Dens_histo_size
             r=Dens_histo_bin*Dble(n-1)
             
             sum=0.0D0
             Do k=1,Dens_histo_size_cos
                sum=sum+Dens_Histo_gofr_cos(k,n,u)
             End Do
             sum=sum*Dens_histo_bin_cos

             Do k=1,Dens_histo_size_cos
                rc=Dens_histo_bin_cos*Dble(k-1)-1.0D0
                If(sum /= 0.0D0) Then
                   y=Dens_Histo_gofr_cos(k,n,u)/sum
                Else
                   y=0.0D0
                End If
                Write(kdenhisto_cos,'(2f14.4,e15.8)') r,rc,y
             End Do
          End Do
       End Do

       If(NoHydrogen) Then
          Write(kdenhisto,'(''# No Hydrogens gofr histogram at step = '',f12.3)') fstep
       Else
          Write(kdenhisto,'(''# gofr histogram at step = '',f12.3)') fstep
       End If
       Do n=1,Dens_histo_size
          r=Dens_histo_bin*Dble(n-1)
          y=(Dens_histo_svol(n,0)+Dens_histo_svol(n,1)+Dens_histo_svol(n,2))/Dble(Counter)
          z0=Dens_histo_gofr(n)/Dble(Counter)
          z=(Dens_histo_gofr_ac(n,0)+Dens_histo_gofr_ac(n,1)+Dens_histo_gofr_ac(n,2))/Dble(Counter)
          If(y /= 0.0D0) Then
             Write(kdenhisto,'(f12.4,4e16.8)') r,y,z,z0,z/y
          Else
             Write(kdenhisto,'(f12.4,4e16.8)') r,y,z,0.0D0,0.0D0
          End If
       End Do
       Write(kdenhisto,'(''&'')')
       Write(kdenhisto,'(''# gofr histogram at step = '',f12.3,''Type 0'')') fstep
       Do n=1,Dens_histo_size
          r=Dens_histo_bin*Dble(n-1)
          y=Dens_histo_svol(n,0)/Dble(Counter)
          z=Dens_histo_gofr_ac(n,0)/Dble(Counter)
          If(y /= 0.0D0) Then
             Write(kdenhisto,'(f12.4,3e16.8)') r,y,z,z/y
          Else
             Write(kdenhisto,'(f12.4,3e16.8)') r,y,z,0.0D0
          End If
       End Do
       Write(kdenhisto,'(''&'')')
       Write(kdenhisto,'(''# gofr histogram at step = '',f12.3,''Type 1'')') fstep
       Do n=1,Dens_histo_size
          r=Dens_histo_bin*Dble(n-1)
          y=Dens_histo_svol(n,1)/Dble(Counter)
          z=Dens_histo_gofr_ac(n,1)/Dble(Counter)
          If(y /= 0.0D0) Then
             Write(kdenhisto,'(f12.4,3e16.8)') r,y,z,z/y
          Else
             Write(kdenhisto,'(f12.4,3e16.8)') r,y,z,0.0D0
          End If
       End Do
       Write(kdenhisto,'(''&'')')
       Write(kdenhisto,'(''# gofr histogram at step = '',f12.3,''Type 2'')') fstep
       Do n=1,Dens_histo_size
          r=Dens_histo_bin*Dble(n-1)
          y=Dens_histo_svol(n,2)/Dble(Counter)
          z=Dens_histo_gofr_ac(n,2)/Dble(Counter)
          If(y /= 0.0D0) Then
             Write(kdenhisto,'(f12.4,3e16.8)') r,y,z,z/y
          Else
             Write(kdenhisto,'(f12.4,3e16.8)') r,y,z,0.0D0
          End If
       End Do
       Write(kdenhisto,'(''&'')')
       Write(kdenhisto,'(''# gofr histogram at step = '',f12.3,''Type 1+2'')') fstep
       Do n=1,Dens_histo_size
          r=Dens_histo_bin*Dble(n-1)
          y=(Dens_histo_svol(n,1)+Dens_histo_svol(n,2))/Dble(Counter)
          z=(Dens_histo_gofr_ac(n,1)+Dens_histo_gofr_ac(n,2))/Dble(Counter)
          If(y /= 0.0D0) Then
             Write(kdenhisto,'(f12.4,3e16.8)') r,y,z,z/y
          Else
             Write(kdenhisto,'(f12.4,3e16.8)') r,y,z,0.0D0
          End If
       End Do
    End If
    IF(kpdb /= 0) THEN
       na=0
       count=0
       DO n=1,natoms_Tot
          m=atoms(na+1)
          DO ii=1,m
             i=atoms(na+1+ii)
             count=count+1
             coords(count) % x = xc(i)/DBLE(counter)
             coords(count) % y = yc(i)/DBLE(counter)
             coords(count) % z = zc(i)/DBLE(counter)
          END DO
          na=na+1+m
       END DO
       CLOSE(kpdb)
       OPEN(unit=kpdb,file=filename_pdb,form='FORMATTED',status='OLD')
       IF(PDB_) THEN
          ntap=SIZE(Coords)
          CALL PDB_Out(fstep,kpdb, beta, Coords(:) % x,  Coords(:) %&
               & y, Coords(:) % z, ntap, Coords(:) % res) 
       ELSE
          CALL Write_PDB
       END IF
    END IF
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

    CONTAINS
      SUBROUTINE Write_PDB
        IMPLICIT NONE
        INTEGER :: k,n,m
        REAL(8) :: x1,y1,z1,charge
        CHARACTER(5) :: bet,bet2,rsd

        
        charge=0.0D0
        WRITE(kpdb,2) fstep
        k=0
        DO n=1,SIZE(Coords)
           x1=Coords(n) % x
           y1=Coords(n) % y
           z1=Coords(n) % z
           m=Coords(n) % res
           bet=ADJUSTL(Coords (n) % beta)
           CALL low_up(bet,5)
           rsd=ADJUSTL(Coords(n) % symb)
           CALL low_up(rsd,5)
           WRITE(kpdb,1)'ATOM  ',n,bet,rsd,m,x1,y1,z1,charge,DFLOAT(k)
       END DO
       WRITE(kpdb,'(a)')'TER  '

1      FORMAT(a5,i6,1x,a5,a3,1x,i5,4x,3f8.3,2f6.2)
2      FORMAT('COMMENT 1 Configuration at time step ',f11.2,'      '&
            &,'                 ')

      END SUBROUTINE Write_PDB
      SUBROUTINE Statistics(fact,Density,avg,sigma)
        IMPLICIT NONE 
        REAL(8), DIMENSION (:,:,:) :: Density
        REAL(8) :: avg, sigma,fact
        INTEGER :: i,j,k
        REAL(8) :: avg2,aux
        avg=0.0D0
        avg2=0.0D0
        DO k=1,ncz
           DO j=1,ncy
              DO i=1,ncx
                 aux=fact*Density(i,j,k)
                 avg=avg+aux
                 avg2=avg2+aux**2
              END DO
           END DO
        END DO

        avg=avg/DBLE(ncx*ncy*ncz)
        avg2=avg2/DBLE(ncx*ncy*ncz)
        sigma=avg2-avg**2
      END SUBROUTINE Statistics
  END SUBROUTINE Write_it
  INCLUDE 'DENSITY_Read.f90'
END MODULE DENSITY_Mod
