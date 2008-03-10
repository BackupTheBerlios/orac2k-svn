MODULE DENSITY_Mod

!!$***********************************************************************
!!$   Time-stamp: <2007-12-03 16:40:58 marchi>                           *
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

  USE PBC_Mod
  USE INPUT_Mod, ONLY: Read_String, Parser, err_open,err_end,err_unr&
       &,err_fnf,err_args
  USE PDB, ONLY: PDB_out=>Write_it, PDB_
  IMPLICIT None

  PRIVATE
  PUBLIC Initialize, Initialize_Ar, Initialize_Par, Compute, Write_it&
       &, Read_It, Density_Calc, n_write
  REAL(8), DIMENSION (:,:,:), ALLOCATABLE, SAVE :: Density,Dens2,Loc_Dens
  REAL(8), DIMENSION (:,:,:), ALLOCATABLE, SAVE :: Ndens
  REAL(8), DIMENSION (:), ALLOCATABLE, SAVE :: atmass
  REAL(8), DIMENSION (:), ALLOCATABLE, SAVE :: xc,yc,zc
  REAL(8), SAVE :: a=0.0D0,b=0.0D0,c=0.0D0,alph,bet,gamm,co(3,3),oc(3&
       &,3),volume,xcm,ycm,zcm,Dvolume,rho=0.0D0,unitfluc,UnitDens,Total_Mass
  INTEGER, SAVE :: natom,n_write=1,ncx=0,ncy=0,ncz=0,ntot=1,nats=0&
       &,natoms_Tot=0
  INTEGER, SAVE :: kdensity=0,kdens2=0,kpdb=0
  INTEGER, SAVE :: counter=0,node=0
  LOGICAL, SAVE :: Density_Calc=.FALSE.,Density_Avg=.FALSE.
  REAL(8), SAVE :: x_ad=0.0D0,y_ad=0.0D0,z_ad=0.0D0
  LOGICAL, DIMENSION (:), ALLOCATABLE, SAVE :: mask
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
  LOGICAL, SAVE :: Other_Density=.FALSE.
CONTAINS
!!$---- This subroutine is part of the program ORAC ----*

  SUBROUTINE Initialize(nodea,prsymb,betab,res1,res2,mass,ntap)

!!$======================== DECLARATIONS ================================*

    IMPLICIT none
    INCLUDE 'unit.h'

!!$----------------------------- ARGUMENTS ------------------------------*

    CHARACTER(8) :: prsymb(*)
    CHARACTER(7) :: betab(*)
    INTEGER :: res1(*),res2(*),ntap,nodea
    REAL(8) :: mass(*)

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: i,na,count,ii,m,n
    LOGICAL :: ok
    REAL(8) :: kT

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    node=nodea
    kT=300.0D0*gascon/1000.0D0
    unitfluc=efact/unitp/kT
    unitdens=unitm/unitl**3/1000.0D0
    co=0.0D0
    oc=0.0D0

    co(1,1)=a/2.0D0
    co(2,2)=b/2.0D0
    co(3,3)=c/2.0D0

    oc(1,1)=1.0D0/co(1,1)
    oc(2,2)=1.0D0/co(2,2)
    oc(3,3)=1.0D0/co(3,3)

    volume=a*b*c

    natom=ntap

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
       count=0
       DO n=1,natoms_Tot
          m=atoms(na+1)
          count=count+m
          na=na+1+m
       END DO
       ALLOCATE(coords(count),xc(ntap),yc(ntap),zc(ntap))
       xc=0.0D0
       yc=0.0D0
       zc=0.0D0
       na=0
       count=0
       DO n=1,natoms_Tot
          m=atoms(na+1)
          DO ii=1,m
             i=atoms(na+1+ii)
             count=count+1
             coords(count) % beta = TRIM(betab(i))
             coords(count) % symb = TRIM(prsymb(res2(i)))
             coords(count) % res = res1(i)
          END DO
          na=na+1+m
       END DO
    END IF
    ALLOCATE(beta(SIZE(Coords)))
    beta=ADJUSTL(Coords(:) % beta)

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
  SUBROUTINE Compute(xp0,yp0,zp0,Real_Volume,My_Density)

!!$======================== DECLARATIONS ================================*

    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

    REAL(8) :: xp0(*),yp0(*),zp0(*),Real_Volume
    REAL(8), OPTIONAL :: My_Density(*)

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: nx,ny,nz,na,m,ii,i,qq,j,k
    REAL(8) :: x1,y1,z1
    INTEGER :: n,nc2x,nc2y,nc2z,p,l
    REAL(8) :: xpa,ypa,zpa,aux1,aux2,fract

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
          xpa=oc(1,1)*xp0(n)+oc(1,2)*yp0(n)+oc(1,3)*zp0(n)
          ypa=oc(2,1)*xp0(n)+oc(2,2)*yp0(n)+oc(2,3)*zp0(n)
          zpa=oc(3,1)*xp0(n)+oc(3,2)*yp0(n)+oc(3,3)*zp0(n)
          x1=xpa*0.5D0
          y1=ypa*0.5D0
          z1=zpa*0.5D0
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
       DO k=1,ncz
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

  END SUBROUTINE Compute
  SUBROUTINE Write_it(fstep)

!!$======================== DECLARATIONS ================================*

    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

    REAL(8) :: fstep

!!$------------------------- LOCAL VARIABLES ----------------------------*

    REAL(8) :: Dvolume,avogad=6.0225d23,unitcm=1.0D-8,fact,coa(3,3)&
         &,dummy=0.0D0,a0=0.529177249D0,a,b,c,alf,bet,gamm,avg,sigma&
         &,fact2,Mass_Avg

    INTEGER :: one=1,i,j,k,count,ncx2,ncx_s,ncx_e,ncy2,ncy_s,ncy_e&
         &,ncz2,ncz_s,ncz_e,na,n,m,ii,ntap
    REAL(8), DIMENSION (:,:,:), ALLOCATABLE, SAVE :: Fluct

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    IF(.NOT. ALLOCATED(Fluct)) ALLOCATE(fluct(ncx,ncy,ncz))
    coa=co
    coa(1,1)=(2.0D0/a0)*coa(1,1)/DBLE(ncx)
    coa(1,2)=(2.0D0/a0)*coa(1,2)/DBLE(ncy)
    coa(2,2)=(2.0D0/a0)*coa(2,2)/DBLE(ncy)
    coa(1,3)=(2.0D0/a0)*coa(1,3)/DBLE(ncz)
    coa(2,3)=(2.0D0/a0)*coa(2,3)/DBLE(ncz)
    coa(3,3)=(2.0D0/a0)*coa(3,3)/DBLE(ncz)
    xcm=-(co(1,1)+co(1,2)+co(1,3))/a0
    ycm=-(co(2,1)+co(2,2)+co(2,3))/a0
    zcm=-(co(3,1)+co(3,2)+co(3,3))/a0

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
       
       WRITE(kdensity,'(i5,3f12.6)') ncx,coa(1,1),coa(2,1),coa(3,1)
       WRITE(kdensity,'(i5,3f12.6)') ncy,coa(1,2),coa(2,2),coa(3,2)
       WRITE(kdensity,'(i5,3f12.6)') ncz,coa(1,3),coa(2,3),coa(3,3)
       
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
       
       WRITE(kdens2,'(i5,3f12.6)') ncx,coa(1,1),coa(2,1),coa(3,1)
       WRITE(kdens2,'(i5,3f12.6)') ncy,coa(1,2),coa(2,2),coa(3,2)
       WRITE(kdens2,'(i5,3f12.6)') ncz,coa(1,3),coa(2,3),coa(3,3)
       
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

       CALL rotb(a,b,c,alf,bet,gamm,co)
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
