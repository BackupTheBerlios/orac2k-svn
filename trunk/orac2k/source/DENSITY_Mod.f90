MODULE DENSITY_Mod

!!$***********************************************************************
!!$   Time-stamp: <2006-04-26 11:42:49 marchi>                           *
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

  USE INPUT_Mod, ONLY: Read_String, Parser, err_open,err_end,err_unr&
       &,err_fnf,err_args
  REAL(8), DIMENSION (:,:,:), ALLOCATABLE, SAVE :: Density
  REAL(8), DIMENSION (:), ALLOCATABLE, SAVE :: atmass
  REAL(8), SAVE :: xcm,ycm,zcm
  INTEGER, SAVE :: natom,n_write=1,ncx=0,ncy=0,ncz=0
  INTEGER, SAVE :: kdensity=0
  INTEGER, SAVE :: counter=0
  LOGICAL, SAVE :: Density_Calc=.FALSE.,Density_Avg=.FALSE.
  LOGICAL, DIMENSION (:), ALLOCATABLE, SAVE :: mask
  CHARACTER(80), SAVE :: filename
  CHARACTER(8), SAVE :: Target_Res,file_format='cube'
CONTAINS
!!$---- This subroutine is part of the program ORAC ----*

  SUBROUTINE Initialize(prsymb,res,mass,ntap)

!!$======================== DECLARATIONS ================================*

    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

    CHARACTER(8) :: prsymb(*)
    INTEGER :: res(*),ntap
    REAL(8) :: mass(*)

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: i

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    natom=ntap
    ALLOCATE(atmass(ntap),mask(ntap))
    mask=.FALSE.

    DO i=1,ntap
       IF(Target_Res == prsymb(res(i))) THEN
          mask(i)=.TRUE.
       END IF
       atmass(i)=mass(i)
    END DO
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*
  END SUBROUTINE Initialize

  SUBROUTINE Initialize_Ar

!!$======================== DECLARATIONS ================================*

    IMPLICIT none

!!$------------------------- LOCAL VARIABLES ----------------------------*


!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    ALLOCATE(Density(ncx,ncy,ncz))
    Density=0.0D0

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*
  END SUBROUTINE Initialize_Ar
  SUBROUTINE Compute(xp0,yp0,zp0,volume)

!!$======================== DECLARATIONS ================================*

    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

    REAL(8) :: xp0(*),yp0(*),zp0(*),volume

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: nx,ny,nz
    REAL(8) :: x1,y1,z1,xc,yc,zc,Dvolume
    INTEGER :: n,nc2x,nc2y,nc2z

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

!!$*=======================================================================
!!$*     Compute chain list for system
!!$*=======================================================================

    Dvolume=Volume/DBLE(ncx*ncy*ncz)
    DO n=1,natom
       IF(Mask(n)) THEN
          x1=xp0(n)*0.5D0
          y1=yp0(n)*0.5D0
          z1=zp0(n)*0.5D0
          nx=INT(DBLE(ncx-1)*(x1-ANINT(x1)+0.5D0)+0.5D0)
          ny=INT(DBLE(ncy-1)*(y1-ANINT(y1)+0.5D0)+0.5D0)
          nz=INT(DBLE(ncz-1)*(z1-ANINT(z1)+0.5D0)+0.5D0)
          Density(nx+1,ny+1,nz+1)=Density(nx+1,ny+1,nz+1)+atmass(n)/Dvolume
       END IF
    END DO
    counter=counter+1

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*
  END SUBROUTINE Compute
  SUBROUTINE Write_it(Volume,co)

!!$======================== DECLARATIONS ================================*

    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

    REAL(8) :: co(3,3),volume

!!$------------------------- LOCAL VARIABLES ----------------------------*

    REAL(8) :: Dvolume,avogad=6.0225d23,unitcm=1.0D-8,fact,coa(3,3)&
         &,dummy=0.0D0,a0=0.529177249D0,a,b,c,alf,bet,gamm,avg,sigma

    INTEGER :: one=1,i,j,k,count,ncx2,ncx_s,ncx_e,ncy2,ncy_s,ncy_e&
         &,ncz2,ncz_s,ncz_e 

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

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
    CLOSE(kdensity)
    OPEN(unit=kdensity,file=filename,form='FORMATTED',status='OLD')


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
       CALL Statistics(avg,sigma)
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
          WRITE(kdensity,'(6e12.5)') ((fact*Density(i,j,k),i=1&
               &,ncx),j=1,ncy)
       END DO
       WRITE(kdensity,'(i8)') -9999
       WRITE(kdensity,'(2(e12.4,1x))') avg,sigma
    END SELECT
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

    CONTAINS
      SUBROUTINE Statistics(avg,sigma)
        IMPLICIT NONE 
        REAL(8) :: avg, sigma
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
