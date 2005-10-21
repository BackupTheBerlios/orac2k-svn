MODULE RMS_Subtract_Mod

!!$***********************************************************************
!!$   Time-stamp: <2005-10-21 17:33:32 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Oct 20 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*
  USE INPUT_Mod, ONLY: Read_String, Parser, err_open,err_end,err_unr&
       &,err_fnf,err_args

  TYPE subtract
     CHARACTER(8) :: target
     CHARACTER(7), DIMENSION(:), POINTER :: beta
     INTEGER :: n0
  END TYPE subtract
  TYPE rms_units
     INTEGER :: n0
     CHARACTER(8) :: target
     INTEGER, DIMENSION(:), POINTER :: ind
  END TYPE rms_units

!!$========================================================================
!!$-- Statics 
!!$========================================================================

  LOGICAL, SAVE :: RMS_Subtract=.FALSE.,RMS_test=.FALSE.
  INTEGER, SAVE :: n_write,krms_sub
  CHARACTER(80), SAVE :: filename
  INTEGER, SAVE :: counter=0,natom,nbunits

  INTEGER :: n_Res_u,n_Res
  PARAMETER(n_Res_u=25,n_Res=100)
  INTEGER, SAVE :: units
  TYPE(Subtract), DIMENSION(:), ALLOCATABLE, SAVE :: Res_u
  TYPE(rms_units), DIMENSION(:), ALLOCATABLE, SAVE :: atoms

!!$========================================================================
!!$-- Allocatables
!!$========================================================================

  REAL(8), DIMENSION(:), ALLOCATABLE, SAVE :: xpt,ypt,zpt
  REAL(8), DIMENSION(:), ALLOCATABLE :: xpr,ypr,zpr
  REAL(8), DIMENSION(:), ALLOCATABLE, SAVE :: xp_avg,yp_avg,zp_avg
  REAL(8), DIMENSION(:), ALLOCATABLE, SAVE :: xp_avg2,yp_avg2,zp_avg2
CONTAINS

!!$======================== DECLARATIONS ================================*

  SUBROUTINE Initialize(atres,nbun,res,prsymb,beta,ntap)
    
    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

    INTEGER :: atres(2,*),nbun,res(*),ntap
    CHARACTER(8) :: prsymb(*)
    CHARACTER(7) :: beta(*)

!!$----------------------- VARIABLES IN COMMON --------------------------*

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: i,k1,k2,n,nn
    INTEGER :: count,count2
    INTEGER, DIMENSION (:), ALLOCATABLE :: index
    LOGICAL, DIMENSION(:), ALLOCATABLE :: mask


!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    nbunits=nbun
    natom=ntap
    ALLOCATE(atoms(nbun),index(nbun),mask(ntap))
    ALLOCATE(xpt(nbun),ypt(nbun),zpt(nbun))
    ALLOCATE(xpr(nbun),ypr(nbun),zpr(nbun))
    ALLOCATE(xp_avg(ntap),yp_avg(ntap),zp_avg(ntap))
    ALLOCATE(xp_avg2(ntap),yp_avg2(ntap),zp_avg2(ntap))
    xpt=0.0D0
    ypt=0.0D0
    zpt=0.0D0
    xpr=0.0D0
    ypr=0.0D0
    zpr=0.0D0
    xp_avg=0.0D0
    yp_avg=0.0D0
    zp_avg=0.0D0
    xp_avg2=0.0D0
    yp_avg2=0.0D0
    zp_avg2=0.0D0
    DO k1=1,units
       count=0
       index=0
       DO n=1,nbun
          i=atres(1,n)
          IF(prsymb(res(i))(1:3) == Res_u(k1)%target(1:3)) THEN
             count=count+1
             index(count)=n
          END IF
       END DO
       DO nn=1,count
          n=index(nn)             
          mask(atres(1,n):atres(2,n))=.FALSE.
          DO i=atres(1,n),atres(2,n)
             DO k2=1,Res_u(k1)%n0
                IF(Res_u(k1)%beta(k2) == beta(i)) THEN
                   mask(i)=.TRUE.
                END IF
             END DO
          END DO
          count2=0
          DO i=atres(1,n),atres(2,n)
             IF(mask(i)) count2=count2+1
          END DO
          atoms(n)%n0=count2
          atoms(n)%target=prsymb(res(atres(1,n)))
          ALLOCATE(atoms(n)%ind(count2))
          count2=0
          DO i=atres(1,n),atres(2,n)
             IF(mask(i)) THEN
                count2=count2+1
                atoms(n)%ind(count2)=i
             END IF
          END DO
       END DO
    END DO
    DEALLOCATE(index,mask)

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

  END SUBROUTINE Initialize
  SUBROUTINE Compute(xp0,yp0,zp0)
    IMPLICIT NONE 

!!$----------------------------- ARGUMENTS ------------------------------*

    REAL(8) :: xp0(*),yp0(*),zp0(*)

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: i,n,nn
    INTEGER :: count
    REAL(8) :: xc,yc,zc

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    counter=counter+1
    IF(counter == 1) THEN
       DO n=1,nbunits
          count=0
          IF(atoms(n)%n0 /= 0) THEN
             DO nn=1,atoms(n)%n0
                count=count+1
                i=atoms(n)%ind(nn)
                xpt(n)=xpt(n)+xp0(i)
                ypt(n)=ypt(n)+yp0(i)
                zpt(n)=zpt(n)+zp0(i)
             END DO
             xpt(n)=xpt(n)/DBLE(count)
             ypt(n)=ypt(n)/DBLE(count)
             zpt(n)=zpt(n)/DBLE(count)
          END IF
       END DO
    END IF
    DO n=1,nbunits
       count=0
       xpr(n)=0.0D0
       ypr(n)=0.0D0
       zpr(n)=0.0D0
       IF(atoms(n)%n0 /= 0) THEN
          DO nn=1,atoms(n)%n0
             count=count+1
             i=atoms(n)%ind(nn)
             xpr(n)=xpr(n)+xp0(i)
             ypr(n)=ypr(n)+yp0(i)
             zpr(n)=zpr(n)+zp0(i)
          END DO
          xpr(n)=xpr(n)/DBLE(count)
          ypr(n)=ypr(n)/DBLE(count)
          zpr(n)=zpr(n)/DBLE(count)

          DO nn=1,atoms(n)%n0
             i=atoms(n)%ind(nn)
             IF(RMS_test) THEN
                xc=xp0(i)
                yc=yp0(i)
                zc=zp0(i)
             ELSE
                xc=xp0(i)+xpt(n)-xpr(n)
                yc=yp0(i)+ypt(n)-ypr(n)
                zc=zp0(i)+zpt(n)-zpr(n)
             END IF
             xp_avg(i)=xp_avg(i)+xc
             yp_avg(i)=yp_avg(i)+yc
             zp_avg(i)=zp_avg(i)+zc
             xp_avg2(i)=xp_avg2(i)+xc**2
             yp_avg2(i)=yp_avg2(i)+yc**2
             zp_avg2(i)=zp_avg2(i)+zc**2
          END DO
       END IF
    END DO

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

  END SUBROUTINE Compute
  SUBROUTINE Write_it(fstep)
    IMPLICIT NONE 

!!$----------------------------- ARGUMENTS ------------------------------*
    
    REAL(8) :: fstep

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: i,n,nn
    INTEGER :: count
    REAL(8) :: xc,yc,zc,xc2,yc2,zc2,fcounter,rms_tot,rms

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    REWIND(krms_sub)

    WRITE(krms_sub,'(''Ring RMS at '',f12.3,'' fs'')') fstep 
    fcounter=1.0D0/DBLE(counter)
    DO n=1,nbunits
       IF(atoms(n)%n0 /= 0) THEN
          rms_tot=0.0D0
          DO nn=1,atoms(n)%n0
             i=atoms(n)%ind(nn)
             xc=xp_avg(i)*fcounter
             yc=yp_avg(i)*fcounter
             zc=zp_avg(i)*fcounter
             xc2=xp_avg2(i)*fcounter
             yc2=yp_avg2(i)*fcounter
             zc2=zp_avg2(i)*fcounter
             rms_tot=rms_tot+xc2+yc2+zc2-xc**2-yc**2-zc**2
          END DO
          rms_tot=rms_tot/DBLE(atoms(n)%n0)
          WRITE(krms_sub,'(a7,2x,e16.7)') atoms(n)%target,rms_tot
          DO nn=1,atoms(n)%n0
             i=atoms(n)%ind(nn)
             xc=xp_avg(i)*fcounter
             yc=yp_avg(i)*fcounter
             zc=zp_avg(i)*fcounter
             xc2=xp_avg2(i)*fcounter
             yc2=yp_avg2(i)*fcounter
             zc2=zp_avg2(i)*fcounter
             rms=xc2+yc2+zc2-xc**2-yc**2-zc**2
             WRITE(krms_sub,'(i8,2x,e16.7)') i,rms
          END DO
       END IF
    END DO
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

  END SUBROUTINE Write_it
  INCLUDE 'RMS_Subtract_Read.f90'
END MODULE RMS_Subtract_Mod
