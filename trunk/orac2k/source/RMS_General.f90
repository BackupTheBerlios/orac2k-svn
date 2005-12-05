MODULE RMS_General

!!$***********************************************************************
!!$   Time-stamp: <2005-10-24 20:31:24 marchi>                           *
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

  TYPE rms
     CHARACTER(8) :: target
     CHARACTER(7), DIMENSION(:), POINTER :: beta
     INTEGER :: n0
  END TYPE rms
  TYPE rms_units
     INTEGER :: n0
     CHARACTER(8) :: target
     INTEGER, DIMENSION(:), POINTER :: ind
  END TYPE rms_units

!!$========================================================================
!!$-- Statics 
!!$========================================================================

  LOGICAL, SAVE :: General=.FALSE.,matrix=.FALSE.
  INTEGER, SAVE :: n_write,krms_sub
  CHARACTER(80), SAVE :: filename
  INTEGER, SAVE :: N_Calls=0,natom,nbunits

  INTEGER :: n_Res_u,n_Res
  PARAMETER(n_Res_u=30,n_Res=100)
  INTEGER, SAVE :: units
  TYPE(rms), DIMENSION(:), ALLOCATABLE, SAVE :: Res_u
  TYPE(rms_units), DIMENSION(:), ALLOCATABLE, SAVE :: atoms

!!$========================================================================
!!$-- Allocatables
!!$========================================================================

  REAL(8), DIMENSION(:), ALLOCATABLE, SAVE :: xpt,ypt,zpt
  REAL(8), DIMENSION(:), ALLOCATABLE :: xpr,ypr,zpr
  REAL(8), DIMENSION(:), ALLOCATABLE, SAVE :: xp_avg,yp_avg,zp_avg
  REAL(8), DIMENSION(:), ALLOCATABLE, SAVE :: xp_avg2,yp_avg2,zp_avg2
  REAL(8), DIMENSION(:), ALLOCATABLE, SAVE :: wc
  REAL(8), DIMENSION(:,:), ALLOCATABLE, SAVE :: rms
CONTAINS

!!$======================== DECLARATIONS ================================*

  SUBROUTINE Initialize(atres,nbun_slt,res,prsymb,beta,nat_slt)
    
    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

    INTEGER :: atres(2,*),nbun_slt,res(*),nat_slt
    CHARACTER(8) :: prsymb(*)
    CHARACTER(7) :: beta(*)

!!$----------------------- VARIABLES IN COMMON --------------------------*

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: i,k1,k2,n,nn
    INTEGER :: count,count2
    INTEGER, DIMENSION (:), ALLOCATABLE :: index
    LOGICAL, DIMENSION(:), ALLOCATABLE :: mask
    REAL(8) :: summa

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    nbunits=nbun_slt
    natom=nat_slt
    ALLOCATE(atoms(nbunits),index(nbunits),mask(natom))
    ALLOCATE(xpt(nbunits),ypt(nbunits),zpt(nbunits))
    ALLOCATE(xpr(nbunits),ypr(nbunits),zpr(nbunits))
    ALLOCATE(xp_avg(natom),yp_avg(natom),zp_avg(natom))
    ALLOCATE(xp_avg2(ntap),yp_avg2(ntap),zp_avg2(ntap))
    ALLOCATE(wc(natom))
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

    DO n=1,natom
       wc(n)=0.0D0
       IF(beta(n)(1:1) /= 'h') THEN
          wc(n)=1.0D0
       END IF
    END IF
    summa=SUM(wc)
    wc=wc/summa

    DO k1=1,units
       count=0
       index=0
       DO n=1,nbunits
          i=atres(1,n)
          IF(Res_u(k1)%target(1:3) == 'all') THEN
             count=count+1
             index(count)=n
          ELSE IF(prsymb(res(i))(1:3) == Res_u(k1)%target(1:3)) THEN
             count=count+1
             index(count)=n
          END IF
       END DO
       DO nn=1,count
          n=index(nn)             
          mask(atres(1,n):atres(2,n))=.FALSE.
          DO i=atres(1,n),atres(2,n)
             DO k2=1,Res_u(k1)%n0
                IF(Res_u(k1)%beta(k2)(1:3) == 'all') THEN
                   mask(i)=.TRUE.
                ELSE IF(Res_u(k1)%beta(k2) == beta(i)) THEN
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
    count=SUM(atoms(1:nbunits)%n0)
    ALLOCATE(rms(count,count))
    rms=0.0D0

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

  END SUBROUTINE Initialize
  SUBROUTINE Compute(xp0,yp0,zp0)
    IMPLICIT NONE 

!!$----------------------------- ARGUMENTS ------------------------------*

    REAL(8) :: xp0(*),yp0(*),zp0(*)

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: i,n,nn
    INTEGER :: count,c1,c2
    REAL(8) :: xc,yc,zc,xci,yci,zci,xcj,ycj,zcj

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    N_Calls=N_Calls+1
    IF(matrix) THEN
       c1=0
       DO n=1,nbunits
          DO nn=1,atoms(n)%n0
             i=atoms(n)%ind(nn)
             xci=xp0(i)
             yci=yp0(i)
             zci=zp0(i)
             xp_avg(i)=xp_avg(i)+xci
             yp_avg(i)=yp_avg(i)+yci
             zp_avg(i)=zp_avg(i)+zci
             c1=c1+1
             c2=0
             DO m=1,nbunits
                DO mm=1,atoms(m)%n0
                   j=atoms(m)%ind(mm)
                   xcj=xp0(j)
                   ycj=yp0(j)
                   zcj=zp0(j)
                   c2=c2+1
                   rms(c1,c2)=rms(c1,c2) + (xci*xcj + yci*ycj + zci*zcj)
                END DO
             END DO
          END DO
       END DO
    ELSE
       DO n=1,nbunits
          IF(atoms(n)%n0 /= 0) THEN
             DO nn=1,atoms(n)%n0
                i=atoms(n)%ind(nn)
                xc=xp0(i)
                yc=yp0(i)
                zc=zp0(i)
                xp_avg(i)=xp_avg(i)+xc
                yp_avg(i)=yp_avg(i)+yc
                zp_avg(i)=zp_avg(i)+zc
                xp_avg2(i)=xp_avg2(i)+xc**2
                yp_avg2(i)=yp_avg2(i)+yc**2
                zp_avg2(i)=zp_avg2(i)+zc**2
             END DO
          END IF
       END DO
    END IF

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
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: rms_avg

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    count=LEN(rms,1)
    ALLOCATE(rms_avg(count,count))
    ALLOCATE(xxii(count))


    REWIND(krms_sub)
    WRITE(krms_sub,'(''Ring RMS at '',f12.3,'' fs'')') fstep 
    fcounter=1.0D0/DBLE(N_Calls)
    IF(matrix) THEN
       WRITE(krms_sub,'(''Averages After '',i7,'' Calls '')') N_Calls
       c1=0
       DO n=1,nbunits
          DO nn=1,atoms(n)%n0
             i=atoms(n)%ind(nn)
             xci=xp_avg(i)*fcounter
             yci=yp_avg(i)*fcounter
             zci=zp_avg(i)*fcounter
             c1=c1+1
             xixj=rms(c1,c1)*fcounter
             xxii(c1)=SQRT(xixj-(xci*xci+yci*yci+zci*zci))
          END DO
       END DO
       c1=0
       DO n=1,nbunits
          DO nn=1,atoms(n)%n0
             i=atoms(n)%ind(nn)
             xci=xp_avg(i)*fcounter
             yci=yp_avg(i)*fcounter
             zci=zp_avg(i)*fcounter
             c1=c1+1
             c2=0
             DO m=1,nbunits
                DO mm=1,atoms(m)%n0
                   j=atoms(m)%ind(mm)
                   xcj=xp_avg(j)*fcounter
                   ycj=yp_avg(j)*fcounter
                   zcj=zp_avg(j)*fcounter
                   c2=c2+1
                   avg=rms(c1,c2)*fcounter
                   rms_avg(c1,c2)=(xci*xcj+yci*ycj+zci*zcj) + avg
                   WRITE(krms_sub,'(2i5,e20.12)') i,j,rms_avg(c1,c2)&
                        &/xxii(c1)/xxii(c2)
                END DO
             END DO
          END DO
       END DO
       IF(N_Calls /= 1) THEN
          lda=LEN(rms,1)
          CALL DPOTRF('U',lda,rms,lda,info)
          CALL DPOTRI('U',lda,rms,lda,info)
       END IF
       WRITE(krms_matrix,'(''Inverse of RMS Matrix'')')
       c1=0
       DO n=1,nbunits
          summb=0.0D0
          DO nn=1,atoms(n)%n0
             i=atoms(n)%ind(nn)
             c1=c1+1
             summb=summb+rms(c1,c1)
          END DO
          WRITE(krms_sub,'(a7,e15.6)') atoms(n)%target,summb/DBLE(atoms(n)%n0)
          DO nn=1,atoms(n)%n0
             i=atoms(n)%ind(nn)
             c1=c1+1
             WRITE(krms_sub,'(2i5,e20.12)') i,i,rms(c1,c1)
          END DO
       END DO
    ELSE
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
    END IF
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

  END SUBROUTINE Write_it
  SUBROUTINE Rotate(xp0,yp0,zp0)
    IMPLICIT NONE
    REAL(8) :: xp0(*),yp0(*),zp0(*)

    REAL(8), DIMENSION (:,:), ALLOCATABLE :: xyz,xyzfit
    REAL(8), DIMENSION (:), ALLOCATABLE :: work
    
    REAL(8) :: qt(4),err,dr1,dr2
    INTEGER :: i,ii

    ALLOCATE(xyz(3,natom),xyzfit(3,natom),work(natom+1))
    
    DO i=1,natom 
       xyz(1,i)=xp0(i)
       xyz(2,i)=yp0(i)
       xyz(3,i)=zp0(i)
    END DO
    CALL xfit(xyz0,xyz,xyzfit,qt,wc,work,natom,err)
    DO i=1,natom
       xp0(i)=xyzfit(1,i)
       yp0(i)=xyzfit(2,i)
       zp0(i)=xyzfit(3,i)
    END DO
    DEALLOCATE(xyz,xyzfit,work)
  END SUBROUTINE Rotate
  INCLUDE 'RMS_Genaral_Read.f90'
END MODULE RMS_General
