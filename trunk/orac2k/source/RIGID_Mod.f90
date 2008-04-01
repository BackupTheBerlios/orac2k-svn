MODULE RIGID_Mod

!!$***********************************************************************
!!$   Time-stamp: <2008-03-28 11:07:09 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Sat Oct 15 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

  LOGICAL, SAVE :: Rigid_object=.FALSE.
  CHARACTER(7), SAVE :: target
  CHARACTER(3), DIMENSION(22), SAVE :: protein=(/ 'ala','arg','asn'&
       &,'asp','cys','gln','glu','gly','hsd','hse','hsp','ile'&
       &,'leu','lys','met','phe','pro','ser','thr','trp','tyr'&
       &,'val'/)
  LOGICAL, DIMENSION(:), ALLOCATABLE, SAVE :: mask
  REAL(8), ALLOCATABLE, SAVE :: wca(:)
  REAL(8), ALLOCATABLE, SAVE :: xp_ini(:),yp_ini(:),zp_ini(:)
  INTEGER, SAVE :: nato_slt,ntap
CONTAINS
  SUBROUTINE Initialize(prsymb,beta,res,ss_index,ntapa)

!!$======================== DECLARATIONS ================================*

    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*
    
    INTEGER :: ntapa,ss_index(*),res(*)
    CHARACTER(8) :: prsymb(*)
    CHARACTER(7) :: beta(*)

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: i,j
    CHARACTER(3) :: Residue

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    ntap=ntapa
    ALLOCATE(mask(ntap))
    mask=.FALSE.
    SELECT CASE (target)
    CASE('PROTEIN')
       DO i=1,ntap
          Residue=prsymb(res(i))(1:3)
          IF(beta(i)(1:1) == 'h') CYCLE
          DO j=1,LEN(protein)
             IF(Residue == protein(j)) mask(i)=.TRUE.
             IF(Residue == upper(protein(j))) mask(i)=.TRUE.
          END DO
       END DO
    CASE DEFAULT
       DO i=1,ntap
          IF(beta(i)(1:1) == 'h') CYCLE
          IF(ss_index(i) == 1) mask(i)=.TRUE.
       END DO
    END SELECT

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

  END SUBROUTINE Initialize
  SUBROUTINE Initialize_xp(xp,yp,zp,nato_slta)
    INTEGER :: nato_slta
    REAL(8) :: xp(*),yp(*),zp(*)
    REAL(8) :: norm
    INTEGER :: i

    nato_slt=nato_slta
    ALLOCATE(xp_ini(nato_slt),yp_ini(nato_slt),zp_ini(nato_slt))
    xp_ini=xp(1:nato_slt)
    yp_ini=yp(1:nato_slt)
    zp_ini=zp(1:nato_slt)
    ALLOCATE(wca(nato_slt))
    DO i=1,nato_slt
       IF(mask(i)) THEN
          wca(i)=1.0D0
       END IF
    END DO
    norm=SUM(wca)
    wca=wca/norm
  END SUBROUTINE Initialize_xp
  SUBROUTINE Compute(xp0,yp0,zp0)

!!$======================== DECLARATIONS ================================*

    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*
    
    REAL(8) :: xp0(*),yp0(*),zp0(*)

!!$------------------------- LOCAL VARIABLES ----------------------------*

    REAL(8) :: xcm,ycm,zcm,totmass
    INTEGER :: i,n,m
    REAL(8), POINTER :: xyz(:,:),xyz0(:,:),xyzfit(:,:)
    REAL(8) :: RotMat(3,3),q(0:7),xc,yc,zc,xd,yd,zd,xt_cm,yt_cm,zt_cm&
         &,err,dcm(3),xpp,ypp,zpp,xpcc,ypcc,zpcc
    INTEGER :: iret
    CHARACTER(80) :: msg
    INTEGER, SAVE :: last0=0

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    
    last0=last0+1
    ALLOCATE(xyz(3,nato_slt),xyz0(3,nato_slt),xyzfit(3,nato_slt))
    m=nato_slt
    xyz=0.0D0
    xyz0=0.0D0
    xyzfit=0.0D0    
    DO n=1,m
       xyz(1,n)=xp0(n)
       xyz(2,n)=yp0(n)
       xyz(3,n)=zp0(n)
       xyz0(1,n)=xp_ini(n)
       xyz0(2,n)=yp_ini(n)
       xyz0(3,n)=zp_ini(n)
    END DO
    CALL rigfit(0,m,xyz0,xyz,wca,wca,q,dcm,xyzfit,err,iret,msg)
    RotMat(1,1)=-2.d0*q(2)**2-2.d0*q(3)**2+1.d0
    RotMat(1,2)=2.d0*(-q(0)*q(3)+q(1)*q(2))
    RotMat(1,3)=2.d0*(q(0)*q(2)+q(1)*q(3))
    RotMat(2,1)=2.d0*(q(0)*q(3)+q(1)*q(2))
    RotMat(2,2)=-2.d0*q(1)**2-2.d0*q(3)**2+1.d0
    RotMat(2,3)=2.d0*(-q(0)*q(1)+q(2)*q(3))
    RotMat(3,1)=2.d0*(-q(0)*q(2)+q(1)*q(3))
    RotMat(3,2)=2.d0*(q(0)*q(1)+q(2)*q(3))
    RotMat(3,3)=-2.d0*q(1)**2-2.d0*q(2)**2+1.d0
    xcm=0.0D0
    ycm=0.0D0
    zcm=0.0D0
    DO i=1,m
       xcm=xcm+wca(i)*xyz0(1,i)
       ycm=ycm+wca(i)*xyz0(2,i)
       zcm=zcm+wca(i)*xyz0(3,i)
    END DO

    xt_cm=xcm
    yt_cm=ycm
    zt_cm=zcm
    xpcc=xt_cm-dcm(1)
    ypcc=yt_cm-dcm(2)
    zpcc=zt_cm-dcm(3)

    DO i=1,ntap
       xd=xp0(i)-xpcc
       yd=yp0(i)-ypcc
       zd=zp0(i)-zpcc
       xc=RotMat(1,1)*xd+RotMat(1,2)*yd+RotMat(1,3)*zd
       yc=RotMat(2,1)*xd+RotMat(2,2)*yd+RotMat(2,3)*zd
       zc=RotMat(3,1)*xd+RotMat(3,2)*yd+RotMat(3,3)*zd
       xp0(i)=xc+xt_cm
       yp0(i)=yc+yt_cm
       zp0(i)=zc+zt_cm
    END DO
    DEALLOCATE(xyz,xyz0,xyzfit)
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

  END SUBROUTINE Compute
  FUNCTION upper(string) RESULT(out)
    IMPLICIT NONE 
    CHARACTER(3) :: string,out
    INTEGER :: upper_to_lower,len_string,i

    upper_to_lower = IACHAR("a") - IACHAR("A")
    out=string
    len_string=LEN_TRIM(out)
    DO i=1,len_string
       SELECT CASE (out(i:i))
       CASE ('a':'z')
          out(i:i) = ACHAR(IACHAR(out(i:i)) - upper_to_lower)
       END SELECT
    END DO
  END FUNCTION upper
END MODULE RIGID_Mod
