MODULE RMS_Subtract_Mod

!!$***********************************************************************
!!$   Time-stamp: <2006-02-24 12:26:38 marchi2>                           *
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
  USE Xerror_Mod

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

  LOGICAL, SAVE :: RMS_Subtract=.FALSE.,RMS_test=.FALSE.,RMS_Invert&
       &=.TRUE.
  INTEGER, SAVE :: n_write,krms_sub=0,krms_dip=0
  CHARACTER(80), SAVE :: filename,filename_dip,rms_type='cm'
  INTEGER, SAVE :: counter=0,natom,nbunits

  INTEGER :: n_Res_u,n_Res
  PARAMETER(n_Res_u=40,n_Res=100)
  INTEGER, SAVE :: units,n_dip
  TYPE(Subtract), DIMENSION(:), ALLOCATABLE, SAVE :: Res_u
  TYPE(rms_units), DIMENSION(:), ALLOCATABLE, SAVE :: atoms
  REAL(8), ALLOCATABLE, SAVE :: wc(:),xyz0(:,:)

!!$========================================================================
!!$-- Allocatables
!!$========================================================================

  REAL(8), DIMENSION(:), ALLOCATABLE, SAVE :: charge
  REAL(8), DIMENSION(:,:), ALLOCATABLE, SAVE :: dip
  REAL(8), DIMENSION(:,:), ALLOCATABLE, SAVE :: dip_corr
  REAL(8), DIMENSION(:,:), ALLOCATABLE, SAVE :: dists
CONTAINS

!!$======================== DECLARATIONS ================================*

  SUBROUTINE Initialize(atres,nbun,res,prsymb,beta,charges,sunit,ntap)
    
    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

    INTEGER :: atres(2,*),nbun,res(*),ntap
    CHARACTER(8) :: prsymb(*)
    CHARACTER(7) :: beta(*)
    REAL(8) :: charges(*),sunit

!!$----------------------- VARIABLES IN COMMON --------------------------*

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: i,k1,k2,n,nn
    INTEGER :: count,count2,count_res,i_f,nlist
    CHARACTER(80) :: errmsg
    INTEGER, DIMENSION (:), ALLOCATABLE :: index,trace
    LOGICAL, DIMENSION(:), ALLOCATABLE :: mask
    LOGICAL :: ok


!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    natom=ntap
    ALLOCATE(wc(natom),xyz0(3,natom))
    wc=0.0D0

    count_res=0
    DO n=1,nbun
       i_f=atres(1,n)
       DO k1=1,units
          IF(prsymb(res(i_f))(1:3) == Res_u(k1)%target(1:3)) THEN
             ok=.FALSE.
             DO k2=1,Res_u(k1)%n0
                ok=.FALSE.
                DO i=atres(1,n),atres(2,n)
                   IF(Res_u(k1)%beta(k2) == beta(i)) THEN
                      ok=.TRUE.
                   END IF
                END DO
                IF(.NOT. ok) EXIT
             END DO
             IF(ok) THEN
                count_res=count_res+1
             ELSE
                WRITE(errmsg,'(''In unit '',a8,'' wrong atom in the li&
                     &st. Abort. '')') Res_u(k1)%target(1:3) 
                CALL abort_now(errmsg)
             END IF
          END IF
       END DO
    END DO
    nbunits=count_res

    ALLOCATE(atoms(nbunits),index(units),mask(ntap))
    ALLOCATE(charge(ntap))
    SELECT CASE(rms_type)
    CASE('cm')
       charge(1:ntap)=1.0D0
    CASE('dipole')
       charge(1:ntap)=charges(1:ntap)*sunit
    CASE('mechanic')
       DO i=1,ntap
          charge(i)=SIGN(1.0D0,charges(i))
       END DO
    CASE DEFAULT
       STOP
    END SELECT
    
    nlist=0
    DO n=1,nbun
       i_f=atres(1,n)
       count=0
       index=0
       DO k1=1,units
          IF(prsymb(res(i_f))(1:3) == Res_u(k1)%target(1:3)) THEN
             ok=.FALSE.
             DO k2=1,Res_u(k1)%n0
                ok=.FALSE.
                DO i=atres(1,n),atres(2,n)
                   IF(Res_u(k1)%beta(k2) == beta(i)) THEN
                      ok=.TRUE.
                   END IF
                END DO
                IF(.NOT. ok) EXIT
             END DO
             IF(ok) THEN
                count=count+1
                index(count)=k1
             END IF
          END IF
       END DO
       DO nn=1,count
          nlist=nlist+1
          k1=index(nn) 
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
          atoms(nlist)%n0=count2
          atoms(nlist)%target=prsymb(res(atres(1,n)))
          ALLOCATE(atoms(nlist)%ind(count2))
          count2=0
          DO i=atres(1,n),atres(2,n)
             IF(mask(i)) THEN
                count2=count2+1
                atoms(nlist)%ind(count2)=i
             END IF
          END DO
       END DO
    END DO
    DEALLOCATE(index,mask)
    count=0
    DO n=1,nbunits
       IF(atoms(n)%n0 /= 0) THEN
          count=count+1
          DO i=1,atoms(n)%n0
             k1=atoms(n)%ind(i)
             WRITE(*,'(2i6,a5,2x,a5)') n,k1,atoms(n)%target,beta(k1)
          END DO
       END IF
    END DO
    IF(count /= nbunits) STOP
    n_dip=count
    ALLOCATE(dip(3,n_dip),dip_corr(n_dip,n_dip),dists(n_dip,n_dip))
    dip=0.0D0
    dip_corr=0.0D0
    dists=0.0D0

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

  END SUBROUTINE Initialize

  SUBROUTINE Start(xp0,yp0,zp0)
    IMPLICIT NONE 
    REAL(8) :: xp0(*),yp0(*),zp0(*)
    INTEGER :: i,n,nn,count
    REAL(8) :: sum_wc
    
    DO n=1,nbunits
       IF(atoms(n)%n0 /= 0) THEN
          DO nn=1,atoms(n)%n0
             i=atoms(n)%ind(nn)
             wc(i)=1.0D0
             xyz0(1,i)=xp0(i)
             xyz0(2,i)=yp0(i)
             xyz0(3,i)=zp0(i)
          END DO
       END IF
    END DO
    sum_wc=SUM(wc)
    wc=wc/sum_wc
  END SUBROUTINE Start

  SUBROUTINE Compute(co,oc,xp0,yp0,zp0)
    IMPLICIT NONE 

!!$----------------------------- ARGUMENTS ------------------------------*

    REAL(8) :: xp0(*),yp0(*),zp0(*),oc(3,3),co(3,3)

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: i,n,nn,m,c1,c2
    INTEGER :: count
    REAL(8) :: xc,yc,zc,sumx,sumy,sumz,fatoms,xd,yd,zd,xdn,ydn,zdn&
         &,xdm,ydm,zdm,sum_wc
    REAL(8), DIMENSION(:), ALLOCATABLE :: dipx,dipy,dipz
    
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    ALLOCATE(dipx(n_dip),dipy(n_dip),dipz(n_dip))
    counter=counter+1

    c1=0
    DO n=1,nbunits
       IF(atoms(n)%n0 /= 0) THEN
          c1=c1+1
          fatoms=1.0D0/DBLE(atoms(n)%n0)
          dipx(c1)=0.0D0
          dipy(c1)=0.0D0
          dipz(c1)=0.0D0
          DO nn=1,atoms(n)%n0
             i=atoms(n)%ind(nn)
             dipx(c1)=dipx(c1)+charge(i)*xp0(i)
             dipy(c1)=dipy(c1)+charge(i)*yp0(i)
             dipz(c1)=dipz(c1)+charge(i)*zp0(i)
          END DO
          dipx(c1)=dipx(c1)*fatoms
          dipy(c1)=dipy(c1)*fatoms
          dipz(c1)=dipz(c1)*fatoms
          dip(1,c1)=dip(1,c1)+dipx(c1)
          dip(2,c1)=dip(2,c1)+dipy(c1)
          dip(3,c1)=dip(3,c1)+dipz(c1)
       END IF
    END DO

    DO n=1,n_dip
       xdn=oc(1,1)*dipx(n)+oc(1,2)*dipy(n)+oc(1,3)*dipz(n)
       ydn=oc(2,1)*dipx(n)+oc(2,2)*dipy(n)+oc(2,3)*dipz(n)
       zdn=oc(3,1)*dipx(n)+oc(3,2)*dipy(n)+oc(3,3)*dipz(n)
       DO m=1,n_dip
          xdm=oc(1,1)*dipx(m)+oc(1,2)*dipy(m)+oc(1,3)*dipz(m)
          ydm=oc(2,1)*dipx(m)+oc(2,2)*dipy(m)+oc(2,3)*dipz(m)
          zdm=oc(3,1)*dipx(m)+oc(3,2)*dipy(m)+oc(3,3)*dipz(m)
          xd=xdm-xdn
          yd=ydm-ydn
          zd=zdm-zdn
          xd=xd-2.0D0*PBC(xd)
          yd=yd-2.0D0*PBC(yd)
          zd=zd-2.0D0*PBC(zd)
          xc=co(1,1)*xd+co(1,2)*yd+co(1,3)*zd
          yc=co(2,1)*xd+co(2,2)*yd+co(2,3)*zd
          zc=co(3,1)*xd+co(3,2)*yd+co(3,3)*zd
          dists(n,m)=dists(n,m)+xc*xc+yc*yc+zc*zc
          dip_corr(n,m)=dip_corr(n,m)+(dipx(n)*dipx(m)+dipy(n)&
               &*dipy(m)+dipz(n)*dipz(m)) 
       END DO
    END DO


    DEALLOCATE(dipx,dipy,dipz)
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

  END SUBROUTINE Compute
  SUBROUTINE Write_it(fstep)
    IMPLICIT NONE 

!!$----------------------------- ARGUMENTS ------------------------------*
    
    REAL(8) :: fstep

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: i,n,nn,m,k
    INTEGER :: count
    INTEGER(8) :: lda,info,lwork
    REAL(8) :: xc,yc,zc,xc2,yc2,zc2,fcounter,rms_tot,rms,snm,summo
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: mrirj
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: mrirj_or
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: Multiply
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: ddip

    REAL(8), DIMENSION(:,:), ALLOCATABLE :: U,VT,sigma
    REAL(8), DIMENSION(:), ALLOCATABLE :: sing,work
    
    
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    ALLOCATE(mrirj(n_dip,n_dip),mrirj_or(n_dip,n_dip),Multiply(n_dip&
         &,n_dip),ddip(3,n_dip))
    ALLOCATE(U(n_dip,n_dip),VT(n_dip,n_dip),sigma(n_Dip,n_Dip)&
         &,sing(n_dip)) 

    REWIND(krms_dip)

    fcounter=1.0D0/DBLE(counter)
    WRITE(krms_dip,'(''Ring Dipole Correlation at '',f12.3,'' fs'')')&
         & fstep 

    DO n=1,n_dip
       ddip(1,n)=dip(1,n)*fcounter
       ddip(2,n)=dip(2,n)*fcounter
       ddip(3,n)=dip(3,n)*fcounter
    END DO
    DO n=1,n_dip
       xc=ddip(1,n)
       yc=ddip(2,n)
       zc=ddip(3,n)
       DO m=1,n_dip
          snm=dip_corr(n,m)*fcounter
          mrirj(n,m)=snm-(xc*ddip(1,m)+yc*ddip(2,m)+zc*ddip(3,m))
          mrirj_or(n,m)=mrirj(n,m)
          WRITE(krms_dip,'(2i7,e15.7)') n,m,mrirj(n,m)
       END DO
    END DO
    IF(counter /= 1 .AND. RMS_Invert) THEN

       lda=SIZE(mrirj,1)
       lwork=5*lda
       ALLOCATE(work(lwork))
       CALL DGESVD('A','A',lda,lda,mrirj,lda,sing,U,lda,vt,lda,work,lwork,info)
       DEALLOCATE(work)
       sigma=0.0D0
       DO i=1,n_dip
          IF(DABS(sing(i)) > 1.0D-3) THEN
             sigma(i,i)=1.0D0/sing(i)
          END IF
       END DO

       vt=TRANSPOSE(vt)
       U=TRANSPOSE(U)
       U=MATMUL(sigma,U)
       mrirj=MATMUL(vt,U)

       Multiply=MATMUL(mrirj,mrirj_or)
       WRITE(*,*) Multiply(10,10),Multiply(2,4),Multiply(4,2)
    END IF
    WRITE(krms_dip,'(''Dipolar fluctuations inverse'')')
    DO n=1,n_dip
       DO m=n,n_dip
          snm=dists(n,m)*fcounter
          WRITE(krms_dip,'(2i7,e15.7,2x,f15.7)') n,m,mrirj(n,m),DSQRT(snm)
       END DO
    END DO
    DEALLOCATE(mrirj,mrirj_or,ddip)
    DEALLOCATE(U,VT,sigma,sing) 


!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

  END SUBROUTINE Write_it
  FUNCTION PBC(x)
    IMPLICIT NONE 
    REAL(8) :: x,PBC

    PBC=DNINT(0.5D0*x)
  END FUNCTION PBC
  SUBROUTINE Rotate(xp0,yp0,zp0)
    IMPLICIT NONE
    REAL(8) :: xp0(*),yp0(*),zp0(*)

    REAL(8), DIMENSION (:,:), ALLOCATABLE :: xyz,xyzfit
    REAL(8), DIMENSION (:), ALLOCATABLE :: work
    
    REAL(8) :: qt(4),errca,dr1,dr2
    INTEGER :: i,ii

    ALLOCATE(xyz(3,natom),xyzfit(3,natom),work(natom+1))
    
    
    DO i=1,natom 
       xyz(1,i)=xp0(i)
       xyz(2,i)=yp0(i)
       xyz(3,i)=zp0(i)
    END DO
    CALL xfit(xyz0,xyz,xyzfit,qt,wc,work,natom,errca)
    DO i=1,natom
       xp0(i)=xyzfit(1,i)
       yp0(i)=xyzfit(2,i)
       zp0(i)=xyzfit(3,i)
    END DO
    DEALLOCATE(xyz,xyzfit,work)
  END SUBROUTINE Rotate
  INCLUDE 'RMS_Subtract_Read.f90'
END MODULE RMS_Subtract_Mod
