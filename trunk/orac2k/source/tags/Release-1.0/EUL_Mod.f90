MODULE EUL_Mod
  USE xerror_mod
!!$***********************************************************************
!!$   Time-stamp: <2005-06-08 13:46:59 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Tue Jun  7 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************
  LOGICAL, SAVE :: eul_angles=.FALSE.
  TYPE lista
     INTEGER :: no
     INTEGER, DIMENSION (:), POINTER :: nb
  END TYPE lista

  INTEGER :: Nindex_Max
  PARAMETER(Nindex_Max=100)
  TYPE(LISTA), SAVE :: index(Nindex_Max)
  INTEGER, DIMENSION (:), ALLOCATABLE, SAVE :: index_i
  INTEGER, DIMENSION (:), ALLOCATABLE, SAVE :: index_l
  REAL(8), ALLOCATABLE, SAVE :: wc(:),xyz0(:,:)
  REAL(8), ALLOCATABLE :: xyz(:,:),xyzfit(:,:),work(:)
  INTEGER, SAVE :: natom,N_Eul,Eul_domain,kdump

!!$---- This MODULE is part of the program ORAC ----*

CONTAINS
  SUBROUTINE Initialize(wca,xpt,ypt,zpt)
    IMPLICIT NONE 
    REAL(8) :: wca(:),xpt(:),ypt(:),zpt(:)
    INTEGER :: i,count
    REAL(8) :: total

    natom=SIZE(wca)
    count=0
    DO i=1,natom
       IF(wca(i) /= 0.0D0) count=count+1
    END DO
    N_eul=count
    ALLOCATE(index_i(N_eul))
    ALLOCATE(wc(N_eul),xyz0(3,N_eul))
    ALLOCATE(xyz(3,N_eul),xyzfit(3,N_eul),work(N_eul+1))

    count=0
    DO i=1,natom
       IF(wca(i) /= 0.0D0) THEN
          count=count+1
          index_i(count)=i
       END IF
    END DO
    Total=SUM(wca)
    DO i=1,N_eul
       wc(i)=wca(index_i(i))
       xyz0(1,i)=xpt(index_i(i))
       xyz0(2,i)=ypt(index_i(i))
       xyz0(3,i)=zpt(index_i(i))
    END DO
    Total=SUM(wc)
    wc=wc/Total
  END SUBROUTINE Initialize

  SUBROUTINE Create(m)
    IMPLICIT NONE
    INTEGER :: m
    ALLOCATE(index_l(m+1))
    index(:)%no=0
  END SUBROUTINE Create
  SUBROUTINE Add
    IMPLICIT NONE
    INTEGER, SAVE :: count=0
    INTEGER :: m
    count=count+1
    IF(count > Nindex_Max) CALL &
         & abort_now('Domain No. Larger than expected. Abort now!')
    m=index_l(1)
    index(count)%no=m
    ALLOCATE(index(count)%nb(m))
    index(count)%nb(:)=index_l(2:m+1)
    EUL_Domain=count
  END SUBROUTINE Add
  SUBROUTINE Compute(xp0,yp0,zp0)
    IMPLICIT NONE 
    REAL(8) :: xp0(*),yp0(*),zp0(*)

    REAL(8) :: Angles(3)
    REAL(8), DIMENSION(:), ALLOCATABLE :: wcc_l
    INTEGER, SAVE :: N_Calls=0
    INTEGER :: i,jj,j,k
    REAL(8) :: Total
    REAL(8) :: qt(4),errca,dr1,dr2,v(3),Theta,dct
    REAL(8), SAVE :: pi
    REAL(8), ALLOCATABLE :: xyz00(:,:)
    
    N_Calls=N_Calls+1
    ALLOCATE (wcc_l(N_eul))
    DO i=1,N_Eul
       xyz(1,i)=xp0(index_i(i))
       xyz(2,i)=yp0(index_i(i))
       xyz(3,i)=zp0(index_i(i))
    END DO
       
    IF(N_Calls == 1) THEN
       pi=4.0*ATAN(1.0D0)
       DO i=1,eul_domain
          wcc_l=0.0D0
          DO k=1,index(i)%no
             jj=index(i)%nb(k)
             wcc_l(jj)=wc(jj)
          END DO
          Total=SUM(wcc_l)
          wcc_l=wcc_l/Total

          CALL xfitta(xyz,xyz0,xyzfit,dct,qt,wcc_l,work,N_eul,errca)

          DO j=1,N_eul 
             IF(wcc_l(j) /= 0.0D0) THEN
                xyz0(1,j)=xyzfit(1,j)
                xyz0(2,j)=xyzfit(2,j)
                xyz0(3,j)=xyzfit(3,j)
             END IF
          END DO
       END DO
    END IF
    
    DO i=1,eul_domain
       wcc_l=0.0D0
       DO k=1,index(i)%no
          jj=index(i)%nb(k)
          wcc_l(jj)=wc(jj)
       END DO
       Total=SUM(wcc_l)
       wcc_l=wcc_l/Total
       
       CALL xfitta(xyz,xyz0,xyzfit,dct,qt,wcc_l,work,N_eul,errca)
       Theta=2.0D0*ACOS(qt(1))
       v=0.0D0
       IF(ABS(Theta)  > 1.0D-10) THEN
          v(1)=2.0D0*qt(1)*qt(2)/SIN(Theta)
          v(2)=2.0D0*qt(1)*qt(3)/SIN(Theta)
          v(3)=2.0D0*qt(1)*qt(4)/SIN(Theta)
       END IF
       Theta=Theta*180.0D0/pi
       WRITE(kdump,'(i8,i5,3x,6(f14.8,1x))') N_Calls,i,dct,errca,Theta,v(1:3)
    END DO
    DEALLOCATE(wcc_l)
  END SUBROUTINE Compute
  SUBROUTINE xfitta(xyz,xyz0,xyzfit,dct,qt,w1,w2,ntap,error)
    IMPLICIT none
    REAL(8) ::  w1(*),w2(*),dct
    REAL(8) ::  xyz(3,*),xyz0(3,*),xyzfit(3,*),qt(4),error
    INTEGER :: ntap

    INTEGER :: i,iret
    REAL(8) ::  dcm(3),q(0:7)
    CHARACTER(80) ::  msg
      
    DO i=1,ntap
       w2(i)=w1(i)
    END DO
    
    CALL rigfit(0,ntap,xyz,xyz0,w1,w2,q,dcm,xyzfit,error,iret,msg)
!!$    dct=DSQRT(SUM(dcm(1:3)**2))
    dct=dcm(3)
    DO i=1,4
       qt(i)=q(i-1)
    END DO
  END SUBROUTINE xfitta
END MODULE EUL_Mod
