MODULE VORONOI_Mod

!!$***********************************************************************
!!$   Time-stamp: <2006-04-06 12:05:33 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Wed Apr  5 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This MODULE is part of the program ORAC ----*

#include "config.h"

  USE Xerror_Mod

  USE INPUT_Mod, ONLY: Read_String, Parser, err_open,err_end,err_unr&
       &,err_fnf,err_args

  INTEGER :: max_neigh,pig_nnl,pnnlpp_vor,maxpla,maxver
  PARAMETER ( &
       & max_neigh=_MAX_NEIGH_, &
       & pnnlpp_vor=MAXPLA+1, &
       & maxpla=MAXPLA, &
       & maxver=MAXVER)

!!$======================== Input Parameters=============================*

  LOGICAL, SAVE :: voronoi=.FALSE.,heavy=.FALSE.,access=.FALSE.&
       &,volume=.FALSE.,fluct=.FALSE.,neighbor=.FALSE.,rewind_vor=.FALSE.
  INTEGER, SAVE :: nvoronoi=1,kvoronoi=0,nfluct=2
  REAL(8), SAVE :: cutoff
  CHARACTER(80), SAVE :: filename


!!$======================== Input Parameters end =========================*

!!$============================= Arrays ==================================*

  INTEGER, SAVE :: ntap,nbun,nunits
  INTEGER, DIMENSION (:,:), ALLOCATABLE, SAVE :: nnlpp_vor
  INTEGER, DIMENSION (:,:), ALLOCATABLE, SAVE :: res
  REAL(8), DIMENSION (:,:), ALLOCATABLE, SAVE :: area_vor
  REAL(8), DIMENSION (:), ALLOCATABLE, SAVE :: volume_vor
  REAL(8), SAVE :: vol_slv
  REAL(8), ALLOCATABLE, SAVE :: vol_res(:),vol_slt(:),area_slt(:)&
       &,area_slv(:),area_frac(:),area_tot(:),vol_type(:)
  INTEGER, ALLOCATABLE, SAVE :: index(:),vol_typep(:)
  
!!$============================= Arrays End ==============================*

!!$============================= Fluctations =============================*

  REAL(8), DIMENSION(:,:), ALLOCATABLE, SAVE :: dvij
  REAL(8), DIMENSION(:,:), ALLOCATABLE, SAVE :: dist

  REAL(8), DIMENSION(:), ALLOCATABLE, SAVE :: dvi
  INTEGER, SAVE :: fluct_calls=0
  
!!$============================= Fluctations End =========================*


CONTAINS
  SUBROUTINE Init(mma,natom,nbuna,mres,grppt,ss_index)

!!$======================== DECLARATIONS ================================*

    IMPLICIT none
    
    INTEGER :: mma,natom,nbuna,mres(2,*), grppt(2,*),ss_index(*)

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: n,m,i,count,nslv
    LOGICAL :: ok

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*
    

    nbun=nbuna
    ntap=natom

    ALLOCATE(nnlpp_vor(pnnlpp_vor,mma))
    ALLOCATE(area_vor(pnnlpp_vor,mma),volume_vor(ntap))    
    count=0
    DO n=1,nbun
       ok=.FALSE.
       DO m=mres(1,n),mres(2,n)
          DO i=grppt(1,m),grppt(2,m)
             IF(ss_index(i) == 1) ok=.TRUE.               
          END DO
       END DO
       IF(ok) THEN
          count=count+1
       END IF
    END DO
    ALLOCATE(res(2,count))
    ALLOCATE(vol_res(nbun),vol_slt(count),area_slt(nbun)&
         &,area_slv(nbun),area_frac(nbun),area_tot(nbun)&
         &,vol_type(nbun),index(nbun),vol_typep(nbun))
    ALLOCATE(dvij(count+1,count+1),dvi(count+1),dist(count+1,count+1))
    dvij=0.0D0
    dvi=0.0D0
    nunits=count
    

    count=0
    DO n=1,nbun
       ok=.FALSE.
       DO m=mres(1,n),mres(2,n)
          DO i=grppt(1,m),grppt(2,m)
             IF(ss_index(i) == 1) ok=.TRUE.               
          END DO
       END DO
       IF(ok) THEN
          count=count+1
          res(1,count)=grppt(1,mres(1,n))
          res(2,count)=grppt(2,mres(2,n))
       END IF
    END DO


    IF(fluct) THEN
       access=.TRUE.
       volume=.TRUE.
    END IF
  END SUBROUTINE Init
  SUBROUTINE Fluctuations(xp0,yp0,zp0)

!!$----------------------------- ARGUMENTS ------------------------------*

    IMPLICIT NONE 

    REAL(8) :: xp0(*),yp0(*),zp0(*)

!!$----------------------- VARIABLES IN COMMON --------------------------*

!!$------------------------- LOCAL VARIABLES ----------------------------*

    REAL(8), DIMENSION(:,:), ALLOCATABLE, SAVE :: Kij,Kij_orig
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: U,VT,sigma,Multiply
    REAL(8), DIMENSION(:), ALLOCATABLE, SAVE :: vols,work,sing
    
    INTEGER(8) :: lda,info,lwork
    INTEGER :: n_U,i,j
    REAL(8) :: f_Calls,xpi,xpj

    IF(rewind_vor) REWIND(kvoronoi)

    n_U=nunits+1
    
    ALLOCATE(vols(n_U))
         

    vols(1:nunits)=vol_res(1:nunits)
    vols(n_U)=vol_slv

    fluct_calls=fluct_calls+1
    f_Calls=1.0D0/DBLE(fluct_calls)

    DO i=1,n_U
       xpi=vols(i)
       dvi(i)=dvi(i)+xpi
       DO j=i,n_U
          xpj=vols(j)
          dvij(i,j)=dvij(i,j)+xpi*xpj
       END DO
    END DO

    CALL Get_Dist

    DEALLOCATE(vols)

    IF(MOD(fluct_calls,nfluct) == 0) THEN
       ALLOCATE(U(n_U,n_U),VT(n_U,n_U),sigma(n_U,n_U),sing(n_U)&
            &,kij_orig(n_U,n_U),kij(n_U,n_U),Multiply(n_U,n_U))
       Kij=0.0D0
       DO i=1,n_U
          xpi=dvi(i)*f_Calls
          DO j=i,n_U
             xpj=dvi(j)*f_Calls
             Kij(i,j)=dvij(i,j)*f_Calls-xpi*xpj
             Kij_orig(i,j)=Kij(i,j)

             Kij(j,i)=Kij(i,j)
             Kij_orig(j,i)=Kij_orig(i,j)
          END DO
       END DO


       WRITE(kvoronoi,'(''Average Volumes '',i7)') fluct_calls
       DO i=1,n_U
          WRITE(kvoronoi,'(i7,e15.7)') i,dvi(i)*f_Calls
       END DO

       WRITE(kvoronoi,'(''Direct Fluctuations at step '',i7)') fluct_calls

       DO i=1,n_U
          DO j=1,n_U
             WRITE(kvoronoi,'(2i7,e15.7)') i,j,Kij(i,j)             
          END DO
       END DO


       lda=SIZE(Kij,1)
       lwork=5*lda

       ALLOCATE(work(lwork))
       CALL DGESVD('A','A',lda,lda,Kij,lda,sing,U,lda,vt,lda,work,lwork,info)
       DEALLOCATE(work)
       sigma=0.0D0
       DO i=1,n_U
          IF(DABS(sing(i)) > 1.0D-6) THEN
             sigma(i,i)=1.0D0/sing(i)
          END IF
       END DO

       vt=TRANSPOSE(vt)
       U=TRANSPOSE(U)
       U=MATMUL(sigma,U)
       Kij=MATMUL(vt,U)

       Multiply=MATMUL(Kij,Kij_orig)
       WRITE(*,*) Multiply(10,10),Multiply(2,4),Multiply(4,2)

       WRITE(kvoronoi,'(''Inverse Fluctuations at step '',i7)') fluct_calls

       DO i=1,n_U
          DO j=i,n_U
             WRITE(kvoronoi,'(2i7,e15.7,2x,f15.7)') i,j,Kij(i,j)&
                  &,dist(i,j)*f_Calls
          END DO
       END DO

       DEALLOCATE(U,VT,sigma,sing,kij_orig,Kij,Multiply)
    END IF

  CONTAINS
    SUBROUTINE Get_Dist
      IMPLICIT NONE 
      INTEGER :: n,m
      REAL(8), DIMENSION(:), ALLOCATABLE :: rx,ry,rz
      REAL(8) :: fcount,xpi,ypi,zpi,xc,yc,zc

      ALLOCATE(rx(nunits),ry(nunits),rz(nunits))
      
      DO n=1,nunits
         rx(n)=0.0D0
         ry(n)=0.0D0
         rz(n)=0.0D0
         DO m=res(1,n),res(2,n)
            rx(n)=rx(n)+xp0(m)
            ry(n)=ry(n)+yp0(m)
            rz(n)=rz(n)+zp0(m)
         END DO
         fcount=1.0D0/DBLE(res(2,n)-res(1,n)+1)
         rx(n)=rx(n)*fcount
         ry(n)=ry(n)*fcount
         rz(n)=rz(n)*fcount
      END DO
      DO n=1,nunits
         xpi=rx(n)
         ypi=ry(n)
         zpi=rz(n)
         DO m=n,nunits
            xc=xpi-rx(m)
            yc=ypi-ry(m)
            zc=zpi-rz(m)
            dist(n,m)=dist(n,m)+SQRT(xc**2+yc**2+zc**2)
         END DO
      END DO
      DEALLOCATE(rx,ry,rz)

    END SUBROUTINE Get_Dist
  END SUBROUTINE Fluctuations

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*
  INCLUDE 'VORONOI_Read.f90'
END MODULE VORONOI_Mod
