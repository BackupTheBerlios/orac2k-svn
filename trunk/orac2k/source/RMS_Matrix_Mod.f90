MODULE RMS_Matrix_Mod
!!$***********************************************************************
!!$   Time-stamp: <2005-08-11 09:10:46 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Jun  3 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

  REAL(8), POINTER, SAVE :: rms_xi(:,:),rms_xixj(:,:),rms_sxi(:,:)&
       &,rms_sxixj(:,:)
  REAL(8), ALLOCATABLE, SAVE :: wc(:),xyz0(:,:)
  INTEGER, SAVE :: N_rms, krms_matrix,krms_matrix_plot,natom
  LOGICAL, SAVE :: rms_matrix=.FALSE.,rms_matrix_avg=.FALSE.&
       &,rms_matrix_plot=.FAlSE.
  INTEGER, DIMENSION (:), POINTER, SAVE :: index
  INTEGER, SAVE :: N_calls=0,N_Calls_avg=0
  INTEGER, SAVE :: Write_Freq=100
  INTEGER, SAVE :: Navg=4166

!!$======================== DECLARATIONS ================================*

CONTAINS
  SUBROUTINE Initialize(wca,xpt,ypt,zpt)
    IMPLICIT NONE 
    REAL(8) :: wca(:),xpt(:),ypt(:),zpt(:)
    INTEGER :: i,count

    natom=SIZE(wca)
    ALLOCATE(wc(natom),xyz0(3,natom))
    wc=wca/SUM(wca)
    count=0
    DO i=1,natom
       IF(wca(i) /= 0.0D0) count=count+1
       xyz0(1,i)=xpt(i)
       xyz0(2,i)=ypt(i)
       xyz0(3,i)=zpt(i)
    END DO
    ALLOCATE(index(count+1))

    count=0
    DO i=1,natom
       IF(wca(i) /= 0.0D0) THEN
          count=count+1
          index(count)=i
       END IF
    END DO
    N_rms=INT(SUM(wca))    
    ALLOCATE(rms_xi(3,N_rms),rms_xixj(N_rms,N_rms))
    ALLOCATE(rms_sxi(3,N_rms),rms_sxixj(N_rms,N_rms))
    rms_xi=0.0D0
    rms_xixj=0.0D0
    rms_sxi=0.0D0
    rms_sxixj=0.0D0
  END SUBROUTINE Initialize
  SUBROUTINE Compute(xp0,yp0,zp0)
    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*
  
    REAL(8) :: xp0(*),yp0(*),zp0(*)

!!$----------------------- VARIABLES IN COMMON --------------------------*

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: i,j
    REAL(8) :: xpi,ypi,zpi,xpj,ypj,zpj

    N_Calls=N_calls+1
    DO i=1,N_rms
       xpi=xp0(index(i))
       ypi=yp0(index(i))
       zpi=zp0(index(i))
       rms_xi(1,i)=rms_xi(1,i)+xpi
       rms_xi(2,i)=rms_xi(2,i)+ypi
       rms_xi(3,i)=rms_xi(3,i)+zpi
       DO j=i,N_rms
          xpj=xp0(index(j))
          ypj=yp0(index(j))
          zpj=zp0(index(j))
          rms_xixj(i,j)=rms_xixj(i,j)+xpi*xpj+ypi*ypj+zpi*zpj
       END DO
    END DO
  END SUBROUTINE Compute
  SUBROUTINE Compute_avg(xp0,yp0,zp0)
    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*
  
    REAL(8) :: xp0(*),yp0(*),zp0(*)

!!$------------------------- LOCAL VARIABLES ----------------------------*

    REAL(8), DIMENSION (:,:,:), ALLOCATABLE, SAVE :: xi_s
    INTEGER :: i,j,k
    REAL(8) :: xpi,ypi,zpi,xpj,ypj,zpj,f_avg,xixj,dxidxj

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    IF(.NOT. ALLOCATED(xi_s)) THEN
       ALLOCATE(xi_s(3,N_rms,Navg))
    END IF
    f_avg=1.0D0/DFLOAT(Navg)
    IF(N_Calls < Navg) THEN
       N_Calls=N_calls+1
       DO i=1,N_rms
          xpi=xp0(index(i))
          ypi=yp0(index(i))
          zpi=zp0(index(i))
          xi_s(1,i,N_Calls)=xpi
          xi_s(2,i,N_Calls)=ypi
          xi_s(3,i,N_Calls)=zpi
       END DO
    ELSE
       N_Calls_avg=N_Calls_avg+1
       IF(N_Calls_avg == 1) THEN
          DO k=1,Navg
             DO i=1,N_rms
                xpi=xi_s(1,i,k)
                ypi=xi_s(2,i,k)
                zpi=xi_s(3,i,k)
                rms_sxi(1,i)=rms_sxi(1,i)+xpi
                rms_sxi(2,i)=rms_sxi(2,i)+ypi
                rms_sxi(3,i)=rms_sxi(3,i)+zpi
                DO j=i,N_rms
                   xpj=xi_s(1,j,k)
                   ypj=xi_s(2,j,k)
                   zpj=xi_s(3,j,k)
                   rms_sxixj(i,j)=rms_sxixj(i,j)+(xpi*xpj+ypi*ypj+zpi*zpj)
                END DO
             END DO
          END DO
       ELSE
          DO i=1,N_rms
             xpi=xi_s(1,i,1)
             ypi=xi_s(2,i,1)
             zpi=xi_s(3,i,1)
             rms_sxi(1,i)=rms_sxi(1,i)-xpi
             rms_sxi(2,i)=rms_sxi(2,i)-ypi
             rms_sxi(3,i)=rms_sxi(3,i)-zpi
             DO j=i,N_rms
                xpj=xi_s(1,j,1)
                ypj=xi_s(2,j,1)
                zpj=xi_s(3,j,1)
                rms_sxixj(i,j)=rms_sxixj(i,j)-(xpi*xpj+ypi*ypj+zpi*zpj)
             END DO
          END DO
          DO i=1,N_rms
             xpi=xp0(index(i))
             ypi=yp0(index(i))
             zpi=zp0(index(i))
             rms_sxi(1,i)=rms_sxi(1,i)+xpi
             rms_sxi(2,i)=rms_sxi(2,i)+ypi
             rms_sxi(3,i)=rms_sxi(3,i)+zpi
             DO j=i,N_rms
                xpj=xp0(index(j))
                ypj=yp0(index(j))
                zpj=zp0(index(j))
                rms_sxixj(i,j)=rms_sxixj(i,j)+(xpi*xpj+ypi*ypj+zpi*zpj)
             END DO
          END DO
          
       END IF
       
       DO k=1,Navg-1
          xi_s(:,:,k)=xi_s(:,:,k+1)
       END DO
       DO i=1,N_rms
          xpi=xp0(index(i))
          ypi=yp0(index(i))
          zpi=zp0(index(i))
          xi_s(1,i,Navg)=xpi
          xi_s(2,i,Navg)=ypi
          xi_s(3,i,Navg)=zpi
       END DO

       DO i=1,N_rms
          xpi=rms_sxi(1,i)*f_avg
          ypi=rms_sxi(2,i)*f_avg
          zpi=rms_sxi(3,i)*f_avg
          DO j=i,N_rms
             xpj=rms_sxi(1,j)*f_avg
             ypj=rms_sxi(2,j)*f_avg
             zpj=rms_sxi(3,j)*f_avg
             xixj=rms_sxixj(i,j)*f_avg
             dxidxj=xixj-(xpi*xpj+ypi*ypj+zpi*zpj)
             rms_xixj(i,j)=rms_xixj(i,j)+dxidxj
          END DO
       END DO
    END IF
  END SUBROUTINE Compute_avg
  SUBROUTINE Write_it
    IMPLICIT NONE 
    INTEGER :: i,j
    REAL(8) :: xpi,ypi,zpi,xpj,ypj,zpj,xixj,dxidxj,f_Calls,xxi,xxj
    REAL(8), DIMENSION (:), ALLOCATABLE, SAVE :: xxii

    IF(MOD(N_Calls,Write_Freq) == 0) THEN
       REWIND(krms_matrix)
       ALLOCATE(xxii(N_rms))
       f_Calls=1.0D0/DFLOAT(N_Calls)
       WRITE(krms_matrix,'(''Averages After '',i7,'' Calls '')') N_Calls
       DO i=1,N_rms
          xpi=rms_xi(1,i)*f_Calls
          ypi=rms_xi(2,i)*f_Calls
          zpi=rms_xi(3,i)*f_Calls
          xixj=rms_xixj(i,i)*f_Calls
          xxii(i)=xixj-(xpi*xpi+ypi*ypi+zpi*zpi)
          WRITE(krms_matrix,'(2i5,e20.12)') i,i,xxii(i)
          xxii(i)=DSQRT(xxii(i))
       END DO

       DO i=1,N_rms
          xpi=rms_xi(1,i)*f_Calls
          ypi=rms_xi(2,i)*f_Calls
          zpi=rms_xi(3,i)*f_Calls
          xxi=xxii(i)
          DO j=i,N_rms
             xxj=xxii(j)*xxi
             xpj=rms_xi(1,j)*f_Calls
             ypj=rms_xi(2,j)*f_Calls
             zpj=rms_xi(3,j)*f_Calls
             xixj=rms_xixj(i,j)*f_Calls
             dxidxj=xixj-(xpi*xpj+ypi*ypj+zpi*zpj)
             WRITE(krms_matrix,'(2i5,e20.12)') i,j,dxidxj/xxj
          END DO
       END DO
       DEALLOCATE(xxii)
    END IF
  END SUBROUTINE Write_it
  SUBROUTINE Write_it_avg
    IMPLICIT NONE 
    INTEGER :: i,j
    REAL(8) :: xpi,ypi,zpi,xpj,ypj,zpj,xixj,dxidxj,f_Calls,xxi,xxj
    REAL(8), DIMENSION (:), ALLOCATABLE, SAVE :: xxii

    IF(N_Calls_avg /=0 .AND. MOD(N_Calls_avg,Write_Freq) == 0) THEN
       REWIND(krms_matrix)
       f_Calls=1.0D0/DFLOAT(N_Calls_avg)
       WRITE(krms_matrix,'(''Averages After '',i7,'' Calls '')') N_Calls_avg
       ALLOCATE(xxii(N_rms))
       DO i=1,N_rms
          xxii(i)=rms_xixj(i,i)*f_Calls
          WRITE(krms_matrix,'(2i5,e20.12)') i,i,xxii(i)
          xxii(i)=DSQRT(xxii(i))
       END DO

       DO i=1,N_rms
          xxi=xxii(i)
          DO j=i,N_rms
             xxj=xxii(j)*xxi
             dxidxj=rms_xixj(i,j)*f_Calls
             WRITE(krms_matrix,'(2i5,e20.12)') i,j,dxidxj/xxj
          END DO
       END DO
    END IF
  END SUBROUTINE Write_it_avg
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
  SUBROUTINE Plot
    IMPLICIT NONE 
    INTEGER :: i,k
    REAL(8) :: xpi,ypi,zpi,dr,f_Calls

    IF(N_Calls /=0 .AND. MOD(N_Calls,Write_Freq) == 0) THEN
       f_Calls=1.0D0/DFLOAT(N_Calls)
       REWIND(krms_matrix_plot)
       dr=0.0D0
       k=0
       DO i=1,N_rms
          xpi=rms_xi(1,i)*f_Calls
          ypi=rms_xi(2,i)*f_Calls
          zpi=rms_xi(3,i)*f_Calls
          WRITE(krms_matrix_plot,1)'ATOM  ',i,'CA   ','ALA',i,xpi,ypi&
               &,zpi,dr,DFLOAT(k)
       END DO
       WRITE(krms_matrix_plot,'(''TER'')')
       DO i=1,N_rms
          xpi=xyz0(1,index(i))
          ypi=xyz0(2,index(i))
          zpi=xyz0(3,index(i))
          WRITE(krms_matrix_plot,1)'ATOM  ',i,' CA  ','ALA',i,xpi,ypi&
               &,zpi,dr,DFLOAT(k)
       END DO
    END IF
1   FORMAT(a6,i5,1x,a5,a3,2x,i4,4x,3f8.3,f8.4,f4.1)
  END SUBROUTINE Plot
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*
END MODULE RMS_Matrix_Mod
