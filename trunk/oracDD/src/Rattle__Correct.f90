!!$/---------------------------------------------------------------------\
!!$   Copyright  © 2006-2007 Massimo Marchi <Massimo.Marchi at cea.fr>   |
!!$                                                                      |
!!$    This software is a computer program named oracDD whose            |
!!$    purpose is to simulate and model complex molecular systems.       |
!!$    The code is written in fortran 95 compliant with Technical        |
!!$    Report TR 15581, and uses MPI-1 routines for parallel             |
!!$    coding.                                                           |
!!$                                                                      |
!!$    This software is governed by the CeCILL license under             |
!!$    French law and abiding by the rules of distribution of            |
!!$    free software.  You can  use, modify and/ or redistribute         |
!!$    the software under the terms of the CeCILL icense as              |
!!$    circulated by CEA, CNRS and INRIA at the following URL            |
!!$    "http://www.cecill.info".                                         |
!!$                                                                      |
!!$    As a counterpart to the access to the source code and rights      |
!!$    to copy, modify and redistribute granted by the license,          |
!!$    users are provided only with a limited warranty and the           |
!!$    software's author, the holder of the economic rights, and         |
!!$    the successive licensors have only limited liability.             |
!!$                                                                      |
!!$    The fact that you are presently reading this means that you       |
!!$    have had knowledge of the CeCILL license and that you accept      |
!!$    its terms.                                                        |
!!$                                                                      |
!!$    You should have received a copy of the CeCILL license along       |
!!$    with this program; if not, you can collect copies on the URL's    |
!!$    "http://www.cecill.info/licences/Licence_CeCILL_V2-en.html"       |
!!$    "http://www.cecill.info/licences/Licence_CeCILL_V2-fr.html"       |
!!$                                                                      |
!!$----------------------------------------------------------------------/
FUNCTION Correct_(dt,xp0a,yp0a,zp0a,vpxa,vpya,vpza) RESULT(out)
  LOGICAL :: out
  REAL(8) :: dt,xp0a(:),yp0a(:),zp0a(:),vpxa(:),vpya(:),vpza(:)
  INTEGER, POINTER :: cnstp(:,:)
  REAL(8), POINTER :: dssp(:),coeffp(:),xp1(:),yp1(:),zp1(:),dssip(:)&
       &,coeffip(:) 
  TYPE(Rattle__Type2), ALLOCATABLE, TARGET :: dssi(:),coeffi(:)
  TYPE(Rattle__Type2), ALLOCATABLE, TARGET :: xxp1(:),yyp1(:),zzp1(:)
  INTEGER :: nn,m,n,n0
  INTEGER :: la,lb,k,iter,iox,i,iter1,ka,lab,lbb&
       &,kb,count1
  REAL(8) ::  xab,yab,zab,dpx,dpy,dpz,dpp,dps,gg,amsla,amslb&
       &,dpax,dpay,dpaz,dpbx,dpby,dpbz
  REAL(8) ::  gab,gg1,gg2,dcnt,xk,yk,zk,tol_mim,aux1,aux2,aux3&
       &,det,a11,a12,a13,a21,a22,a23,a31,a32,a33,b11,b12&
       &,b13,b21,b22,b23,b31,b32,b33,gcpu_mm,vfcp_mm,tfcp_mm&
       &,tdelta_mm,elapse,vdpx,vdpy,vdpz,xd1,yd1,zd1
  
  REAL(8), PARAMETER :: tol=1.0D-6
  REAL(8), SAVE ::tol_min=1.0D-5,zero=0.0_8,one=1.0_8&
       &,two=2.0_8,three=3.0_8,four=4.0_8,tolh=tol
  REAL(8) ::   aux,xd,yd,zd
  INTEGER :: info
  LOGICAL, ALLOCATABLE :: mask(:)

  out=.TRUE.
  IF(.NOT. Rattle__Param % switch) RETURN
  CALL Gather_Atoms
  
  ALLOCATE(xxp1(nc),yyp1(nc),zzp1(nc),dssi(nc),coeffi(nc),mask(natom))
  
  DO nn=1,nc
     cnstp=>cnst(nn) % n1
     n0=SIZE(cnstp,2)
     ALLOCATE(xxp1(nn) % g(n0))
     ALLOCATE(yyp1(nn) % g(n0))
     ALLOCATE(zzp1(nn) % g(n0))
     ALLOCATE(dssi(nn) % g(n0))
     ALLOCATE(coeffi(nn) % g(n0))
  END DO
  DO nn=1,nc
     cnstp=>cnst(nn) % n1
     dssp=>dss(nn) % g
     coeffp=>coeff(nn) % g
     xp1=>xxp1(nn) % g
     yp1=>yyp1(nn) % g
     zp1=>zzp1(nn) % g
     n0=SIZE(cnstp,2)
     coeffi(nn) % g=1.0D0/coeffp
     dssi(nn) % g=1.0D0/dssp
     DO m=1,n0
        la=cnstp(1,m)
        lb=cnstp(2,m)
        xp1(m)=xp0(la)-xp0(lb)
        yp1(m)=yp0(la)-yp0(lb)
        zp1(m)=zp0(la)-zp0(lb)
     END DO
  END DO

  iter1=0
  DO nn=1,nc
     cnstp=>cnst(nn) % n1
     dssip=>dssi(nn) % g
     coeffip=>coeffi(nn) % g
     xp1=>xxp1(nn) % g
     yp1=>yyp1(nn) % g
     zp1=>zzp1(nn) % g
     n0=SIZE(cnstp,2)
     IF(n0 > Rattle__Param % mim_Max .OR. (.NOT. Rattle__Param % mim)) THEN
        iter=0
1000    CONTINUE
        iox=0
        DO k=1,n0
           la=cnstp(1,k)
           lb=cnstp(2,k)
           dpx=xp1(k)
           dpy=yp1(k)
           dpz=zp1(k)
           vdpx=vpx(la)-vpx(lb)
           vdpy=vpy(la)-vpy(lb)
           vdpz=vpz(la)-vpz(lb)
           dpp=dpx*vdpx+dpy*vdpy+dpz*vdpz
           IF(DABS(dpp) > tolh) THEN
              gab=-dpp*dssip(k)*coeffip(k)
              iox=1
              gg1=mass0(la)*gab
              gg2=-mass0(lb)*gab
              vpx(la)=vpx(la)+gg1*dpx
              vpx(lb)=vpx(lb)+gg2*dpx
              
              vpy(la)=vpy(la)+gg1*dpy
              vpy(lb)=vpy(lb)+gg2*dpy
              
              vpz(la)=vpz(la)+gg1*dpz
              vpz(lb)=vpz(lb)+gg2*dpz
           END IF
        END DO
        IF(iox.ne.0) THEN
           iter1=iter1+1
           iter=iter+1
           IF(iter.GT.5000)THEN
              errmsg_f='While RATTLEing: The iteration procedure did n&
                   &ot converge.'
              CALL Add_Errors(-1,errmsg_f)
              out=.FALSE.
              RETURN
           END IF
           GOTO 1000
        END IF
     ELSE
        
!!$=======================================================================
!!$---  Use Matrix Inversion Method for Constraints ---------------------- 
!!$=======================================================================
!!$
!!$=========== Build constraint matrix ====================================

        DO ka=1,n0
           la=cnstp(1,ka)
           lb=cnstp(2,ka)
           xc(ka)=vpx(la)-vpx(lb)
           yc(ka)=vpy(la)-vpy(lb)
           zc(ka)=vpz(la)-vpz(lb)
        END DO
        
!!$======== First loop on constraints k1 =================================
        
        
        count1=0
        DO ka=1,n0
           la=cnstp(1,ka)
           lb=cnstp(2,ka)
           dd(ka)=xc(ka)*xp1(ka)+yc(ka)*yp1(ka)+zc(ka)*zp1(ka)
           xd1=xp1(ka)
           yd1=yp1(ka)
           zd1=zp1(ka)
           
!!$============= Second loop on constraints k2 ===========================
                     
           DO kb=1,n0
              lab=cnstp(1,kb)
              lbb=cnstp(2,kb)
                  
              xd=xp1(kb)
              yd=yp1(kb)
              zd=zp1(kb)
                  
              aux1=mass0(la)
              aux2=mass0(lb)
                  
              IF(la /= lab .AND. la /= lbb .AND. lb /= lab.AND. lb /= lbb) THEN
                 mat(ka,kb)=zero
              ELSE
                 dpax=zero
                 dpay=zero
                 dpaz=zero
                 dpbx=zero
                 dpby=zero
                 dpbz=zero
                 IF(la .EQ. lab) THEN
                    dpax=-xd
                    dpay=-yd
                    dpaz=-zd
                 ELSE IF(la .EQ. lbb) THEN
                    dpax=xd
                    dpay=yd
                    dpaz=zd
                 END IF
                 
                 IF(lb .EQ. lab) THEN
                    dpbx=-xd
                    dpby=-yd
                    dpbz=-zd
                 ELSE IF(lb .EQ. lbb) THEN
                    dpbx=xd
                    dpby=yd
                    dpbz=zd
                 END IF
                 dpx=-dpax
                 dpy=-dpay
                 dpz=-dpaz
                 
                 count1=count1+1
                 matx(count1)=dpax
                 maty(count1)=dpay
                 matz(count1)=dpaz
                 aux3=(xd1*dpx+yd1*dpy+zd1*dpz)*aux1
                 
                 dpx=dpbx
                 dpy=dpby
                 dpz=dpbz
                 count1=count1+1
                 matx(count1)=dpbx
                 maty(count1)=dpby
                 matz(count1)=dpbz
                 
                 aux3=aux3+(xd1*dpx+yd1*dpy+zd1*dpz)*aux2
                 aux3=aux3*two
                 mat(ka,kb)=-aux3
              END IF
           END DO
        END DO
        
!!$=======================================================================
!!$---- Solve the linear system ------------------------------------------
!!$=======================================================================

        info=0
        DO ka=1,n0
           gam(ka)=-dd(ka)
        END DO
        IF(n0 .EQ. 1) THEN
           gam(1)=gam(1)/mat(1,1)
        ELSE IF(n0 .EQ. 2) THEN
           a11=mat(1,1)
           a12=mat(1,2)
           a21=mat(2,1)
           a22=mat(2,2)
           b11=a22
           b12=-a21
           b21=-a12
           b22=a11
           det=a11*a22-a21*a12
           det=1.0D0/det
           xd=b11*gam(1)+b21*gam(2)
           yd=b12*gam(1)+b22*gam(2)
           gam(1)=xd*det
           gam(2)=yd*det
           
        ELSEIF(n0 .EQ. 3) THEN
           a11=mat(1,1)
           a12=mat(1,2)
           a13=mat(1,3)
           a21=mat(2,1)
           a22=mat(2,2)
           a23=mat(2,3)
           a31=mat(3,1)
           a32=mat(3,2)
           a33=mat(3,3)
           b11= a22*a33-a32*a23
           b12=-a21*a33+a31*a23
           b13= a21*a32-a22*a31
           b21=-a12*a33+a32*a13
           b22= a11*a33-a31*a13
           b23=-a11*a32+a31*a12
           b31= a12*a23-a22*a13
           b32=-a11*a23+a21*a13
           b33= a11*a22-a21*a12
           
           det=a11*b11+a21*b21+a31*b31
           det=1.0D0/det
           xd=b11*gam(1)+b21*gam(2)+b31*gam(3)
           yd=b12*gam(1)+b22*gam(2)+b32*gam(3)
           zd=b13*gam(1)+b23*gam(2)+b33*gam(3)
           gam(1)=xd*det
           gam(2)=yd*det
           gam(3)=zd*det
        ELSE
           CALL dgesv(n0,1,mat,nmat,ipiv,gam,nmat,info)
           IF(info .NE. 0) THEN
              errmsg_f='In Rattle__Correct: While constraining with MIM&
                   &: matrix inversion has failed.' 
              CALL Add_Errors(-1,errmsg_f)
              out=.FALSE.
              RETURN
           END IF
        END IF
        
!!$=======================================================================
!!$---- Compute corrected positions and velocities -----------------------
!!$=======================================================================

        count1=0
        DO ka=1,n0
           la=cnstp(1,ka)
           lb=cnstp(2,ka)
           xd=mass0(la)*two
           yd=mass0(lb)*two
           DO kb=1,n0
              count1=count1+1
              IF(mask(la)) THEN
                 aux=gam(kb)
                 aux1=matx(count1)*aux
                 aux2=maty(count1)*aux
                 aux3=matz(count1)*aux
                 vpx(la)=vpx(la)+aux1*xd
                 vpy(la)=vpy(la)+aux2*xd
                 vpz(la)=vpz(la)+aux3*xd
              END IF
              count1=count1+1
              IF(mask(lb)) THEN
                 aux=gam(kb)
                 aux1=matx(count1)*aux
                 aux2=maty(count1)*aux
                 aux3=matz(count1)*aux
                 vpx(lb)=vpx(lb)+aux1*yd
                 vpy(lb)=vpy(lb)+aux2*yd
                 vpz(lb)=vpz(lb)+aux3*yd
              END IF
           END DO
           mask(la)=.FALSE.
           mask(lb)=.FALSE.
        END DO
     END IF
  END DO
  CALL Scatter_Atoms
CONTAINS
  SUBROUTINE Gather_Atoms
    INTEGER :: n,nn
    DO nn=1,natom
       n=IndBox_a_p(nn)
       xp0(nn)=xp0a(n)
       yp0(nn)=yp0a(n)
       zp0(nn)=zp0a(n)
       vpx(nn)=vpxa(n)
       vpy(nn)=vpya(n)
       vpz(nn)=vpza(n)
    END DO
  END SUBROUTINE Gather_Atoms
  SUBROUTINE Scatter_Atoms
    INTEGER :: n,nn
    DO nn=1,natom
       n=IndBox_a_p(nn)
       xp0a(n)=xp0(nn)
       yp0a(n)=yp0(nn)
       zp0a(n)=zp0(nn)
       vpxa(n)=vpx(nn)
       vpya(n)=vpy(nn)
       vpza(n)=vpz(nn)
    END DO
  END SUBROUTINE Scatter_Atoms
END FUNCTION Correct_
