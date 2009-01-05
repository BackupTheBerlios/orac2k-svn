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
FUNCTION Verlet_(dt,xp0a,yp0a,zp0a,vpxa,vpya,vpza) RESULT(out)
  LOGICAL :: out
  REAL(8) :: dt,xp0a(:),yp0a(:),zp0a(:),vpxa(:),vpya(:),vpza(:)
  INTEGER, POINTER :: cnstp(:,:)
  REAL(8), POINTER :: dssp(:),coeffp(:),x0(:),y0(:),z0(:)
  TYPE(Rattle__Type2), ALLOCATABLE, TARGET :: xx0(:),yy0(:),zz0(:)
  INTEGER :: nn,m,n,n0
  INTEGER :: la,lb,k,iter,iox,i,iter1,count,ka,lab,lbb,k1,k2&
       &,kb,count1
  REAL(8) ::  xab,yab,zab,dpx,dpy,dpz,dpp,dps,gg,amsla,amslb&
       &,dpax,dpay,dpaz,dpbx,dpby,dpbz
  REAL(8) ::  gg1,gg2,dcnt,dti,xk,yk,zk,aux1,aux2,aux3&
       &,det,a11,a12,a13,a21,a22,a23,a31,a32,a33,b11,b12&
       &,b13,b21,b22,b23,b31,b32,b33,gcpu_mm,vfcp_mm,tfcp_mm&
       &,tdelta_mm,elapse
  
  REAL(8), SAVE :: tol=1.0D-7,tol_mim=1.0D-5,zero=0.0_8,one=1.0_8&
       &,two=2.0_8,three=3.0_8,four=4.0_8
  REAL(8) ::   aux,xd,yd,zd
  INTEGER :: info
  LOGICAL, ALLOCATABLE :: mask(:)

  out=.TRUE.
  IF(.NOT. Rattle__Param % switch) RETURN  
  dti=1.0D0/dt
  
  CALL Gather_Atoms
  
  ALLOCATE(xx0(nc),yy0(nc),zz0(nc),mask(natom))
  
  DO n=1,nc
     cnstp=>cnst(n) % n1
     n0=SIZE(cnstp,2)
     ALLOCATE(xx0(n) % g(n0))
     ALLOCATE(yy0(n) % g(n0))
     ALLOCATE(zz0(n) % g(n0))
  END DO
  
  DO nn=1,nc
     cnstp=>cnst(nn) % n1
     x0=>xx0(nn) % g
     y0=>yy0(nn) % g
     z0=>zz0(nn) % g
     n0=SIZE(cnstp,2)
     DO m=1,n0
        la=cnstp(1,m)
        lb=cnstp(2,m)
        x0(m)=xp0(la)-xp0(lb)
        y0(m)=yp0(la)-yp0(lb)
        z0(m)=zp0(la)-zp0(lb)
     END DO
  END DO
  
!!$
!!$---- SHAKE loop -------------------------------------------------------
!!$
  
  iter1=0
  count=0
  DO nn=1,nc
     cnstp=>cnst(nn) % n1
     dssp=>dss(nn) % g
     coeffp=>coeff(nn) % g
     x0=>xx0(nn) % g
     y0=>yy0(nn) % g
     z0=>zz0(nn) % g
     n0=SIZE(cnstp,2)
     IF(n0 > Rattle__Param % mim_Max) THEN
        iter=0
1000    CONTINUE
        iox=0
        DO k=1,n0
           la=cnstp(1,k)
           lb=cnstp(2,k)
           dpx=xp1(la)-xp1(lb)
           dpy=yp1(la)-yp1(lb)
           dpz=zp1(la)-zp1(lb)
           dpp=dpx**2+dpy**2+dpz**2
           dcnt=dpp-dssp(k)
           IF(DABS(dcnt) >= tol) THEN
              iox=1
              xk=x0(k)
              yk=y0(k)
              zk=z0(k)
              amsla=mass0(la)
              amslb=mass0(lb)
              dps=dpx*xk+dpy*yk+dpz*zk
              gg=dcnt/(dps*coeffp(k))
              gg1=-gg*amsla
              gg2=gg*amslb
              xp1(la)=xp1(la)+xk*gg1
              xp1(lb)=xp1(lb)+xk*gg2
              
              yp1(la)=yp1(la)+yk*gg1
              yp1(lb)=yp1(lb)+yk*gg2
              
              zp1(la)=zp1(la)+zk*gg1
              zp1(lb)=zp1(lb)+zk*gg2
              
              gg1=gg1*dti
              gg2=gg2*dti
              
              vpx(la)=vpx(la)+xk*gg1
              vpx(lb)=vpx(lb)+xk*gg2
              
              vpy(la)=vpy(la)+yk*gg1
              vpy(lb)=vpy(lb)+yk*gg2
              
              vpz(la)=vpz(la)+zk*gg1
              vpz(lb)=vpz(lb)+zk*gg2
           END IF
        END DO
        IF(iox.ne.0) THEN
           iter1=iter1+1
           iter=iter+1
           IF(iter.GT.5000)THEN
              errmsg_f=' While SHAKEing : The iteration procedure &
                   &did not converge.'
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
        
        DO ka=1,n0
           gam(ka)=0.0D0
        END DO
        
!!$=========== Build constraint matrix ====================================
        
        DO ka=1,n0
           la=cnstp(1,ka)
           lb=cnstp(2,ka)
           xc(ka)=xp1(la)-xp1(lb)
           yc(ka)=yp1(la)-yp1(lb)
           zc(ka)=zp1(la)-zp1(lb)
        END DO
        
!!$======== First loop on constraints ka =================================
        
        
        count1=0
        DO ka=1,n0
           la=cnstp(1,ka)
           lb=cnstp(2,ka)
           dd(ka)=dssp(ka)-(xc(ka)*xc(ka)+yc(ka)*yc(ka)+zc(ka)*zc(ka))
           
!!$============= Second loop on constraints kb ===========================
           
           DO kb=1,n0
              lab=cnstp(1,kb)
              lbb=cnstp(2,kb)
              
              xd=x0(kb)
              yd=y0(kb)
              zd=z0(kb)
              
              aux1=mass0(la)
              aux2=mass0(lb)
              
              IF(la /= lab .AND. la /= lbb .AND. lb /=&
                   & lab .AND. lb /= lbb) THEN 
                 mat(ka,kb)=0.0D0
              ELSE
                 dpax=0.0D0
                 dpay=0.0D0
                 dpaz=0.0D0
                 dpbx=0.0D0
                 dpby=0.0D0
                 dpbz=0.0D0
                 IF(la == lab) THEN
                    dpax=-xd
                    dpay=-yd
                    dpaz=-zd
                 ELSE IF(la == lbb) THEN
                    dpax=xd
                    dpay=yd
                    dpaz=zd
                 END IF
                 
                 IF(lb == lab) THEN
                    dpbx=-xd
                    dpby=-yd
                    dpbz=-zd
                 ELSE IF(lb == lbb) THEN
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
                 
                 aux3=(dpx*xc(ka)+dpy*yc(ka)+dpz*zc(ka))*aux1
                 
                 dpx=dpbx
                 dpy=dpby
                 dpz=dpbz
                 count1=count1+1
                 matx(count1)=dpbx
                 maty(count1)=dpby
                 matz(count1)=dpbz
                 
                 aux3=aux3+(dpx*xc(ka)+dpy*yc(ka)+dpz*zc(ka))*aux2
                 aux3=aux3*four
                 
                 mat(ka,kb)=-aux3
                 
              END IF
           END DO
        END DO
        
!!$
!!$---- Iteration to solve the non linear system -------------------------
!!$
        
        iox=1
        iter=0
        info=0
        DO WHILE(iox .EQ. 1)
           DO ka=1,n0
              gamo(ka)=gam(ka)
           END DO
           
           IF(iter .EQ. 0) THEN
              DO ka=1,n0
                 gam(ka)=dd(ka)
              END DO
              IF(n0 .EQ. 1) THEN
                 gam(1)=gam(1)/mat(1,1)
              ELSEIF(n0 .EQ. 2) THEN
                 b11=mat(2,2)
                 b12=-mat(2,1)
                 b21=-mat(1,2)
                 b22=mat(1,1)
                 det=b11*b22-b21*b12
                 det=1.0D0/det
                 xd=b11*gam(1)+b21*gam(2)
                 yd=b12*gam(1)+b22*gam(2)
                 gam(1)=xd*det
                 gam(2)=yd*det
              ELSE IF(n0 .EQ. 3) THEN
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
              END IF
           ELSE
              count1=0
              DO ka=1,n0
                 xd=0.0D0
                 yd=0.0D0
                 zd=0.0D0
                 gg1=mass0(la)*two
                 gg2=mass0(lb)*two
                 DO kb=1,n0
                    count1=count1+1
                    aux=gam(kb)*gg1
                    xd=xd+matx(count1)*aux
                    yd=yd+maty(count1)*aux
                    zd=zd+matz(count1)*aux
                    count1=count1+1
                    aux1=-matx(count1)
                    aux2=-maty(count1)
                    aux3=-matz(count1)
                    aux=gam(kb)*gg2
                    xd=xd+aux1*aux
                    yd=yd+aux2*aux
                    zd=zd+aux3*aux
                 END DO
                 xc(ka)=dd(ka)-(xd*xd+yd*yd+zd*zd)
              END DO
              
              DO ka=1,n0
                 gam(ka)=xc(ka)
              END DO
              
              IF(n0 .EQ. 1) THEN
                 gam(1)=gam(1)/mat(1,1)
              ELSE IF(n0 .EQ. 2) THEN
                 xd=b11*gam(1)+b21*gam(2)
                 yd=b12*gam(1)+b22*gam(2)
                 gam(1)=xd*det
                 gam(2)=yd*det
              ELSE IF(n0 .EQ. 3) THEN
                 xd=b11*gam(1)+b21*gam(2)+b31*gam(3)
                 yd=b12*gam(1)+b22*gam(2)+b32*gam(3)
                 zd=b13*gam(1)+b23*gam(2)+b33*gam(3)
                 gam(1)=xd*det
                 gam(2)=yd*det
                 gam(3)=zd*det
              ELSE
                 CALL dgetrs('N',n0,1,mat,nmat,ipiv,gam,nmat,info)
              END IF
           END IF
           IF(info .NE. 0) THEN
              errmsg_f=' While constraining with MIM: matrix inversion failed. '
              CALL Add_Errors(-1,errmsg_f)
              out=.FALSE.
              RETURN
           END IF
           iter=iter+1
           IF(iter.GT.5000)THEN
              errmsg_f=' While constraining with MIM: The iteration &
                   &procedure did not converge.' 
              CALL Add_Errors(-1,errmsg_f)
              out=.FALSE.
              RETURN
           END IF
           iox=0
           DO ka=1,n0
              IF(gam(ka) .NE. 0) THEN
                 aux=(gam(ka)-gamo(ka))/gam(ka)
                 aux=DABS(aux)
                 IF(aux .GT. tol_mim) iox=1
              END IF
           END DO
        END DO
        
        iter1=iter1+iter
        
!!$
!!$---- Compute corrected position and velocities ------------------------
!!$
        
        count1=0
        DO ka=1,n0
           la=cnstp(1,ka)
           lb=cnstp(2,ka)
           xd=mass0(la)*two
           yd=mass0(lb)*two
           DO kb=1,n0
              count1=count1+1
              IF(mask(la)) THEN
                 aux=gam(kb)*xd
                 aux1=matx(count1)
                 aux2=maty(count1)
                 aux3=matz(count1)
                 xp1(la)=xp1(la)+aux1*aux
                 yp1(la)=yp1(la)+aux2*aux
                 zp1(la)=zp1(la)+aux3*aux
                 aux=aux*dti
                 vpx(la)=vpx(la)+aux1*aux
                 vpy(la)=vpy(la)+aux2*aux
                 vpz(la)=vpz(la)+aux3*aux
              END IF
              count1=count1+1
              IF(mask(lb)) THEN
                 aux=gam(kb)*yd
                 aux1=matx(count1)
                 aux2=maty(count1)
                 aux3=matz(count1)
                 xp1(lb)=xp1(lb)+aux1*aux
                 yp1(lb)=yp1(lb)+aux2*aux
                 zp1(lb)=zp1(lb)+aux3*aux
                 aux=aux*dti
                 vpx(lb)=vpx(lb)+aux1*aux
                 vpy(lb)=vpy(lb)+aux2*aux
                 vpz(lb)=vpz(lb)+aux3*aux
              END IF
           END DO
           mask(la)=.FALSE.
           mask(lb)=.FALSE.
        END DO
     END IF
  END DO
  CALL Scatter_Atoms
  WRITE(*,*) 'Iter 1',DBLE(Iter1)/DBLE(nc),Rattle__Param % mim_Max
CONTAINS
  SUBROUTINE Gather_Atoms
    INTEGER :: n,nn
    DO nn=1,natom
       n=IndBox_a_p(nn)
       xp0(nn)=xp0a(n)
       yp0(nn)=yp0a(n)
       zp0(nn)=zp0a(n)
       xp1(nn)=xp0a(n)
       yp1(nn)=yp0a(n)
       zp1(nn)=zp0a(n)
       vpx(nn)=vpxa(n)
       vpy(nn)=vpya(n)
       vpz(nn)=vpza(n)
    END DO
  END SUBROUTINE Gather_Atoms
  SUBROUTINE Scatter_Atoms
    INTEGER :: n,nn
    DO nn=1,natom
       n=IndBox_a_p(nn)
       xp0a(n)=xp1(nn)
       yp0a(n)=yp1(nn)
       zp0a(n)=zp1(nn)
       vpxa(n)=vpx(nn)
       vpya(n)=vpy(nn)
       vpza(n)=vpz(nn)
    END DO
  END SUBROUTINE Scatter_Atoms
END FUNCTION Verlet_
