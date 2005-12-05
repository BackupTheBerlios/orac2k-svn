      SUBROUTINE rattle_verlet(nstart,nend,dt,xp1,yp1,zp1,xp0,yp0
     &     ,zp0,vpx,vpy,vpz,nato,cnst,dss,coeff,ncnstr,mass,dnit,nprot
     &     ,cnst_protl,mim_lim,gcpu,iret,errmsg)

************************************************************************
*   Time-stamp: <01/03/13 22:45:07 marchi>                             *
*                                                                      *
*     It will supply a new set of coordinates for which                *
*     the set of constraints are satisfied.                            *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Feb 13 1997 -                                     *
*                                                                      *
************************************************************************


*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER ncnstr,nato,iret,nprot,mim_lim,nstart,nend
      INTEGER cnst(2,*),cnst_protl(*)
      REAL*8  xp1(*),yp1(*),zp1(*), xp0(*),yp0(*),zp0(*),vpx(*)
     &     ,vpy(*),vpz(*),dss(*),coeff(*),mass(*),dnit,dt,gcpu
      CHARACTER*80 errmsg

*----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'parst.h'
      INTEGER nmat
      PARAMETER (nmat = 4)
      REAL*8  x0(m9),y0(m9),z0(m9),mass0(m1)
      REAL*8  mat(nmat,nmat),gam(nmat),gamo(nmat),matx(2*nmat*nmat)
     &     ,maty(2*nmat*nmat),matz(2*nmat*nmat),xc(nmat),yc(nmat)
     &     ,zc(nmat),dd(nmat),aux,xd,yd,zd
      LOGICAL mask(m1)
      INTEGER ipiv(2*nmat),info

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER la,lb,k,iter,iox,i,cnstp,iter1,count,ka,lab,lbb,k1,k2,kb
     &     ,count1
      REAL*8 tol,xab,yab,zab,dpx,dpy,dpz,dpp,dps,gg,amsla,amslb,dpax
     &     ,dpay,dpaz,dpbx,dpby,dpbz
      REAL*8  gg1,gg2,dcnt,dti,xk,yk,zk,tol_mim,two,aux1,aux2,aux3,det
     &     ,zero,four,a11,a12,a13,a21,a22,a23,a31,a32,a33,b11,b12,b13
     &     ,b21,b22,b23,b31,b32,b33,gcpu_mm,vfcp_mm,tfcp_mm,tdelta_mm
     &     ,elapse
      COMMON /rag1/ xab,yab,zab,gg1,gg2,dcnt,dti,xk,yk,zk,aux1,aux2,aux3
     &     ,det,a11,a12,a13,a21,a22,a23,a31,a32,a33,b11,b12,b13,b21,b22
     &     ,b23,b31,b32,b33,dpx,dpy,dpz,dpp,dps,gg,amsla,amslb,dpax,dpay
     &     ,dpaz,dpbx,dpby,dpbz,x0,y0,z0,mass0,mat,gam,gamo,matx,maty
     &     ,matz,ipiv,mask
      DATA tol/1.0D-7/tol_mim/1.0D-5/two/2.0D0/zero/0.0D0/four/4.0D0/

*==================== EXECUTABLE STATEMENTS ============================

      IF(ncnstr .EQ. 0 .OR. nprot .EQ. 0) RETURN
      IF(nmat .LT. mim_lim) THEN
         errmsg=' While SHAKEing dimension of the cnstr matrix exceeds'
     &        //' physical dimanesions.'
         iret=1
         RETURN
      END IF

      IF(m9 .LT. ncnstr) THEN
          errmsg=' In CNSTRP. Physical dimensions of working arrays'//
     x ' are insufficient. Abort!'
          iret=1
          WRITE(6,'('' M9 = '',i6,'' MCNSTR = '',i6)') m9,ncnstr
          RETURN
      END IF

      CALL timer(vfcp_mm,tfcp_mm,elapse)
      tdelta_mm=tfcp_mm

*=======================================================================
*---- Copy the coordinates of the molecule IA to a temporary array -----
*=======================================================================

      dti=1.0D0/dt
      DO i=nstart,nend
         mass0(i)=1.0D0/mass(i)
         mask(i)=.TRUE.
      END DO

*=======================================================================
*---- Compute the vectors associated with each constraint and store ----
*---- them in a temporary array ----------------------------------------
*=======================================================================

      count=0
      DO i=1,nprot
         cnstp=cnst_protl(1+count)
         DO ka=1,cnstp
            k=cnst_protl(1+count+ka)
            la=cnst(1,k)
            lb=cnst(2,k)
            xab=xp0(la)-xp0(lb)
            yab=yp0(la)-yp0(lb)
            zab=zp0(la)-zp0(lb)
            x0(k)=xab
            y0(k)=yab
            z0(k)=zab
         END DO
         count=count+cnstp+1
      END DO

*=======================================================================
*---- SHAKE loop -------------------------------------------------------
*=======================================================================

      iter1=0
      count=0
      DO i=1,nprot
         cnstp=cnst_protl(1+count)
         IF(cnstp .GT. mim_lim) THEN
            iter=0
1000        CONTINUE
            iox=0
            DO ka=1,cnstp
               k=cnst_protl(1+count+ka)
               la=cnst(1,k)
               lb=cnst(2,k)
               dpx=xp1(la)-xp1(lb)
               dpy=yp1(la)-yp1(lb)
               dpz=zp1(la)-zp1(lb)
               dpp=dpx**2+dpy**2+dpz**2
               dcnt=dpp-dss(k)
               IF(DABS(dcnt).GE.tol) THEN
                  iox=1
                  xk=x0(k)
                  yk=y0(k)
                  zk=z0(k)
                  amsla=mass0(la)
                  amslb=mass0(lb)
                  dps=dpx*xk+dpy*yk+dpz*zk
                  gg=dcnt/(dps*coeff(k))
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
                  iret=1
                  errmsg=
     & ' While SHAKEing : The iteration procedure did not'
     x              //' converge.'
                  RETURN
               END IF
               GOTO 1000
            END IF
         ELSE

*=======================================================================
*---  Use Matrix Inversion Method for Constraints ---------------------- 
*=======================================================================


            DO ka=1,cnstp
               gam(ka)=0.0D0
            END DO

*=========== Build constraint matrix ====================================

            DO ka=1,cnstp
               k1=cnst_protl(1+count+ka)
               la=cnst(1,k1)
               lb=cnst(2,k1)
               xc(ka)=xp1(la)-xp1(lb)
               yc(ka)=yp1(la)-yp1(lb)
               zc(ka)=zp1(la)-zp1(lb)
            END DO
               
*======== First loop on constraints k1 =================================


            count1=0
            DO ka=1,cnstp
               k1=cnst_protl(1+count+ka)
               la=cnst(1,k1)
               lb=cnst(2,k1)
               dd(ka)=dss(k1)-(xc(ka)*xc(ka)+yc(ka)*yc(ka)+zc(ka)
     &              *zc(ka))
               
*============= Second loop on constraints k2 ===========================
                     
               DO kb=1,cnstp
                  k2=cnst_protl(1+count+kb)
                  lab=cnst(1,k2)
                  lbb=cnst(2,k2)
                  
                  xd=x0(k2)
                  yd=y0(k2)
                  zd=z0(k2)
                  
                  aux1=mass0(la)
                  aux2=mass0(lb)
                  
                  IF(la .NE. lab .AND. la .NE. lbb .AND. lb .NE.
     &                 lab.AND. lb .NE. lbb) THEN
                     mat(ka,kb)=0.0D0
                  ELSE
                     dpax=0.0D0
                     dpay=0.0D0
                     dpaz=0.0D0
                     dpbx=0.0D0
                     dpby=0.0D0
                     dpbz=0.0D0
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

*=======================================================================
*---- Iteration to solve the non linear system -------------------------
*=======================================================================

            iox=1
            iter=0
            info=0
            DO WHILE(iox .EQ. 1)
               DO ka=1,cnstp
                  gamo(ka)=gam(ka)
               END DO
               
               IF(iter .EQ. 0) THEN
                  DO ka=1,cnstp
                     gam(ka)=dd(ka)
                  END DO
                  IF(cnstp .EQ. 1) THEN
                     gam(1)=gam(1)/mat(1,1)
                  ELSEIF(cnstp .EQ. 2) THEN
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
                  ELSE IF(cnstp .EQ. 3) THEN
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
                     CALL dgesv(cnstp,1,mat,nmat,ipiv,gam,nmat,info)
                  END IF
               ELSE
                  count1=0
                  DO ka=1,cnstp
                     xd=0.0D0
                     yd=0.0D0
                     zd=0.0D0
                     gg1=mass0(la)*two
                     gg2=mass0(lb)*two
                     DO kb=1,cnstp
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

                  DO ka=1,cnstp
                     gam(ka)=xc(ka)
                  END DO

                  IF(cnstp .EQ. 1) THEN
                     gam(1)=gam(1)/mat(1,1)
                  ELSE IF(cnstp .EQ. 2) THEN
                     xd=b11*gam(1)+b21*gam(2)
                     yd=b12*gam(1)+b22*gam(2)
                     gam(1)=xd*det
                     gam(2)=yd*det
                  ELSE IF(cnstp .EQ. 3) THEN
                     xd=b11*gam(1)+b21*gam(2)+b31*gam(3)
                     yd=b12*gam(1)+b22*gam(2)+b32*gam(3)
                     zd=b13*gam(1)+b23*gam(2)+b33*gam(3)
                     gam(1)=xd*det
                     gam(2)=yd*det
                     gam(3)=zd*det
                  ELSE
                     CALL dgetrs('N',cnstp,1,mat,nmat,ipiv,gam,nmat,info
     &                    )
                  END IF
               END IF
               IF(info .NE. 0) THEN
                  iret=1
                  errmsg=
     & ' While constraining with MIM: matrix inversion failed. '
                  RETURN
               END IF
               iter=iter+1
               IF(iter.GT.5000)THEN
                  iret=1
                  errmsg=
     & ' While constraining with MIM: The iteration procedure did not'
     x              //' converge.'
                  RETURN
               END IF
               iox=0
               DO ka=1,cnstp
                  IF(gam(ka) .NE. 0) THEN
                     aux=(gam(ka)-gamo(ka))/gam(ka)
                     aux=DABS(aux)
                     IF(aux .GT. tol_mim) iox=1
                  END IF
               END DO
            END DO

            iter1=iter1+iter

*=======================================================================
*---- Compute corrected position and velocities ------------------------
*=======================================================================

            count1=0
            DO ka=1,cnstp
               k1=cnst_protl(1+count+ka)
               la=cnst(1,k1)
               lb=cnst(2,k1)
               xd=mass0(la)*two
               yd=mass0(lb)*two
               DO kb=1,cnstp
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
         count=count+cnstp+1
      END DO

      dnit=dnit+DFLOAT(iter1)/DFLOAT(nprot)
      CALL timer(vfcp_mm,tfcp_mm,elapse)
      tdelta_mm=tfcp_mm-tdelta_mm
      gcpu=gcpu+tdelta_mm

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
