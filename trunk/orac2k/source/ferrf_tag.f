      SUBROUTINE ferrf_tag(ss_index,alphal,charge,x0,y0,z0,list
     &     ,nlist,sp,fscnstr_slt,fscnstr_slv,fpx,fpy,fpz,phi,ma,tags
     &     ,erf_corr,erf_arr_corr,delew,rlew)

************************************************************************
*   Time-stamp: <2007-11-27 16:32:11 marchi>                             *
*                                                                      *
*                                                                      *
*     Compute the Ewald correction from intramolecular interactions    *
*     using the tag vector                                             *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Jan 10 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nlist,ma,list(2,*),ss_index(*),tags(*),sp(*)
      REAL*8  charge(*),x0(*),y0(*),z0(*),fpx(ma,*),fpy(ma,*),fpz(ma,*)
     &     ,phi(*),fscnstr_slt,fscnstr_slv,erf_arr_corr(4,*),delew,rlew 
      REAL*8  alphal
      LOGICAL erf_corr

*----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'unit.h'

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,ia,ib,type,p1,ic
      INTEGER ii,mma,iam
      REAL*8  furpar,qforce,alphar,twrtpi,gsrtal,xab,yab,zab,rsq,rsp
      REAL*8  a1,a2,a3,a4,a5,qp,qt,expcst,erfcst,erfst,fsrtal(2),h,c1,c2
     &     ,c3,c4,qq,corr,dcorr
      REAL*8 fact,half,factor,furpar1,furpar2,corr1,corr2
      DATA half/0.5D0/
      DATA a1,a2,a3/0.2548296d0,-0.28449674d0,1.4214137d0/
      DATA a4,a5/-1.453152d0,1.0614054d0/
      DATA qp/0.3275911d0/
      FACTOR(i)=(1.0D0+DFLOAT(i/IABS(i)))*half

*----------------------- EXECUTABLE STATEMENTS ------------------------*

*=======================================================================
*---- Compute the intramolecular term ----------------------------------
*=======================================================================

      twrtpi=2.0d0/DSQRT(pi)
      fsrtal(1)=0.0D0
      fsrtal(2)=0.0D0
      mma=sp(1)
      IF(erf_corr) THEN
         DO ii=1,mma
            iam=sp(ii+1)
            i=IABS(iam)
            fact=FACTOR(iam)
            ia=list(1,i)
            ib=list(2,i)
            p1=tags(i)
            
*-----------------------------------------------------------------------
*------------ No bonded interaction exists between solvent and solute --
*-----------------------------------------------------------------------
            
            type=ss_index(ia)
            furpar=charge(ia)*charge(ib)
            furpar1=charge(ib)
            furpar2=charge(ia)
            qq=furpar
            
            xab=x0(ia)-x0(ib)
            yab=y0(ia)-y0(ib)
            zab=z0(ia)-z0(ib)
            rsq=xab*xab+yab*yab+zab*zab
            rsp=DSQRT(rsq)
            ic = 1+ (rsp-rlew)/delew  
            h = rsp-(rlew+(ic-1)*delew)
            c1=erf_arr_corr(1,ic) 
            c2=erf_arr_corr(2,ic) 
            c3=erf_arr_corr(3,ic) 
            c4=erf_arr_corr(4,ic) 
            corr = (c1+h*(c2+h*(c3+h*c4/3)/2))
            corr1=charge(ib)*corr
            corr2=charge(ia)*corr
            corr=qq*corr

            dcorr = qq*(c2+h*(c3+h*c4/2) )/rsp
            alphar=alphal*rsp
            qt=1.0D0/(1.0e0+qp*alphar)
            expcst=exp(-alphar*alphar)
            erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*
     x           expcst
            erfst=1.0D0-erfcst
            qforce=-dcorr-furpar*(erfst-twrtpi*alphar*expcst)/(rsp*rsq)
            fsrtal(type)=fsrtal(type)-fact*(furpar*erfst/rsp-corr)
            phi(ia)=phi(ia)-fact*(furpar1*erfst/rsp-corr1)
            phi(ib)=phi(ib)-fact*(furpar2*erfst/rsp-corr2)
            fpx(ia,p1)=fpx(ia,p1)+qforce*xab
            fpy(ia,p1)=fpy(ia,p1)+qforce*yab
            fpz(ia,p1)=fpz(ia,p1)+qforce*zab
            fpx(ib,p1)=fpx(ib,p1)-qforce*xab
            fpy(ib,p1)=fpy(ib,p1)-qforce*yab
            fpz(ib,p1)=fpz(ib,p1)-qforce*zab
         END DO
      ELSE
         DO ii=1,mma
            iam=sp(ii+1)
            i=IABS(iam)
            fact=FACTOR(iam)
            ia=list(1,i)
            ib=list(2,i)
            p1=tags(i)
            
*-----------------------------------------------------------------------
*------------ No bonded interaction exists between solvent and solute --
*-----------------------------------------------------------------------
            
            type=ss_index(ia)
            furpar=charge(ia)*charge(ib)
            furpar1=charge(ib)
            furpar2=charge(ia)

            xab=x0(ia)-x0(ib)
            yab=y0(ia)-y0(ib)
            zab=z0(ia)-z0(ib)
            rsq=xab*xab+yab*yab+zab*zab
            rsp=DSQRT(rsq)
            alphar=alphal*rsp
            qt=1.0D0/(1.0e0+qp*alphar)
            expcst=exp(-alphar*alphar)
            erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*
     x           expcst
            erfst=1.0D0-erfcst
            qforce=-furpar*(erfst-twrtpi*alphar*expcst)/(rsp*rsq)
            fsrtal(type)=fsrtal(type)-fact*furpar*erfst/rsp
            phi(ia)=phi(ia)-fact*furpar1*erfst/rsp
            phi(ib)=phi(ib)-fact*furpar2*erfst/rsp

            fpx(ia,p1)=fpx(ia,p1)+qforce*xab
            fpy(ia,p1)=fpy(ia,p1)+qforce*yab
            fpz(ia,p1)=fpz(ia,p1)+qforce*zab
            fpx(ib,p1)=fpx(ib,p1)-qforce*xab
            fpy(ib,p1)=fpy(ib,p1)-qforce*yab
            fpz(ib,p1)=fpz(ib,p1)-qforce*zab
         END DO
      END IF
      fscnstr_slt=fsrtal(1)
      fscnstr_slv=fsrtal(2)
      
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*
      
      RETURN
      END
