      SUBROUTINE ferrf(ss_index,alphal,charge,fudge,x0,y0,z0,iz
     &     ,list,nlist,sp,fscnstr_slt,fscnstr_slv,fpx,fpy,fpz,erf_corr
     &     ,erf_arr_corr,delew,rlew)

************************************************************************
*                                                                      *
*                                                                      *
*     Compute the Ewald correction deriving from intramolecular        *
*     interaction.                                                     *
*                                                                      *
*---- Last update 06/12/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNALS  NONE                                                  *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER iz,nlist,list(2,*),ss_index(*),sp(*)
      REAL*8  charge(*),x0(*),y0(*),z0(*),fpx(*),fpy(*),fpz(*)
     &     ,fscnstr_slt,fscnstr_slv,erf_arr_corr(4,*),delew,rlew 
      REAL*8  alphal,fudge
      LOGICAL erf_corr

*------------------ VARIABLES IN COMMON --------------------------------

      INCLUDE 'unit.h'

*-------------------- LOCAL VARIABLES ----------------------------------
      INTEGER i,ia,ib,ic,type
      INTEGER ii,mma,iam
      REAL*8  furpar,qforce,alphar,twrtpi,gsrtal,xab,yab,zab,rsq,rsp
      REAL*8  a1,a2,a3,a4,a5,qp,qt,expcst,erfcst,erfst,fsrtal(2),h,c1,c2
     &     ,c3,c4,qq,corr,dcorr
      REAL*8 fact,half,factor
      DATA half/0.5D0/
      DATA a1,a2,a3/0.2548296d0,-0.28449674d0,1.4214137d0/
      DATA a4,a5/-1.453152d0,1.0614054d0/
      DATA qp/0.3275911d0/
      FACTOR(i)=(1.0D0+DFLOAT(i/IABS(i)))*half

*==================== EXECUTABLE STATEMENTS ============================


*=======================================================================
*---- Compute the intramolecular term ----------------------------------
*=======================================================================

      twrtpi=2.0d0/DSQRT(pi)
      fsrtal(1)=0.0D0
      fsrtal(2)=0.0D0
      mma=sp(1)
      IF(iz.EQ.1) THEN
         IF(erf_corr) THEN
            DO ii=1,mma
               iam=sp(ii+1)
               i=IABS(iam)
               fact=FACTOR(iam)
               ia=list(1,i)
               ib=list(2,i)

*-----------------------------------------------------------------------
*------------ No bonded interaction exists between solvent and solute --
*-----------------------------------------------------------------------
            
               type=ss_index(ia)
               qq = charge(ia)*charge(ib)
               furpar=fudge*qq
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
               corr = qq*(c1+h*(c2+h*(c3+h*c4/3)/2))
               dcorr = qq*(c2+h*(c3+h*c4/2) )/rsp
               alphar=alphal*rsp
               qt=1.0D0/(1.0e0+qp*alphar)
               expcst=exp(-alphar*alphar)
               erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*
     x              expcst
               erfst=1.0D0-erfcst
               qforce= -dcorr*fudge - furpar*(erfst-twrtpi*alphar*expcst
     &              )/(rsp*rsq)
               fsrtal(type)=fsrtal(type)-fact*(furpar*erfst/rsp-corr
     &              *fudge)
               fpx(ia)=fpx(ia)+qforce*xab 
               fpy(ia)=fpy(ia)+qforce*yab
               fpz(ia)=fpz(ia)+qforce*zab
               fpx(ib)=fpx(ib)-qforce*xab
               fpy(ib)=fpy(ib)-qforce*yab
               fpz(ib)=fpz(ib)-qforce*zab
            END DO
         ELSE
            DO ii=1,mma
               iam=sp(ii+1)
               i=IABS(iam)
               fact=FACTOR(iam)
               ia=list(1,i)
               ib=list(2,i)
               
*-----------------------------------------------------------------------
*------------ No bonded interaction exists between solvent and solute --
*-----------------------------------------------------------------------
               
               type=ss_index(ia)
               furpar=fudge*charge(ia)*charge(ib)
               xab=x0(ia)-x0(ib)
               yab=y0(ia)-y0(ib)
               zab=z0(ia)-z0(ib)
               rsq=xab*xab+yab*yab+zab*zab
               rsp=DSQRT(rsq)
               alphar=alphal*rsp
               qt=1.0D0/(1.0e0+qp*alphar)
               expcst=exp(-alphar*alphar)
               erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*
     x              expcst
               erfst=1.0D0-erfcst
               qforce=-furpar*(erfst-twrtpi*alphar*expcst)/(rsp*rsq)
               fsrtal(type)=fsrtal(type)-fact*furpar*erfst/rsp
               fpx(ia)=fpx(ia)+qforce*xab
               fpy(ia)=fpy(ia)+qforce*yab
               fpz(ia)=fpz(ia)+qforce*zab
               fpx(ib)=fpx(ib)-qforce*xab
               fpy(ib)=fpy(ib)-qforce*yab
               fpz(ib)=fpz(ib)-qforce*zab
            END DO
         END IF
      ELSE
         IF(erf_corr) THEN
            DO ii=1,mma
               iam=sp(ii+1)
               i=IABS(iam)
               fact=FACTOR(iam)
               ia=list(1,i)
               ib=list(2,i)
               
*-----------------------------------------------------------------------
*------------ No bonded interaction exists between solvent and solute --
*-----------------------------------------------------------------------
               
               type=ss_index(ia)
               furpar=charge(ia)*charge(ib)
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
               corr = (c1+h*(c2+h*(c3+h*c4/3)/2))*furpar
               alphar=alphal*rsp
               qt=1.0e0/(1.0e0+qp*alphar)
               expcst=exp(-alphar*alphar)
               erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*expcst
               erfst=1.0D0-erfcst
               fsrtal(type)=fsrtal(type)-fact*(furpar*erfst/rsp
     &              -corr)
            END DO
         ELSE
            DO ii=1,mma
               iam=sp(ii+1)
               i=IABS(iam)
               fact=FACTOR(iam)
               ia=list(1,i)
               ib=list(2,i)
               
*-----------------------------------------------------------------------
*------------ No bonded interaction exists between solvent and solute --
*-----------------------------------------------------------------------
               
               type=ss_index(ia)
               furpar=charge(ia)*charge(ib)*fudge
               xab=x0(ia)-x0(ib)
               yab=y0(ia)-y0(ib)
               zab=z0(ia)-z0(ib)
               rsq=xab*xab+yab*yab+zab*zab
               rsp=DSQRT(rsq)
               alphar=alphal*rsp
               qt=1.0e0/(1.0e0+qp*alphar)
               expcst=exp(-alphar*alphar)
               erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*
     x              expcst
               erfst=1.0D0-erfcst
               fsrtal(type)=fsrtal(type)-fact*furpar*erfst/
     x              rsp
            END DO
            
         END IF
      END IF
      fscnstr_slt=fsrtal(1)
      fscnstr_slv=fsrtal(2)

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
