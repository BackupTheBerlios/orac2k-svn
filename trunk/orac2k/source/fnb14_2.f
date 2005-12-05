      SUBROUTINE fnb14_2(ss_index,xp0,yp0,zp0,charge,nato
     &     ,alphal,int14,int14p,type14,sp,fudge,fpx
     &     ,fpy,fpz,ucoul_slt,ucoul_slv)

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER nato,int14p,int14(2,*),type14(*),ss_index(*),sp(*)
      REAL*8  xp0(*),yp0(*),zp0(*),fpx(*),fpy(*),
     x        fpz(*),alphal,charge(*),fudge,ucoul_slt,ucoul_slv

*------------------- VARIABLES IN COMMON -------------------------------

      INCLUDE 'unit.h'

*-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8 xpi,ypi,zpi,xc,yc,zc,rsq,rsp,rsqi,qforce
      REAL*8 xpj,ypj,zpj,chrgei,chrgej,aux1
      REAL*8 a1,a2,a3,a4,a5,qp,qt,expcst,erfcst
      REAL*8 rspqi,alphar,furpar,twrtpi,ucoul(2)
      INTEGER i,i1,j,k,l,lij,type
      INTEGER ii,mma,iam
      REAL*8 fact,quarter,half,factor
      DATA a1,a2,a3/0.2548296d0,-0.28449674d0,1.4214137d0/
      DATA a4,a5/-1.453152d0,1.0614054d0/
      DATA qp/0.3275911d0/
      DATA quarter,half/0.25D0,0.5D0/
      FACTOR(i)=(1.0D0+DFLOAT(i/IABS(i)))*half


*==================== EXECUTABLE STATEMENTS ============================

      DO i=1,2
         ucoul(i)=0.0D+0
      END DO

      ucoul_slt=0.0D0
      ucoul_slv=0.0D0

      IF(int14p.EQ.0) RETURN
      twrtpi=2.0d0/DSQRT(pi)

*=======================================================================
*----- Take Ewald for the electrostatic interaction --------------------
*=======================================================================

      mma=sp(1)
      DO ii=1,mma
         iam=sp(1+ii)
         i=IABS(iam)
         fact=FACTOR(iam)
         i1=int14(1,i)
         j=int14(2,i)

*-----------------------------------------------------------------------
*------------ No bonded interaction exists between solvent and solute --
*-----------------------------------------------------------------------
         
         type=ss_index(i1)
         lij=type14(i)
         xpi=xp0(i1)
         ypi=yp0(i1)
         zpi=zp0(i1)
         xpj=xp0(j)
         ypj=yp0(j)
         zpj=zp0(j)
         chrgei=charge(i1)
         chrgej=charge(j)
         xc=xpi-xpj
         yc=ypi-ypj
         zc=zpi-zpj
         rsq=xc**2+yc**2+zc**2
         rsp=DSQRT(rsq)
         rsqi=1.0d0/rsq
         rspqi=rsqi/rsp
         furpar=fudge*chrgei*chrgej

         alphar=alphal*rsp
         qt=1.0d0/(1.0d0+qp*alphar)
         expcst=dexp(-alphar*alphar)
         erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*expcst
         aux1=furpar*erfcst/rsp
         ucoul(type)=ucoul(type)+aux1*fact

         qforce=furpar*(erfcst+twrtpi*alphar*expcst)*rspqi
         fpx(i1)=fpx(i1)+qforce*xc
         fpy(i1)=fpy(i1)+qforce*yc
         fpz(i1)=fpz(i1)+qforce*zc
         fpx(j)=fpx(j)-qforce*xc
         fpy(j)=fpy(j)-qforce*yc
         fpz(j)=fpz(j)-qforce*zc

      END DO

      ucoul_slt=ucoul(1)
      ucoul_slv=ucoul(2)

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
