      SUBROUTINE fnbgrp2(ss_index,xp0,yp0,zp0,charge,ecc12,ecc6,nbtype
     &     ,type,m6,alphal,ingrp,ingrpp,sp,fpx,fpy,fpz,do_LJ,flj_x,flj_y
     &     ,flj_z,ucoul_slt,ucoul_slv,uconf_slt,uconf_slv,plrzbij,Utotal
     &     ,U_Thole,Edx,Edy,Edz,dipole)

************************************************************************
*                                                                      *
*     Compute the contribution from third neighbour (1-4) non          *
*     bonded interactions to forces and energies.                      *
*                                                                      *
*     XP0     :  Coordinates of the macromolecule.                (I)  *
*     YP0        >> real*8 XP0(NATO), YP0(NATO), ZP0(NATO) <<          *
*     ZP0                                                              *
*                                                                      *
*     CHARGE  :  List of atomic charges for the solute.           (I)  *
*                >> real*8 CHARGE(NATO) <<                             *
*     ALPHAL  :  Ewald sum exponential parameter.                 (I)  *
*     FPX     :  Forces for each atom of the macromolecule.      (I/O) *
*     FPY        >> real*8 FPX(NATO), FPY(NATO), FPZ(NATO) <<          *
*     FPZ                                                              *
*                                                                      *
*     UCOUL   :  Coulombic energy.                                 (O) *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER m6,sp(*)
      INTEGER nato,ingrpp,type(m6,*),ingrp(2,*),ss_index(*),nbtype(*)
      REAL*8  xp0(*),yp0(*),zp0(*),fpx(*),fpy(*),fpz(*),flj_x(*),flj_y(
     &     *),flj_z(*),alphal,charge(*),ucoul_slt,ucoul_slv,uconf_slt
     &     ,uconf_slv,ecc12(*),ecc6(*)
      REAL*8 plrzbij(*),Utotal,Edx(*),Edy(*),Edz(*),dipole(3,*),U_Thole
      REAL*8 cgi,cgj,cgij,dphij_dx,dphij_dy,dphij_dz,epol,echarge,term
      LOGICAL do_LJ

*------------------- VARIABLES IN COMMON -------------------------------

      INCLUDE 'unit.h'

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,i1,j,lij,typei,nbti,nbtj
      REAL*8 xpi,ypi,zpi,xc,yc,zc,rsq,rsp,rsqi,qforce,rspi
      REAL*8 xpj,ypj,zpj,chrgei,chrgej,aux1
      REAL*8 a1,a2,a3,a4,a5,qp,qt,expcst,erfcst,ucoul(2),uconf(2)
      REAL*8 rspqi,alphar,twrtpi,fac,switch,d_switch,r6,r12,ssvir,conf
     &     ,udirect
      DATA a1,a2,a3/0.2548296d0,-0.28449674d0,1.4214137d0/
      DATA a4,a5/-1.453152d0,1.0614054d0/
      DATA qp/0.3275911d0/
      real*8 arg,fact,cexp,B0,B1,B2,B3,Gij0,Gij11,Gij1,Gij2,B1_0,B2_0,B4
      REAL*8 BB1,BB2,BB3,BB4,d_B2,d_B3
      real*8 mui_x,mui_y,mui_z,muj_x,muj_y,muj_z,dotir,dotjr,dotij
      real*8 termi,termj,Uddd,dphii_dx,dphii_dy,dphii_dz,ucoula,ULJ
      REAL*8 D1,D2,d_D1,d_D2
      integer k,l

      INTEGER ii,mma,iam

*==================== EXECUTABLE STATEMENTS ============================

      DO i=1,2
         ucoul(i)=0.0D+0
         uconf(i)=0.0D+0
      END DO
      ucoul_slt=0.0D0
      ucoul_slv=0.0D0
      uconf_slt=0.0D0
      uconf_slv=0.0D0
      IF(ingrpp.EQ.0) RETURN
      IF(alphal == 0) THEN
         twrtpi=0.0D0
      ELSE
         twrtpi=1.0d0/DSQRT(pi)/alphal
      END IF
      fac=2.0D0*alphal*alphal
      
*=======================================================================
*----- Take Ewald for the electrostatic interaction --------------------
*=======================================================================
      mma=sp(1)
      DO ii=1,mma
         iam = sp(1+ii)
         i = IABS(iam)
         i1=ingrp(1,i)
         j=ingrp(2,i)

*-----------------------------------------------------------------------
*------------ No bonded interaction exists between solvent and solute --
*-----------------------------------------------------------------------
            
         typei=ss_index(i1)
         nbti = nbtype(i1)
         nbtj = nbtype(j)
         lij=type(nbti,nbtj)
         xpi=xp0(i1)
         ypi=yp0(i1)
         zpi=zp0(i1)
         xpj=xp0(j)
         ypj=yp0(j)
         zpj=zp0(j)
         cgi=charge(i1)

         mui_x=dipole(1,i1)
         mui_y=dipole(2,i1)
         mui_z=dipole(3,i1)

         xc=xpi-xpj
         yc=ypi-ypj
         zc=zpi-zpj
         rsq=xc**2+yc**2+zc**2
         ucoula=0.0D0
         Uddd=0.0D0
         ULJ=0.0D0
         INCLUDE 'dipole_interaction.f'
         IF(do_LJ) THEN
            rsqi=1.0D0/rsq
            r6=rsqi*rsqi*rsqi
            r12=r6*r6
            ssvir=12.0d0*ecc12(lij)*r12-6.0d0*ecc6(lij
     &           )*r6
            qforce=ssvir*rsqi
            
            conf=ecc12(lij)*r12-ecc6(lij)*r6
            
            flj_x(i1)=flj_x(i1)+qforce*xc
            flj_y(i1)=flj_y(i1)+qforce*yc
            flj_z(i1)=flj_z(i1)+qforce*zc
            flj_x(j )=flj_x(j )-qforce*xc
            flj_y(j )=flj_y(j )-qforce*yc
            flj_z(j )=flj_z(j )-qforce*zc
            uconf(typei)=uconf(typei)+ conf
         END IF
         ucoul(typei)=ucoul(typei)+ucoula
         U_Thole=U_Thole+Uddd
      END DO

      ucoul_slt=ucoul(1)
      ucoul_slv=ucoul(2)
      uconf_slt=uconf(1)
      uconf_slv=uconf(2)
      Utotal=Utotal+ucoul_slt+ucoul_slv

*================= END OF EXECUTABLE STATEMENTS ========================
      
      RETURN
      END
