      SUBROUTINE mts_forpp2(ss_index,xp0,yp0,zp0,xpg,ypg,zpg,charge
     &     ,nbtype,type,ma,nato,atomg,xpcm,ypcm,zpcm,groupp,atomp,co
     &     ,alphal,mapnl,ngrp,grppt,ucoul_slt,ucoul_slv,ucoul_ss
     &     ,uconf_slt,uconf_slv,uconf_ss,U_thole,ULJ,fpx,fpy,fpz,do_LJ
     &     ,fLJ_x,fLJ_y,fLJ_z,ecc12,ecc6,Rcut_EL,Rcut_LJ
     &     ,npp_loc,ncount1,worka,nstart,nend ,iret,errmsg,plrzbij,Edx
     &     ,Edy,Edz,dipole,skip_direct)

*****MultipleTimeScale Version*****P.Procacci-CECAM*********************
*                                                                      *
*     Compute the contribution from non-bonded interaction to          *
*     the forces and energies of the macromolecule. It does not        *
*     include code for hydrogen bonded interactions. Interactions      *
*     are switched with a group-group 3-spline function. Works both    *
*     for Ewald and NoEwald.                                           *
*                                                                      *
*     XP0     :  Coordinates of the macromolecule.                (I)  *
*     YP0        >> real*8 XP0(NATO), YP0(NATO), ZP0(NATO) <<          *
*     ZP0                                                              *
*                                                                      *
*     CHARGE  :  List of atomic charges for the macromolecule.    (I)  *
*                >> real*8 CHARGE(NATO) <<                             *
*     NBTYPE  :  List of atomic types for the macromolecule.      (I)  *
*                >> real*8 NBTYPE(NATO) <<                             *
*     NATO    :  Number of atoms forming the macromolecule.       (I)  *
*     CO      :  Transformation matrix from box coordinates       (I)  *
*                to orthogonal frame.                                  *
*                >> real*8 CO(3,3) <<                                  *
*     ALPHAL  :  Ewald sum exponential parameter.                 (I)  *
*     MAPNL   :  Integer 1-2 and 1-3 list.                        (I)  *
*                >> integer*2 MAPNL(*) <<                              *
*     UCOUL   :  Coulombic energy.                               (I/O) *
*                                                                      *
*====================== WORK ARRAYS ===================================*
*                                                                      *
*     XMAP0   :  >> real*8  XMAP1(group)<<       Non-switched          *
*     YMAP0   :  >> real*8  YMAP1(group)<<       pbc vectors           *
*     ZMAP0   :  >> real*8  ZMAP1(group)<<                             *
*                                                                      *
*     EXTERNAL 1                                                    *
*                                                                      *
************************************************************************
C======================= DECLARATIONS ==================================
      USE Module_Neighbors, ONLY: neigh=>neigha
      IMPLICIT none
C----------------------- ARGUMENTS -------------------------------------
      INTEGER nato,ngrp,ma,nstart,nend
      INTEGER nbtype(*),type(ma,*),grppt(2,*),ss_index(*),groupp(*)
     &     ,atomp(*),atomg(*),worka(*),iret
      CHARACTER*80 errmsg
      INTEGER mapnl(*),ncount1,npp_loc
      REAL*8  xp0(*),yp0(*),zp0(*),co(3,3),charge(*)
     &     ,xpg(*),ypg(*),zpg(*),xpcm(*),ypcm(*),zpcm(*),fpx(*),fpy(*)
     &     ,fpz(*),fLJ_x(*),fLJ_y(*),fLJ_z(*),ecc12(*),ecc6(*)
      REAL*8  alphal,Rcut_LJ,Rcut_EL, U_thole,ULJ,ucoul_slt,ucoul_slv
     &     ,ucoul_ss ,uconf_slt,uconf_slv,uconf_ss
      REAL*8 plrzbij(*),Utotal,Edx(*),Edy(*),Edz(*),dipole(3,*)
      LOGICAL do_LJ,skip_direct
C------------------- VARIABLES IN COMMON -------------------------------
      INCLUDE 'unit.h'
C-------------------- LOCAL VARIABLES ----------------------------------
      INTEGER i,j,jj,li,lj,lij,n,m,mm,noff,nbti,mapa,mapb,la,na,map,
     x        j1,i1,mp0,typei,typej,typeij,p1,p2,k
      REAL*8 xpi,ypi,zpi,xc,yc,zc,rsq,rsp,rsqi,qforce,
     x       xg,yg,zg,xpgi,ypgi,zpgi,xgg,ygg,zgg,drj
      REAL*8 cgi,cgj,cgij,dphij_dx,dphij_dy,dphij_dz,epol,echarge,term
     &     ,Ucharge
      REAL*8 a1,a2,a3,a4,a5,qp,qt,expcst,erfcst
      REAL*8 rspqi,alphar,furpar,twrtpi,fac,switch,d_switch
      REAL*8 Rcut_LJ2,Rcut_EL2,xmap0j,ymap0j,zmap0j,ucoul(3),uconf(3)
      REAL*8 aux1,zero,one,two,three,rspi,xd,yd,zd,xd1
     &     ,yd1,zd1,ucoula,twoi,threei,r2out,uconfa
      PARAMETER(zero=0.0D0,one=1.0D0,two=2.0D0,      
     x          three=3.0D0,twoi=0.5D0,threei=1.0D0/3.0D0)

      real*8 arg,fact,cexp
      real*8 B0,B1,B2,B3,Gij0,Gij11,Gij2,B1_0,B2_0,Uddd,Gij1,B4
      REAL*8 BB1,BB2,BB3,BB4,d_B2,d_B3
      real*8 mui_x,mui_y,mui_z,muj_x,muj_y,muj_z,dotir,dotjr,dotij
      real*8 termi,termj,dphii_dx,dphii_dy,dphii_dz,r6,r12,ssvir,conf
      integer l,p_nn
  
C------------------ DEFINITION OF A SCRATCH COMMON ---------------------
      INTEGER, ALLOCATABLE :: index0(:)
      LOGICAL, ALLOCATABLE ::  maplg(:),mapag(:),mapppa(:)
      LOGICAL :: ok
      INTEGER :: lcountb,nlocal
      REAL(8), ALLOCATABLE ::  xmap0(:),ymap0(:),zmap0(:),rgg(:),xmap1(:
     &     ),ymap1(:),zmap1(:)
      
      REAL*8 D1,D2,d_D1,d_D2,Radius
      
      DATA a1,a2,a3/0.2548296d0,-0.28449674d0,1.4214137d0/
      DATA a4,a5/-1.453152d0,1.0614054d0/
      DATA qp/0.3275911d0/
      INCLUDE 'pbc.h'

C==================== EXECUTABLE STATEMENTS ============================
      iret=0
      errmsg=' '

      nlocal=nend-nstart+1
      ALLOCATE(index0(nato),maplg(nato),mapag(nato),mapppa(nato))
      ALLOCATE(xmap0(nlocal),ymap0(nlocal),zmap0(nlocal))
      ALLOCATE(xmap1(nlocal),ymap1(nlocal),zmap1(nlocal))
      ALLOCATE(rgg(nlocal))

*=======================================================================
*---- set up the big bunch of cut-off radii
*=======================================================================

c---  radius for the outer shell
      Rcut_LJ2=Rcut_LJ**2
      Rcut_EL2=Rcut_EL**2
      r2out=MAX(Rcut_LJ2,Rcut_EL2)
      n=0
      ncount1=0
      na=0
      DO i=1,3
         ucoul(i)=0.0D0
         IF(do_LJ) THEN
            uconf(i)=0.0D0
         END IF
      END DO
      DO j=1,nato
         maplg(j)=.TRUE.
      END DO
      DO j=nstart,nend
         mapag(j)=.TRUE.
      END DO
c========mts_forpp_EWALD=================================================
c---       take EWALD SUM ONLY FOR ELECTRIC FIELD
c========================================================================
c==== start outer loop on groups

      twrtpi=1.0d0/DSQRT(pi)/alphal
      fac=2.0D0*alphal*alphal

      epol=0.0D0
      echarge=0.0D0
      p_nn=0
      DO i=nstart,nend
         p_nn=p_nn+1
         xpgi=xpg(i)
         ypgi=ypg(i)
         zpgi=zpg(i)
         m=neigh(p_nn) % no
         map=0
         mp0=0
         typei=ss_index(grppt(1,i))
         rgg(1:ngrp)=0.0D0
         DO jj=1,m
            j=neigh(p_nn) % nb(jj)
            xgg=xpgi-xpg(j)
            ygg=ypgi-ypg(j)
            zgg=zpgi-zpg(j)
            xmap0j=2.0D0*PBC(xgg)
            ymap0j=2.0D0*PBC(ygg)
            zmap0j=2.0D0*PBC(zgg)
            xgg=xgg-xmap0j
            ygg=ygg-ymap0j
            zgg=zgg-zmap0j
            xc=co(1,1)*xgg+co(1,2)*ygg+co(1,3)*zgg
            yc=            co(2,2)*ygg+co(2,3)*zgg
            zc=                        co(3,3)*zgg
            drj=xc**2+yc**2+zc**2
            
c---        map groups with "normal" (non switched) interactions
            IF (drj.LT.r2out) THEN
               mp0=mp0+1
               xmap0(mp0)=xmap0j
               ymap0(mp0)=ymap0j
               zmap0(mp0)=zmap0j
               xmap1(mp0)=xc
               ymap1(mp0)=yc
               zmap1(mp0)=zc
               index0(mp0)=j
               rgg(mp0)=drj
            ENDIF
         END DO
c---     this is the number of normal (nos switched) group interactions
         mapa=mp0
         
         IF(mapa.eq.0) THEN
            DO i1=grppt(1,i),grppt(2,i)
               noff=mapnl(1+na)+1
               na=na+noff
            END DO
            go to 1002
         END IF
         
         DO i1=grppt(1,i),grppt(2,i)
            xpi=xp0(i1)
            ypi=yp0(i1)
            zpi=zp0(i1)
            nbti=nbtype(i1)
            cgi=charge(i1)
            mui_x=dipole(1,i1)
            mui_y=dipole(2,i1)
            mui_z=dipole(3,i1)
            
            la=mapnl(1+na)
            DO j=1,la
               mm=mapnl(j+1+na)
               maplg(mm)=.FALSE.
               mapag(atomg(mm))=.FALSE.
            END DO
            
c----       compute forces
            
            DO jj=1,mapa
               j1=index0(jj)
               Radius=rgg(jj)
               typej=ss_index(grppt(1,j1))
               typeij=typei+typej-1
               xd=-xmap0(jj)
               yd=-ymap0(jj)
               zd=-zmap0(jj)
               xd1=xpi+xd
               yd1=ypi+yd
               zd1=zpi+zd
               ucoula=0.0D0
               Uddd=0.0D0
               uconfa=0.0D0
               IF(mapag(j1)) THEN
                  DO j=grppt(1,j1),grppt(2,j1)
                     lij=type(nbti,nbtype(j))
                     cgij=charge(j)*cgi
                     ok=.TRUE.
                     IF(DABS(cgij) .LT. 1.0D-10 .AND.
     &                    skip_direct) ok=.FALSE.
                     xg=xd1-xp0(j)
                     yg=yd1-yp0(j)
                     zg=zd1-zp0(j)
                     xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
                     yc=           co(2,2)*yg+co(2,3)*zg
                     zc=                      co(3,3)*zg
                     rsq=xc*xc+yc*yc+zc*zc
                     IF(Radius .LT. Rcut_EL2 .AND. ok) THEN
                        INCLUDE 'dipole_interaction.f'
                     END IF
                     IF(Radius .LT. Rcut_LJ2) THEN
                        IF(do_LJ) THEN
                           rsqi=1.0D0/rsq
                           r6=rsqi*rsqi*rsqi
                           r12=r6*r6
                           ssvir=12.0d0*ecc12(lij)*r12-6.0d0*ecc6(lij
     &                          )*r6
                           qforce=ssvir*rsqi
                           
                           conf=ecc12(lij)*r12-ecc6(lij)*r6
                           
                           ULJ=ULJ+conf
                           uconfa=uconfa+conf
                           
                           flj_x(i1)=flj_x(i1)+qforce*xc
                           flj_y(i1)=flj_y(i1)+qforce*yc
                           flj_z(i1)=flj_z(i1)+qforce*zc
                           flj_x(j )=flj_x(j )-qforce*xc
                           flj_y(j )=flj_y(j )-qforce*yc
                           flj_z(j )=flj_z(j )-qforce*zc
                        END IF
                     END IF
                  END DO
               ELSE
                  DO j=grppt(1,j1),grppt(2,j1)
                     lij=type(nbti,nbtype(j))
                     cgij=charge(j)*cgi
                     ok=.TRUE.
                     IF(DABS(cgij) .LT. 1.0D-10 .AND.
     &                    skip_direct) ok=.FALSE.
                     xg=xd1-xp0(j)
                     yg=yd1-yp0(j)
                     zg=zd1-zp0(j)
                     xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
                     yc=           co(2,2)*yg+co(2,3)*zg
                     zc=                      co(3,3)*zg
                     rsq=xc*xc+yc*yc+zc*zc
                     IF(maplg(j)) THEN
                        IF(Radius .LT. Rcut_EL2 .AND. ok) THEN
                           INCLUDE 'dipole_interaction.f'
                        END IF
                        IF(Radius .LT. Rcut_LJ2) THEN
                           IF(do_LJ) THEN
                              rsqi=1.0D0/rsq
                              r6=rsqi*rsqi*rsqi
                              r12=r6*r6
                              ssvir=12.0d0*ecc12(lij)*r12-6.0d0
     &                             *ecc6(lij)*r6
                              qforce=ssvir*rsqi
                              
                              conf=ecc12(lij)*r12-ecc6(lij)*r6
                              
                              ULJ=ULJ+conf
                              uconfa=uconfa+conf
                              
                              flj_x(i1)=flj_x(i1)+qforce*xc
                              flj_y(i1)=flj_y(i1)+qforce*yc
                              flj_z(i1)=flj_z(i1)+qforce*zc
                              flj_x(j )=flj_x(j )-qforce*xc
                              flj_y(j )=flj_y(j )-qforce*yc
                              flj_z(j )=flj_z(j )-qforce*zc
                           END IF
                        END IF
                     END IF
                  END DO
               END IF
               ucoul(typeij)=ucoul(typeij)+ucoula
               U_Thole=U_Thole+Uddd
               IF(do_LJ) THEN
                  uconf(typeij)=uconf(typeij)+uconfa
               END IF
            END DO
            
            la=mapnl(1+na)
            DO j=1,la
               mm=mapnl(j+1+na)
               maplg(mm)=.TRUE.
               mapag(atomg(mm))=.TRUE.
            END DO
            noff=mapnl(1+na)+1
            na=na+noff
            
         END DO
         
1002     CONTINUE
      END DO

      ucoul_slt=ucoul(1)
      ucoul_ss =ucoul(2)
      ucoul_slv=ucoul(3)
      IF(do_LJ) THEN
         uconf_slt=uconf(1)
         uconf_ss =uconf(2)
         uconf_slv=uconf(3)
      END IF

      DEALLOCATE(index0,maplg,mapag,mapppa)
      DEALLOCATE(xmap0,ymap0,zmap0)
      DEALLOCATE(xmap1,ymap1,zmap1)
      DEALLOCATE(rgg)

      RETURN
      END
C=======================================================================
