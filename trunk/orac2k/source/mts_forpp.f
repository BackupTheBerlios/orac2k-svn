#if defined T3E | defined _CRAY_
#define _AXPY_  saxpy
#else
#define _AXPY_  daxpy
#endif
      SUBROUTINE mts_forpp(rshell,ss_index,xp0,yp0,zp0,xpg,ypg,zpg
     &     ,charge,nbtype,type,ma,nato,atomg,xpcm,ypcm,zpcm,groupp,atomp
     &     ,co,ecc12,ecc6,ewald,alphal,erfc_spline,erfc_bin,erfc_arr
     &     ,mapnl,ngrp,grppt,uconf_slt,uconf_slv,uconf_ss,ucoul_slt
     &     ,ucoul_slv,ucoul_ss,fpx,fpy,fpz,phi,stress,rneigh,rinn,rout
     &     ,rtolinn,rtolout,gmass,iz,worka,nstart,nend,Boltz_Group,flag
     &     ,iret,errmsg)
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
*     ECC12   :  List of L-J repulsive parameters.                (I)  *
*                >> real*8 ECC12(*) <<                                 *
*     ECC6    :  List of L-J attractive parameters.               (I)  *
*                >> real*8 ECC6(*) <<                                  *
*     EWALD   :  Logical parameter. If .TRUE. the electrostatic   (I)  *
*                interaction is compute with Ewald.                    *
*                >> logical*4 EWALD <<                                 *
*     ALPHAL  :  Ewald sum exponential parameter.                 (I)  *
*     MAPNL   :  Integer 1-2 and 1-3 list.                        (I)  *
*                >> integer*2 MAPNL(*) <<                              *
*     UCONF   :  Configurational energy.                         (I/O) *
*     UCOUL   :  Coulombic energy.                               (I/O) *
*                                                                      *
*====================== WORK ARRAYS ===================================*
*                                                                      *
*     XMAP0   :  >> real*8  XMAP1(group)<<       Non-switched          *
*     YMAP0   :  >> real*8  YMAP1(group)<<       pbc vectors           *
*     ZMAP0   :  >> real*8  ZMAP1(group)<<                             *
*     XMAP1   :  >> real*8  XMAP1(group)<<       Switched pbc          *
*     YMAP1   :  >> real*8  YMAP1(group)<<       vectors               *
*     ZMAP1   :  >> real*8  ZMAP1(group)<<                             *
*     XMAP3   :  >> real*8  XMAP1(group)<<       switched              *
*     YMAP3   :  >> real*8  YMAP1(group)<<       distance vectors      *
*     ZMAP3   :  >> real*8  ZMAP1(group)<<                             *
*                                                                      *
*     EXTERNAL NONE                                                    *
*                                                                      *
************************************************************************

C======================= DECLARATIONS ==================================

      USE Module_Neighbors, ONLY: Neigh_Start=>Start, Neigh_Delete
     &     =>Delete, Neighbors, neigha, neighb, neighc
      IMPLICIT none

C----------------------- ARGUMENTS -------------------------------------

      CHARACTER*1  rshell
      INTEGER nato,ngrp,ma,nstart,nend
      INTEGER nbtype(*),type(ma,*),grppt(2,*),ss_index(*),groupp(*)
     &     ,atomp(*),atomg(*),worka(*),iret,Boltz_Group(*)
      CHARACTER*80 errmsg
      INTEGER mapnl(*),iz
      REAL*8  xp0(*),yp0(*),zp0(*),co(3,3),charge(*),ecc12(*),ecc6(*)
     &     ,xpg(*),ypg(*),zpg(*),xpcm(*),ypcm(*),zpcm(*),fpx(*),fpy(*)
     &     ,fpz(*),gmass(*),erfc_arr(4,*),phi(*),stress(3,3)
      REAL*8  alphal,uconf_slt,uconf_slv,uconf_ss,ucoul_slt,ucoul_slv
     &     ,ucoul_ss,rneigh,rinn,rout,rtolinn,rtolout,erfc_bin
      LOGICAL ewald,erfc_spline,flag

C------------------- VARIABLES IN COMMON -------------------------------

      INCLUDE 'unit.h'

C-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,jj,li,lj,lij,n,m,mm,noff,nbti,mapa,mapb,la,na,map,
     &     j1,i1,mp0,mp1,mbeg,ind,typei,typej,typeij,k,p1,btypei
     &     ,btypej,btypeij
      REAL*8 xpi,ypi,zpi,xc,yc,zc,rsq,rsp,rsqi,qforce,
     x       xg,yg,zg,xpgi,ypgi,zpgi,xgg,ygg,zgg,drj,massi,massj
      REAL*8 ssvir,r6,r12,chrgei,auxa,auxb
      REAL*8 ucon,a1,a2,a3,a4,a5,qp,qt,expcst,erfcst
      REAL*8 rspqi,alphar,furpar,twrtpi,conf,st(3,3)
      REAL*8 r2neigh,r2inn,r2out,rinn0,r2inn0,rout0,r2out0,arsout1
     &     ,arsout2,arsinn1,arsinn2,xmap0j,ymap0j,zmap0j,uconf(3)
     &     ,ucoul(3),st1,st2,st3,st4,st5,st6,st7,st8,st9,emvir,qfx,qfy
     &     ,qfz
      LOGICAL  lskip_ewald,ok
      REAL*8 xx,dxx,aux1,zero,one,two,three,rspi,derfcst,xd,yd,zd,xd1
     &     ,yd1,zd1,ucoula,uconfa,twoi,threei
      PARAMETER(zero=0.0D0,one=1.0D0,two=2.0D0,      
     x          three=3.0D0,twoi=0.5D0,threei=1.0D0/3.0D0)

C------------------ DEFINITION OF A SCRATCH COMMON ---------------------

      INTEGER nb
      INCLUDE 'parst.h'
      INTEGER lcountb,type_boltz,nlocal
      PARAMETER(nb=m1)

      INTEGER :: p_mapa,p_mapb,p_j,p_nn,ncount1
      INTEGER, DIMENSION(:), ALLOCATABLE :: p_index_jj,p_index_j
      INTEGER, DIMENSION(:), ALLOCATABLE :: index0,index1,ind_a
      REAL(8), DIMENSION(:), ALLOCATABLE :: xmap0,ymap0,zmap0,xmap3
     &     ,ymap3,zmap3,xmap1,ymap1,zmap1,cmap2,swrs,dswrs
     &     ,xmap2,ymap2,zmap2,fppx,fppy,fppz
      LOGICAL, DIMENSION(:), ALLOCATABLE :: maplg,mapag
      REAL*8 erftbdns
      TYPE(Neighbors), DIMENSION (:), POINTER :: neigh
      TYPE(Neighbors), DIMENSION (:), POINTER :: neigh1

      DATA a1,a2,a3/0.2548296d0,-0.28449674d0,1.4214137d0/
      DATA a4,a5/-1.453152d0,1.0614054d0/
      DATA qp/0.3275911d0/
      INCLUDE 'pbc.h'

C==================== EXECUTABLE STATEMENTS ============================

      iret=0
      errmsg=' '
      IF(nb.LT.nato) THEN
          errmsg=' IN FNBOND : Dimensions of the work arrays '/ /
     &        'insufficient. ABORTC '
          iret=1
          RETURN
      END IF

*=======================================================================
*---- Set up dynamic memory 
*=======================================================================


      ALLOCATE(xmap0(ngrp),ymap0(ngrp),zmap0(ngrp),xmap1(ngrp
     &     ),ymap1(ngrp),zmap1(ngrp),xmap2(ngrp),ymap2(ngrp)
     &     ,zmap2(ngrp),xmap3(ngrp),ymap3(ngrp),zmap3(ngrp)
     &     ,cmap2(ngrp),swrs(ngrp),dswrs(ngrp),index0(ngrp)
     &     ,index1(ngrp),maplg(nato),mapag(nato),fppx(nato)
     &     ,fppy(nato),fppz(nato),ind_a(ngrp))
      ALLOCATE(p_index_jj(nato),p_index_j(nato))
      CALL zero3(fppx,fppy,fppz,1,nato)

*=======================================================================
*---- set up the big bunch of cut-off radii
*=======================================================================

      IF(erfc_spline) erftbdns=1.0D0/erfc_bin

c---  radius for the inner shell
      r2inn=rinn**2

c---  radius for the outer shell
      r2out=rout**2

c---  inner shell radius + healing lenght
      rinn0 = rinn+rtolinn
      r2inn0 =rinn0**2
      if(rinn.lt.1.d-5) rtolinn=0.d0

c---  outer shell radius +healing lenght
      rout0 = rout+rtolout
      r2out0 = rout0**2

c---- tolerance radius for the neighbor list of the inner shell
c---- to be used for shorter ranged interactions
      r2neigh=(rinn0+rneigh)**2

      qt=1.0d0/(1.0d0+qp*alphal*rinn)
      expcst=DEXP(-alphal*rinn*alphal*rinn)
      erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)
     x     *qt+a1)*qt*expcst
      lskip_ewald = erfcst.lt.10.d-4
      twrtpi=2.0d0/DSQRT(pi)
      arsout1=rout0-3.0d0*rout
      arsout2=rtolout**3
      arsinn1=rinn0-3.0d0*rinn
      arsinn2=rtolinn**3
      n=0
      ncount1=0
      na=0

      st1=0.0D0
      st2=0.0D0
      st3=0.0D0
      st4=0.0D0
      st5=0.0D0
      st6=0.0D0
      st7=0.0D0
      st8=0.0D0
      st9=0.0D0

      DO i=1,3
         uconf(i)=0.0D0
         ucoul(i)=0.0D0
      END DO

      DO j=1,nato
         maplg(j)=.TRUE.
      END DO

      DO j=nstart,nend
         mapag(j)=.TRUE.
      END DO
      nlocal=nend-nstart+1
      IF(iz == 0) THEN
         CALL Neigh_Delete(neigha)
         CALL Neigh_Start(neigha,nlocal)
         neigh1=>neigha
      ELSE
         IF(rshell == 'h') THEN
            neigh=>neigha
            CALL Neigh_Delete(neighb)
            CALL Neigh_Start(neighb,nlocal)
            neigh1=>neighb
         ELSE IF(rshell == 'l') THEN
            neigh=>neighb
            CALL Neigh_Delete(neighc)
            CALL Neigh_Start(neighc,nlocal)
            neigh1=>neighc
         ELSE IF(rshell == 'm') THEN
            neigh=>neighc
         END IF
      END IF

c========mts_forpp_EWALD=================================================
c---       take EWALD SUM 
c========================================================================

      if(.not.ewald) go to 2005
c==== start outer loop on groups
      twrtpi=2.0d0/DSQRT(pi)
      mbeg=1
      p_nn=0
      DO i=nstart,nend
         
         p_nn=p_nn+1
         xpgi=xpg(i)
         ypgi=ypg(i)
         zpgi=zpg(i)
         m=ngrp
         IF(iz /= 0) m=neigh(p_nn) % no
         map=0
         mp0=0
         mp1=0
         typei=ss_index(grppt(1,i))
         btypei=Boltz_Group(i)
c------------------------------------------------------------------------
c        build up neighbor list for the next inner shell and
c        set up group-group contact maps for:
c        interactions rinn0<r<rout (non switched) [INDEX0]
c                     (rinn<r<rinn0.and.rout<r<rout0) (switched) [INDEX1]
c------------------------------------------------------------------------
         if(iz.eq.0) THEN
            mbeg=0
            IF(flag) THEN
               mbeg=i+1
            ELSE
               mbeg=1
            END IF
         END IF
         DO jj=mbeg,m
            j=jj
            IF(iz /= 0) j=neigh(p_nn) % nb(jj)
            btypej=Boltz_Group(j)
            btypeij=btypei+btypej-1
            IF(btypeij .EQ. 4) GOTO 2003
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
c---        do neighbor list for next inner shell
            if(r2inn.gt.1.0d-5.and.drj.LT.r2neigh) then
               map=map+1
               ind_a(map)=j
            endif   
c---        skip mapping if outer shell is zero
            IF(rinn0.gt.rout) go to 2003
c---        map groups with "normal" (non switched) interactions
            IF (drj.GE.r2inn0.AND.drj.LT.r2out) THEN
               mp0=mp0+1
               xmap0(mp0)=xmap0j
               ymap0(mp0)=ymap0j
               zmap0(mp0)=zmap0j
               index0(mp0)=j
c---        map groups in outer and inner "switched" shell
            ELSE
c----          if (rout < r < rout0 ) the switching function is S(x) 
               IF (drj.GE.r2out.and.drj.lt.r2out0) THEN
                  mp1=mp1+1
                  xmap1(mp1)=xmap0j
                  ymap1(mp1)=ymap0j
                  zmap1(mp1)=zmap0j
                  xmap2(mp1)=xc
                  ymap2(mp1)=yc
                  zmap2(mp1)=zc

                  index1(mp1)=j
                  rsp=DSQRT(drj)
                  auxa=(arsout1+2.0d0*rsp)/arsout2
                  auxb=rout0-rsp
                  swrs(mp1)=auxa*auxb**2
                  dswrs(mp1)=-2.0d0*auxa*auxb+2.0d0*auxb**2/arsout2
                  dswrs(mp1)=dswrs(mp1)/rsp
               ELSE IF(drj.GE.r2inn.and.drj.LT.r2inn0) then
c----          if (rinn < r < rinn0 ) the switching function is 1 - S(x) 
                  mp1=mp1+1
                  xmap1(mp1)=xmap0j
                  ymap1(mp1)=ymap0j
                  zmap1(mp1)=zmap0j
                  xmap2(mp1)=xc
                  ymap2(mp1)=yc
                  zmap2(mp1)=zc

                  index1(mp1)=j
                  rsp=DSQRT(drj)
                  auxa=(arsinn1+2.0d0*rsp)/arsinn2
                  auxb=rinn0-rsp
                  swrs(mp1)=1.d0-auxa*auxb**2
                  dswrs(mp1)=2.0d0*auxa*auxb-2.0d0*auxb**2/arsinn2
                  dswrs(mp1)=dswrs(mp1)/rsp
               ENDIF   
            ENDIF
2003        continue
         END DO

         IF(rinn.gt.1.d-5) then
            neigh1(p_nn) % no = map
            ALLOCATE(neigh1(p_nn) % nb(map))
            neigh1(p_nn) % nb = ind_a(1:map)
            ncount1=ncount1+map+1
         endif

c---     this is the number of normal (nos switched) group interactions
         mapa=mp0

c---     this is the number of switched [S or (1-S)] group interactions
         mapb=mp1
         
c---     set auxiliary array cmap2 for swich to zero
         DO jj=1,mapb
            cmap2(jj)=0.d0
         END DO
         IF(rinn0 .lt. rout) THEN
            IF(mapa.eq.0.and.mapb.eq.0) THEN
               DO i1=grppt(1,i),grppt(2,i)
                  noff=mapnl(1+na)+1
                  na=na+noff
               END DO
               go to 1002
            END IF
         ELSE
            go to 1002
         END IF
         DO i1=grppt(1,i),grppt(2,i)
            xpi=xp0(i1)
            ypi=yp0(i1)
            zpi=zp0(i1)
            nbti=nbtype(i1)
            chrgei=charge(i1)
            la=mapnl(1+na)
            DO j=1,la
               mm=mapnl(j+1+na)
               maplg(mm)=.FALSE.
               mapag(atomg(mm))=.FALSE.
            END DO

c----       compute forces
            p_mapa=0
            DO jj=1,mapa
               j1=index0(jj)
               IF(mapag(j1)) THEN
                  DO j=grppt(1,j1),grppt(2,j1)
                     p_mapa=p_mapa+1
                     p_index_jj(p_mapa)=jj
                     p_index_j(p_mapa)=j
                  END DO
               ELSE
                  DO j=grppt(1,j1),grppt(2,j1)
                     IF(maplg(j)) THEN
                        p_mapa=p_mapa+1
                        p_index_jj(p_mapa)=jj
                        p_index_j(p_mapa)=j
                     END IF
                  END DO
               END IF
            END DO
            IF(.NOT. lskip_ewald) THEN
               IF(erfc_spline) THEN
                  DO p_j=1,p_mapa
                     jj=p_index_jj(p_j)
                     j=p_index_j(p_j)
                     j1=index0(jj)
                     typej=ss_index(grppt(1,j1))
                     typeij=typei+typej-1
                     xd=-xmap0(jj)
                     yd=-ymap0(jj)
                     zd=-zmap0(jj)
                     xd1=xpi+xd
                     yd1=ypi+yd
                     zd1=zpi+zd
                     uconfa=0.0D0
                     ucoula=0.0D0
#include "force_ele_spline.f"
                     ucoul(typeij)=ucoul(typeij)+ucoula
                     uconf(typeij)=uconf(typeij)+uconfa
                  END DO
               ELSE
                  DO p_j=1,p_mapa
                     jj=p_index_jj(p_j)
                     j=p_index_j(p_j)
                     j1=index0(jj)
                     typej=ss_index(grppt(1,j1))
                     typeij=typei+typej-1
                     xd=-xmap0(jj)
                     yd=-ymap0(jj)
                     zd=-zmap0(jj)
                     xd1=xpi+xd
                     yd1=ypi+yd
                     zd1=zpi+zd
                     uconfa=0.0D0
                     ucoula=0.0D0
#include "force_ele_nospline.f"
                     ucoul(typeij)=ucoul(typeij)+ucoula
                     uconf(typeij)=uconf(typeij)+uconfa
                  END DO
               END IF
            ELSE
               DO p_j=1,p_mapa
                  jj=p_index_jj(p_j)
                  j=p_index_j(p_j)
                  j1=index0(jj)
                  typej=ss_index(grppt(1,j1))
                  typeij=typei+typej-1
                  xd=-xmap0(jj)
                  yd=-ymap0(jj)
                  zd=-zmap0(jj)
                  xd1=xpi+xd
                  yd1=ypi+yd
                  zd1=zpi+zd
                  uconfa=0.0D0
#include "force_noele.f"
                  uconf(typeij)=uconf(typeij)+uconfa
               END DO
            END IF
               
c----       compute switched forces:
c-----      calculate S(r)dV/dr and map sum_{k1k2} V_{k1,k2} onto cmap2

            p_mapb=0
            DO jj=1,mapb
               j1=index1(jj)
               IF(mapag(j1)) THEN
                  DO j=grppt(1,j1),grppt(2,j1)
                     p_mapb=p_mapb+1
                     p_index_jj(p_mapb)=jj
                     p_index_j(p_mapb)=j
                  END DO
               ELSE
                  DO j=grppt(1,j1),grppt(2,j1)
                     IF(maplg(j)) THEN
                        p_mapb=p_mapb+1
                        p_index_jj(p_mapb)=jj
                        p_index_j(p_mapb)=j
                     END IF
                  END DO
               END IF
            END DO
            IF(.NOT. lskip_ewald) THEN
               IF(erfc_spline) THEN
                  DO p_j=1,p_mapb
                     jj=p_index_jj(p_j)
                     j=p_index_j(p_j)
                     j1=index1(jj)
                     typej=ss_index(grppt(1,j1))
                     typeij=typei+typej-1
#include "force_ele_spline_sw.f"
                  END DO
               ELSE
                  DO p_j=1,p_mapb
                     jj=p_index_jj(p_j)
                     j=p_index_j(p_j)
                     j1=index1(jj)
                     typej=ss_index(grppt(1,j1))
                     typeij=typei+typej-1
#include "force_ele_nospline_sw.f"
                  END DO
               END IF
            ELSE
               DO p_j=1,p_mapb
                  jj=p_index_jj(p_j)
                  j=p_index_j(p_j)
                  j1=index1(jj)
                  typej=ss_index(grppt(1,j1))
                  typeij=typei+typej-1
#include "force_noele_sw.f"
               END DO
            END IF
            la=mapnl(1+na)
            DO j=1,la
               mm=mapnl(j+1+na)
               maplg(mm)=.TRUE.
               mapag(atomg(mm))=.TRUE.
            END DO
            noff=mapnl(1+na)+1
            na=na+noff

         END DO
         
c===     add the S*V term to the atomic forces

         DO jj=1,mapb
            xmap3(jj)=-dswrs(jj)*cmap2(jj)*xmap2(jj)
            ymap3(jj)=-dswrs(jj)*cmap2(jj)*ymap2(jj)
            zmap3(jj)=-dswrs(jj)*cmap2(jj)*zmap2(jj)
         END DO

         if(mapb.ne.0) then
            DO i1=grppt(1,i),grppt(2,i)
               massi=gmass(i1)
               xpi=xp0(i1)
               ypi=yp0(i1)
               zpi=zp0(i1)
               do jj=1,mapb
                  fppx(i1)=fppx(i1)+massi*xmap3(jj)
                  fppy(i1)=fppy(i1)+massi*ymap3(jj)
                  fppz(i1)=fppz(i1)+massi*zmap3(jj)
                  st1=st1+massi*xmap3(jj)*xpi
                  st2=st2+massi*xmap3(jj)*ypi
                  st3=st3+massi*xmap3(jj)*zpi
                  st4=st4+massi*ymap3(jj)*xpi
                  st5=st5+massi*ymap3(jj)*ypi
                  st6=st6+massi*ymap3(jj)*zpi
                  st7=st7+massi*zmap3(jj)*xpi
                  st8=st8+massi*zmap3(jj)*ypi
                  st9=st9+massi*zmap3(jj)*zpi
               end do
            end do
            do jj=1,mapb
               j=index1(jj)
               xg=xmap1(jj)
               yg=ymap1(jj)
               zg=zmap1(jj)
               do j1=grppt(1,j),grppt(2,j)
                  massj=gmass(j1)
                  fppx(j1)=fppx(j1)-massj*xmap3(jj)
                  fppy(j1)=fppy(j1)-massj*ymap3(jj)
                  fppz(j1)=fppz(j1)-massj*zmap3(jj)
                  xpi=-(xp0(j1)+xg)*massj
                  ypi=-(yp0(j1)+yg)*massj
                  zpi=-(zp0(j1)+zg)*massj
                  st1=st1+xmap3(jj)*xpi
                  st2=st2+xmap3(jj)*ypi
                  st3=st3+xmap3(jj)*zpi
                  st4=st4+ymap3(jj)*xpi
                  st5=st5+ymap3(jj)*ypi
                  st6=st6+ymap3(jj)*zpi
                  st7=st7+zmap3(jj)*xpi
                  st8=st8+zmap3(jj)*ypi
                  st9=st9+zmap3(jj)*zpi
               end do   
            end do
         end if
1002     CONTINUE
      END DO
      ucoul_slt=ucoul(1)
      ucoul_ss =ucoul(2)
      ucoul_slv=ucoul(3)

      uconf_slt=uconf(1)
      uconf_ss =uconf(2)
      uconf_slv=uconf(3)

      IF(DABS(rout-rinn) .GT. 10D-5) THEN
         CALL _AXPY_(nato,1.0D0,fppx,1,fpx,1)
         CALL _AXPY_(nato,1.0D0,fppy,1,fpy,1)
         CALL _AXPY_(nato,1.0D0,fppz,1,fpz,1)
         DO i=1,nato
            p1=atomp(i)
            xg=xp0(i)-xpcm(p1)
            yg=yp0(i)-ypcm(p1)
            zg=zp0(i)-zpcm(p1)
            xc=-fppx(i)
            yc=-fppy(i)
            zc=-fppz(i)
            st1=st1+xc*xg
            st2=st2+xc*yg
            st3=st3+xc*zg
            st4=st4+yc*xg
            st5=st5+yc*yg
            st6=st6+yc*zg
            st7=st7+zc*xg
            st8=st8+zc*yg
            st9=st9+zc*zg
         END DO
         st(1,1)=st1
         st(1,2)=st2
         st(1,3)=st3
         st(2,1)=st4
         st(2,2)=st5
         st(2,3)=st6
         st(3,1)=st7
         st(3,2)=st8
         st(3,3)=st9
         
         DO i=1,3
            DO j=1,3
               DO k=1,3
                  stress(i,j) = stress(i,j)+st(i,k)*co(j,k)
               END DO
            END DO
         END DO
      END IF
      DEALLOCATE(xmap0,ymap0,zmap0,xmap1,ymap1,zmap1,xmap2,ymap2,zmap2
     &     ,xmap3,ymap3,zmap3,cmap2,swrs,dswrs,index0,index1,maplg,mapag
     &     ,fppx,fppy,fppz,ind_a,p_index_jj,p_index_j)
      IF(iz .EQ. 0) THEN
         p_nn=0
         DO i=nstart,nend
            p_nn=p_nn+1
            worka(i) = neigh1(p_nn) % no
         END DO
      END IF
      NULLIFY(neigh,neigh1)

      RETURN

c========END OF EWALD SUM FORCE ROUTINE================================

2005  CONTINUE

c=======mts_forpp NOEWALD================================================
c---       do NOT take EWALD SUM 
c========================================================================

      mbeg=1
      p_nn=0
      DO i=nstart,nend
         p_nn=p_nn+1
         xpgi=xpg(i)
         ypgi=ypg(i)
         zpgi=zpg(i)
         m=ngrp
         IF(iz /= 0) m=neigh(p_nn) % no
         map=0
         mp0=0
         mp1=0
         typei=ss_index(grppt(1,i))

c------------------------------------------------------------------------
c        build up neighbor list for the next inner shell and
c        set up group-group contact maps for:
c        interactions rinn0<r<rout (non switched)
c                     rinn<r<rinn0 and  rout < r < rout0 (switched)
c------------------------------------------------------------------------
         if(iz.eq.0) THEN
            mbeg=0
            IF(flag) THEN
               mbeg=i+1
            ELSE
               mbeg=1
            END IF
         END IF
         DO jj=mbeg,m
            j=jj
            IF(iz /= 0) j=neigh(p_nn) % nb(jj)
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
c---        do neighbor list for next inner shell
            if(r2inn.gt.1.0d-5.and.drj.LT.r2neigh) then
               map=map+1
               ind_a(map)=j
            endif   
c---        skip mapping if outer shell is zero
            IF(rinn0.gt.rout) go to 3003
c---        map groups with "normal" (non switched) interactions
            IF (drj.GT.r2inn0.AND.drj.LT.r2out) THEN
               mp0=mp0+1
               xmap0(mp0)=xmap0j
               ymap0(mp0)=ymap0j
               zmap0(mp0)=zmap0j
               index0(mp0)=j
c---        map groups in outer and inner "switched" shell
            ELSE
c----          if (rout < r < rout0 ) the switching function is S(x) 
               IF (drj.GT.r2out.and.drj.lt.r2out0) THEN
                  mp1=mp1+1
                  xmap1(mp1)=xmap0j
                  ymap1(mp1)=ymap0j
                  zmap1(mp1)=zmap0j
                  xmap2(mp1)=xc
                  ymap2(mp1)=yc
                  zmap2(mp1)=zc

                  index1(mp1)=j
                  rsp=DSQRT(drj)
                  auxa=(arsout1+2.0d0*rsp)/arsout2
                  auxb=rout0-rsp
                  swrs(mp1)=auxa*auxb**2
                  dswrs(mp1)=-2.0d0*auxa*auxb+2.0d0*auxb**2/arsout2
                  dswrs(mp1)=dswrs(mp1)/rsp
               ELSE IF(drj.gt.r2inn.and.drj.lt.r2inn0) then
c----          if (rinn < r < rinn0 ) the switching function is 1 - S(x) 
                  mp1=mp1+1
                  xmap1(mp1)=xmap0j
                  ymap1(mp1)=ymap0j
                  zmap1(mp1)=zmap0j
                  xmap2(mp1)=xc
                  ymap2(mp1)=yc
                  zmap2(mp1)=zc

                  index1(mp1)=j
                  rsp=DSQRT(drj)
                  auxa=(arsinn1+2.0d0*rsp)/arsinn2
                  auxb=rinn0-rsp
                  swrs(mp1)=1.d0-auxa*auxb**2
                  dswrs(mp1)=2.0d0*auxa*auxb-2.0d0*auxb**2/arsinn2
                  dswrs(mp1)=dswrs(mp1)/rsp
               ENDIF   
            ENDIF
3003        continue
         END DO

         IF(rinn.gt.1.d-5) then 
            neigh1(p_nn) % no = map
            ALLOCATE(neigh1(p_nn) % nb(map))
            neigh1(p_nn) % nb = ind_a(1:map)
            ncount1=ncount1+map+1
         endif

c---     this is the number of normal (nos switched) group interactions
         mapa=mp0

c---     this is the number of switched [S or (1-S)] group interactions
         mapb=mp1
         
c------------------------------------------------------------------------
c        start looping on atomic contacs
c------------------------------------------------------------------------

         IF(rinn0 .lt. rout) THEN
            IF(mapa.eq.0.and.mapb.eq.0) THEN
               DO i1=grppt(1,i),grppt(2,i)
                  noff=mapnl(1+na)+1
                  na=na+noff
               END DO
               go to 2002
            END IF
         ELSE
            go to 2002
         END IF

c---     set auxiliary array cmap2 for swich to zero
         DO jj=1,mapb
            cmap2(jj)=0.d0
         END DO
         DO i1=grppt(1,i),grppt(2,i)
            xpi=xp0(i1)
            ypi=yp0(i1)
            zpi=zp0(i1)
            nbti=nbtype(i1)
            chrgei=charge(i1)
            la=mapnl(1+na)
            DO j=1,la
               mm=mapnl(j+1+na)
               maplg(mm)=.FALSE.
            END DO
c----       compute NON switched forces
            DO jj=1,mapa
               j1=index0(jj)
               typej=ss_index(grppt(1,j1))
               typeij=typei+typej-1
               DO j=grppt(1,j1),grppt(2,j1)
                  IF(maplg(j)) THEN
                     lij=type(nbti,nbtype(j))
                     xg=xpi-xp0(j)-xmap0(jj)
                     yg=ypi-yp0(j)-ymap0(jj)
                     zg=zpi-zp0(j)-zmap0(jj)
                     xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
                     yc=           co(2,2)*yg+co(2,3)*zg
                     zc=                      co(3,3)*zg
                     rsq=xc*xc+yc*yc+zc*zc
                     rsp=DSQRT(rsq)
                     rsqi=1.0d0/rsq
                     r6=rsqi*rsqi*rsqi
                     r12=r6*r6
                     ssvir=12.0d0*ecc12(lij)*r12-6.0d0*ecc6(lij)*r6
                     qforce=ssvir*rsqi
                     uconf(typeij)=uconf(typeij)+ecc12(lij)*r12-ecc6(lij
     &                    )*r6
                     ssvir=chrgei*charge(j)/rsp
                     ucoul(typeij)=ucoul(typeij)+ssvir
                     qforce=qforce+ssvir*rsqi
                     emvir=qforce

                     fppx(i1)=fppx(i1)+qforce*xc
                     fppy(i1)=fppy(i1)+qforce*yc
                     fppz(i1)=fppz(i1)+qforce*zc
                     fppx(j)=fppx(j)-qforce*xc
                     fppy(j)=fppy(j)-qforce*yc
                     fppz(j)=fppz(j)-qforce*zc
                     qfx=emvir*xc
                     qfy=emvir*yc
                     qfz=emvir*zc
                     st1 = st1+qfx*xg
                     st2 = st2+qfx*yg
                     st3 = st3+qfx*zg
                     st4 = st4+qfy*xg
                     st5 = st5+qfy*yg
                     st6 = st6+qfy*zg
                     st7 = st7+qfz*xg
                     st8 = st8+qfz*yg
                     st9 = st9+qfz*zg
                  END IF
               END DO
            END DO

c----       compute switched forces
            DO jj=1,mapb
               j1=index1(jj)
               typej=ss_index(grppt(1,j1))
               typeij=typei+typej-1
               DO j=grppt(1,j1),grppt(2,j1)
                  IF(maplg(j)) THEN
                     lij=type(nbti,nbtype(j))
                     xg=xpi-xp0(j)-xmap1(jj)
                     yg=ypi-yp0(j)-ymap1(jj)
                     zg=zpi-zp0(j)-zmap1(jj)
                     xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
                     yc=           co(2,2)*yg+co(2,3)*zg
                     zc=                      co(3,3)*zg
                     rsq=xc*xc+yc*yc+zc*zc
                     rsp=DSQRT(rsq)
                     rsqi=1.0d0/rsq
                     r6=rsqi*rsqi*rsqi
                     r12=r6*r6
                     ssvir=12.0d0*ecc12(lij)*r12-6.0d0*ecc6(lij)*r6
                     ucon=ecc12(lij)*r12-ecc6(lij)*r6
                     cmap2(jj)=cmap2(jj)+ucon
                     uconf(typeij)=uconf(typeij)+ucon*swrs(jj)
                     qforce=ssvir*rsqi*swrs(jj)
                     ssvir=chrgei*charge(j)/rsp
                     cmap2(jj)=cmap2(jj)+ssvir
                     ucoul(typeij)=ucoul(typeij)+ssvir*swrs(jj)
                     qforce=qforce+ssvir*rsqi*swrs(jj)
                     emvir=qforce

                     fppx(i1)=fppx(i1)+qforce*xc
                     fppy(i1)=fppy(i1)+qforce*yc
                     fppz(i1)=fppz(i1)+qforce*zc
                     fppx(j)=fppx(j)-qforce*xc
                     fppy(j)=fppy(j)-qforce*yc
                     fppz(j)=fppz(j)-qforce*zc
                     st1 = st1+emvir*xc*xg
                     st2 = st2+emvir*xc*yg
                     st3 = st3+emvir*xc*zg
                     st4 = st4+emvir*yc*xg
                     st5 = st5+emvir*yc*yg
                     st6 = st6+emvir*yc*zg
                     st7 = st7+emvir*zc*xg
                     st8 = st8+emvir*zc*yg
                     st9 = st9+emvir*zc*zg
                  END IF
               END DO
            END DO
            la=mapnl(1+na)
            DO j=1,la
               mm=mapnl(j+1+na)
               maplg(mm)=.TRUE.
            END DO
            noff=mapnl(1+na)+1
            na=na+noff
         END DO
         
c===     add the S*V term to the atomic forces 

         DO jj=1,mapb
            xmap3(jj)=-dswrs(jj)*cmap2(jj)*xmap2(jj)
            ymap3(jj)=-dswrs(jj)*cmap2(jj)*ymap2(jj)
            zmap3(jj)=-dswrs(jj)*cmap2(jj)*zmap2(jj)
         END DO

         if(mapb.ne.0) then
            DO i1=grppt(1,i),grppt(2,i)
               massi=gmass(i1)
               xpi=xp0(i1)
               ypi=yp0(i1)
               zpi=zp0(i1)
               do jj=1,mapb
                  fppx(i1)=fppx(i1)+massi*xmap3(jj)
                  fppy(i1)=fppy(i1)+massi*ymap3(jj)
                  fppz(i1)=fppz(i1)+massi*zmap3(jj)
                  st1=st1+massi*xmap3(jj)*xpi
                  st2=st2+massi*xmap3(jj)*ypi
                  st3=st3+massi*xmap3(jj)*zpi
                  st4=st4+massi*ymap3(jj)*xpi
                  st5=st5+massi*ymap3(jj)*ypi
                  st6=st6+massi*ymap3(jj)*zpi
                  st7=st7+massi*zmap3(jj)*xpi
                  st8=st8+massi*zmap3(jj)*ypi
                  st9=st9+massi*zmap3(jj)*zpi
               end do
            end do
            do jj=1,mapb
               j=index1(jj)
               xg=xmap1(jj)
               yg=ymap1(jj)
               zg=zmap1(jj)
               do j1=grppt(1,j),grppt(2,j)
                  massj=gmass(j1)
                  fppx(j1)=fppx(j1)-massj*xmap3(jj)
                  fppy(j1)=fppy(j1)-massj*ymap3(jj)
                  fppz(j1)=fppz(j1)-massj*zmap3(jj)
                  xpi=-(xp0(j1)+xg)*massj
                  ypi=-(yp0(j1)+yg)*massj
                  zpi=-(zp0(j1)+zg)*massj
                  st1=st1+xmap3(jj)*xpi
                  st2=st2+xmap3(jj)*ypi
                  st3=st3+xmap3(jj)*zpi
                  st4=st4+ymap3(jj)*xpi
                  st5=st5+ymap3(jj)*ypi
                  st6=st6+ymap3(jj)*zpi
                  st7=st7+zmap3(jj)*xpi
                  st8=st8+zmap3(jj)*ypi
                  st9=st9+zmap3(jj)*zpi
               end do   
            end do
         end if
2002     CONTINUE
      END DO
      ucoul_slt=ucoul(1)
      ucoul_ss =ucoul(2)
      ucoul_slv=ucoul(3)

      uconf_slt=uconf(1)
      uconf_ss =uconf(2)
      uconf_slv=uconf(3)

      IF(DABS(rout-rinn) .GT. 10D-5) THEN
         CALL _AXPY_(nato,1.0D0,fppx,1,fpx,1)
         CALL _AXPY_(nato,1.0D0,fppy,1,fpy,1)
         CALL _AXPY_(nato,1.0D0,fppz,1,fpz,1)
         DO i=1,nato
            p1=atomp(i)
            xg=xp0(i)-xpcm(p1)
            yg=yp0(i)-ypcm(p1)
            zg=zp0(i)-zpcm(p1)
            xc=-fppx(i)
            yc=-fppy(i)
            zc=-fppz(i)
            st1=st1+xc*xg
            st2=st2+xc*yg
            st3=st3+xc*zg
            st4=st4+yc*xg
            st5=st5+yc*yg
            st6=st6+yc*zg
            st7=st7+zc*xg
            st8=st8+zc*yg
            st9=st9+zc*zg
         END DO
         st(1,1)=st1
         st(1,2)=st2
         st(1,3)=st3
         st(2,1)=st4
         st(2,2)=st5
         st(2,3)=st6
         st(3,1)=st7
         st(3,2)=st8
         st(3,3)=st9
         
         DO i=1,3
            DO j=1,3
               DO k=1,3
                  stress(i,j) = stress(i,j)+st(i,k)*co(j,k)
               END DO
            END DO
         END DO
      END IF
      DEALLOCATE(xmap0,ymap0,zmap0,xmap1,ymap1,zmap1,xmap2,ymap2,zmap2
     &     ,xmap3,ymap3,zmap3,cmap2,swrs,dswrs,index0,index1,maplg,mapag
     &     ,fppx,fppy,fppz,ind_a,p_index_jj,p_index_j)
      IF(iz .EQ. 0) THEN
         p_nn=0
         DO i=nstart,nend
            p_nn=p_nn+1
            worka(i) = neigh1(p_nn) % no
         END DO
      END IF
      NULLIFY(neigh,neigh1)
      RETURN
      END
