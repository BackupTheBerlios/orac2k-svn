      SUBROUTINE comp_cavities_vor(nstart,nend,ss_index,nx,ny,nz
     &     ,bin_size,dr_histo,rmax_histo,mqq,jmax,xp0,yp0,zp0,xpg,ypg
     &     ,zpg,ntap,atomg,nbtype,res,mback,nbone,nbun,sigma,co,iret
     &     ,errmsg)

************************************************************************
*   Time-stamp: <2006-04-05 14:27:29 marchi>                             *
*                                                                      *
*   Compute cavities in a system composed of solute and solvent.       *
*   The differential cavity distribution -dp(R)/dR is binned           *
*   in an histogram contained in cavity_h. Only the histogram          *
*   for each solute residue are computed separately. The solvent       *
*   distribution is accumulated in only one histogram.                 *
*                                                                      *
*   The algorithm assigns first the points of the cavity grid          *
*   to Voronoi cells. After this operation the distance from the       *
*   grid points and the neighboring atoms is computed.                 *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Apr 11 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      USE VORONOI_Mod, ONLY: nnlpp_vor,pnnlpp_vor
      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nstart,nend,ntap,nx,ny,nz,mqq,nbun,jmax,nbone,atomg(*)
     &     ,nbtype(*),ss_index(*),res(*),mback(*)
      REAL*8  bin_size,dr_histo,rmax_histo,co(3,3),xp0(*),yp0(*),zp0(*)
     &     ,xpg(*),ypg(*),zpg(*),sigma(*)
      INTEGER iret
      CHARACTER*80 errmsg

*----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'parst.h'
      INCLUDE 'analysis.h'

      LOGICAL ok(maxcav)
      REAL*8  xmap0(maxcav),ymap0(maxcav),zmap0(maxcav),xmap1(maxcav)
     &     ,ymap1(maxcav),zmap1(maxcav),rmap(pnnlpp_vor)
      INTEGER imap(maxcav),jmap(maxcav),kmap(maxcav),index(maxcav)
     &     ,indexa(maxcav)
      COMMON /rag1/ xmap0,ymap0,zmap0,xmap1,ymap1,zmap1,rmap,imap,jmap
     &     ,kmap,index,ok

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8  binx,biny,binz,y,fact,xg,yg,zg,xc,yc,zc,xd,yd,zd,xpi,ypi
     &     ,zpi,xpgi,ypgi,zpgi,dx,dy,dz,xgg,ygg,zgg,xmap,ymap,zmap,rs
     &     ,dot,rsp,sigi,aux
      INTEGER idx,idxm,idy,idym,idz,idzm,count,map,mapp,ia,ib,ic,i,j,ja
     &     ,k,ka,m,typei,typej,tot_map,nmax,ibin,o
      INTEGER GET_IND,groupi,groupj,mm,total_point
      LOGICAL oka
      INCLUDE 'pbc.h'
      GET_IND(y)=INT(y)+(SIGN(1.0D0,y-INT(y))-1.0D0)/2+1


*----------------------- EXECUTABLE STATEMENTS ------------------------*
      
      fact=1.0D0/(2.0D0**(1.0D0/6.0D0))

*=======================================================================
*--- Determines how many grid cells must be considered to include ------ 
*--- each Voronoi cell. All grid point included in a cube of      ------ 
*--- 2*bin_size of side are included                              ------ 
*=======================================================================

      binx=2.0D0/DFLOAT(nx)
      biny=2.0D0/DFLOAT(ny)
      binz=2.0D0/DFLOAT(nz)

      xc=co(1,1)*binx+co(1,2)*biny+co(1,3)*binz
      yc=             co(2,2)*biny+co(2,3)*binz
      zc=                          co(3,3)*binz

      idxm=GET_IND(bin_size/xc)
      idym=GET_IND(bin_size/yc)
      idzm=GET_IND(bin_size/zc)

      map=0
      idx=0
      idy=0
      idz=0

*=======================================================================
*---- See if dimensions of maxcav are sufficients ----------------------
*=======================================================================

      DO ia=idx-idxm,idx+idxm
         DO ib=idy-idym,idy+idym
            DO ic=idz-idzm,idz+idzm
               map=map+1
            END DO
         END DO
      END DO
      IF(map .GT. maxcav) THEN
         iret=1
         WRITE(errmsg,'('' Dimension of _MAX_CAVITIES_ too small.'',
     &        '' Must be dimensioned .GE. to '',i8)') map
         RETURN
      END IF
      total_point=map

*=======================================================================      
*--- First loop on the Voronoi cells of the atoms of the system --------
*=======================================================================      

      count=0
      o=0
      DO i=nstart,nend
         o=o+1
         m=nnlpp_vor(1,o)
         oka=.TRUE.
c$$$         DO ia=1,nbone
c$$$            IF(mback(ia) .EQ. i) oka=.FALSE.
c$$$         END DO
         IF(m .NE. 0 .AND. oka) THEN
            xpi=xp0(i)
            ypi=yp0(i)
            zpi=zp0(i)
            groupi=atomg(i)
            xpgi=xpg(groupi)
            ypgi=ypg(groupi)
            zpgi=zpg(groupi)
            typei=nbtype(i)
            sigi=sigma(typei)*fact
            sigi=sigi**2

            dx=xpi/binx
            dy=ypi/biny
            dz=zpi/binz

            idx=GET_IND(dx)
            idy=GET_IND(dy)
            idz=GET_IND(dz)
            map=0
            tot_map=0

*.......................................................................
*--- Check if any of the grid points are at less than sigma for the ----
*--- atom and eliminate those that are.                             ----
*.......................................................................

            DO ia=idx-idxm,idx+idxm
               DO ib=idy-idym,idy+idym
                  DO ic=idz-idzm,idz+idzm
                     xd=xpi-binx*DFLOAT(ia)
                     yd=ypi-biny*DFLOAT(ib)
                     zd=zpi-binz*DFLOAT(ic)
                     xmap=2.0D0*PBC(xd)
                     ymap=2.0D0*PBC(yd)
                     zmap=2.0D0*PBC(zd)
                     xd=xd-xmap
                     yd=yd-ymap
                     zd=zd-zmap
                     xc=co(1,1)*xd+co(1,2)*yd+co(1,3)*zd
                     yc=           co(2,2)*yd+co(2,3)*zd
                     zc=                      co(3,3)*zd
                     rs=xc*xc+yc*yc+zc*zc
                     IF(rs .GT. sigi) THEN
                        map=map+1
                        xmap0(map)=binx*DFLOAT(ia)
                        ymap0(map)=biny*DFLOAT(ib)
                        zmap0(map)=binz*DFLOAT(ic)
                        xmap1(map)=xc
                        ymap1(map)=yc
                        zmap1(map)=zc
                        ok(map)=.TRUE.
                     END IF
                  END DO
               END DO
            END DO
*.......................................................................
*--- Loop over the Voronoi contact of the atom i-th. Eliminate     -----
*--- all points of the grid which are not in the Voronoi cell.      ----
*.......................................................................
            DO ja=1,m
               j=nnlpp_vor(ja+1,o)
               groupj=atomg(j)
               xgg=xpgi-xpg(groupj)
               ygg=ypgi-ypg(groupj)
               zgg=zpgi-zpg(groupj)
               xmap=2.0D0*PBC(xgg)
               ymap=2.0D0*PBC(ygg)
               zmap=2.0D0*PBC(zgg)
               xg=xpi-xp0(j)-xmap
               yg=ypi-yp0(j)-ymap
               zg=zpi-zp0(j)-zmap
               xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
               yc=           co(2,2)*yg+co(2,3)*zg
               zc=                      co(3,3)*zg
               rsp=DSQRT(xc*xc+yc*yc+zc*zc)
               rs=rsp*0.5D0
               DO k=1,map
                  IF(ok(k)) THEN
                     xd=xmap1(k)
                     yd=ymap1(k)
                     zd=zmap1(k)
                     dot=(xd*xc+yd*yc+zd*zc)/rsp
                     IF(dot .GT. rs) THEN
                        ok(k)=.FALSE.
                     END IF
                  END IF
               END DO
            END DO
            mapp=0
            DO k=1,map
               IF(ok(k)) THEN
                  mapp=mapp+1
                  index(mapp)=k
               END IF
            END DO

*--- Total of point in the Voronoi cell. To be used to compute histogram

            tot_map=total_point+mapp-map

*.......................................................................
*--- Loop over the grid points. Compute distance with the atom i-th ----
*--- and all its Voronoi contact and compare with sigma of each.    ----
*.......................................................................

            DO ka=1,mapp
               k=index(ka)
               xc=xmap1(k)
               yc=ymap1(k)
               zc=zmap1(k)
               rsp=DSQRT(xc**2+yc**2+zc**2)
               rmap(1)=rsp-sigma(typei)*fact
               xg=xmap0(k)
               yg=ymap0(k)
               zg=zmap0(k)
               DO ja=1,m
                  j=nnlpp_vor(ja+1,o)
                  typej=nbtype(j)
                  xd=xg-xp0(j)
                  yd=yg-yp0(j)
                  zd=zg-zp0(j)
                  xmap=2.0D0*PBC(xd)
                  ymap=2.0D0*PBC(yd)
                  zmap=2.0D0*PBC(zd)
                  xd=xd-xmap
                  yd=yd-ymap
                  zd=zd-zmap
                  xc=co(1,1)*xd+co(1,2)*yd+co(1,3)*zd
                  yc=           co(2,2)*yd+co(2,3)*zd
                  zc=                      co(3,3)*zd
                  rsp=DSQRT(xc**2+yc**2+zc**2)
                  rmap(ja+1)=rsp-sigma(typej)*fact
               END DO
               CALL indexx(m+1,rmap,indexa)
               cavity_r(count+ka)=rmap(indexa(1))
            END DO
*--- cavity(1,i) contains the grid points not overlapping with atoms
            cavity_n(1,o)=mapp
*--- cavity(2,i) contains the grid points contained in the Voronoi cell
            cavity_n(2,o)=tot_map
            count=count+mapp
            IF(count .GT. mqq) THEN
               iret=1
               WRITE(errmsg,'('' _MAX_CAVITIES_ATOM_ too small ''
     &              ,'' Grid dimensios must be .GE. '',i8)') count
               RETURN
            END IF
         ELSE
            cavity_n(1,o)=0
            cavity_n(2,o)=0
         END IF
      END DO

*=======================================================================      
*--- Construct the histogram. dr_histo is the bin size and      --------
*--- rmax_histo is the max radial extension of the histogram    --------
*=======================================================================      

      dx=dr_histo
      nmax=GET_IND(rmax_histo/dx)

*=======================================================================      
*--- Accumulate histograms of cavities. For solute save histograms  ----
*--- of each residue.                                               ----
*=======================================================================      

      count=0
      o=0
      DO i=nstart,nend
         o=o+1
         mapp=cavity_n(1,o)
         IF(ss_index(i) .EQ. 1) THEN
            tot_map=cavity_n(2,o)
            j=res(i)
            DO m=1,mapp
               ibin=GET_IND(cavity_r(count+m)/dx)
               IF(ibin .GT. 0)     cavity_h(1+ibin,j)=cavity_h(1+ibin,j)
     &              +1.0D0
            END DO
            cavity_h(1,j)=cavity_h(1,j)+DFLOAT(tot_map)
         END IF
         count=count+mapp
      END DO
      count=0

*=======================================================================
*--  For solvent accumulate only one histogram -------------------------
*=======================================================================

      o=0
      DO i=nstart,nend
         o=o+1
         mapp=cavity_n(1,o)
         IF(ss_index(i) .EQ. 2) THEN
            tot_map=cavity_n(2,o)
            DO m=1,mapp
               ibin=GET_IND(cavity_r(count+m)/dx)
               IF(ibin .GT. 0) cavity_h(1+ibin,jmax)=cavity_h(1+ibin
     &              ,jmax)+1.0D0
            END DO
            cavity_h(1,jmax)=cavity_h(1,jmax)+DFLOAT(tot_map)
         END IF
         count=count+mapp
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
