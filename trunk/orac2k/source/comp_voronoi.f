      SUBROUTINE comp_voronoi(nstart,nend,nato,xpg,ypg,zpg,atomg
     &     ,xp0,yp0,zp0,co,iret,errmsg)

************************************************************************
*   Time-stamp: <2006-10-24 14:42:07 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Jun 21 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      USE VORONOI_Mod, ONLY: nnlpp_vor,area_vor,volume_vor,maxpla,maxver
     &     ,only_water,rsq_nnlpp_vor
      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER nato,iret,nstart,nend
      INTEGER atomg(*)
      REAL*8  xp0(*),yp0(*),zp0(*),xpg(*),ypg(*),zpg(*),co(3,3)
      CHARACTER*80 errmsg

*----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'parst.h'
      REAL(8), ALLOCATABLE :: pla(:,:),vrt(:,:,:),area(:),d2(:)
#if defined _CRAY_ | defined T3E
      INTEGER(8), ALLOCATABLE ::  nver(:)
#else
      INTEGER, ALLOCATABLE ::  nver(:)
#endif
      INTEGER, ALLOCATABLE ::  iver(:)

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8  rs,xgg,ygg,zgg,xpgi,ypgi,zpgi,xpi,ypi,zpi,xc,yc,zc,xg,yg
     &     ,zg,frac,surface,volume,vol,ax,ay,az,bx,by,bz,cx,cy,cz,abcx
     &     ,abcy,abcz,xmap,ymap,zmap
      INTEGER j,k,map,i1,j1,groupi,groupj,i,o
      INCLUDE 'pbc.h'

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      ALLOCATE(nver(maxpla),iver(maxpla))
      ALLOCATE(pla(4,maxpla),vrt(3,maxver,maxpla),area(maxpla),d2(maxpla
     &     ))

      iret=0
      frac=0.5D0
      o=0
      DO i1=nstart,nend
         o=o+1
         surface=0.0D0
         volume=0.0D0
         map=0
         IF(nnlpp_vor(1,o) .EQ. 0) goto 100
         xpi=xp0(i1)
         ypi=yp0(i1)
         zpi=zp0(i1)
         groupi=atomg(i1)
         xpgi=xpg(groupi)
         ypgi=ypg(groupi)
         zpgi=zpg(groupi)
         
         DO j=1,maxpla
            j1=nnlpp_vor(j,o)
            groupj=atomg(j1)
            xgg=xpgi-xpg(groupj)
            ygg=ypgi-ypg(groupj)
            zgg=zpgi-zpg(groupj)
            xmap=2.0D0*PBC(xgg)
            ymap=2.0D0*PBC(ygg)
            zmap=2.0D0*PBC(zgg)
            xg=xpi-xp0(j1)-xmap
            yg=ypi-yp0(j1)-ymap
            zg=zpi-zp0(j1)-zmap
            xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
            yc=           co(2,2)*yg+co(2,3)*zg
            zc=                      co(3,3)*zg
            rs=xc*xc+yc*yc+zc*zc
            d2(j)=rs
            pla(1,j)=frac*xc
            pla(2,j)=frac*yc
            pla(3,j)=frac*zc
            pla(4,j)=frac*frac*rs
            nver(j)=0
            iver(j)=j1
         END DO
         CALL vstart(pla,vrt,nver)

         DO j=1,maxpla
            IF(nver(j) .GT. 0) THEN            
               ax=vrt(1,1,j)
               ay=vrt(2,1,j)
               az=vrt(3,1,J)
               vol = 0.0
               DO k=2,nver(j)-1
                  bx   = vrt(1,k,j)
                  by   = vrt(2,k,j)
                  bz   = vrt(3,k,j)
                  cx   = vrt(1,k+1,j)
                  cy   = vrt(2,k+1,j)
                  cz   = vrt(3,k+1,j)
                  abcx = ax*(by*cz-bz*cy)
                  abcy = ay*(bz*cx-bx*cz)
                  abcz = az*(bx*cy-by*cx)
                  vol  = vol + DABS(abcx+abcy+abcz)
               END DO
               map         = map + 1
               iver(map)   = iver(j)
               nver(map)   = nver(j)
               area(map)   = vol / SQRT(d2(j))
               volume      = volume+vol
               surface     =surface+area(map)
            END IF
         END DO
100      CONTINUE
         nnlpp_vor(1,o)=map
         area_vor(1,o)=surface
         volume_vor(i1)=volume/6.0D0
         DO j=1,map
            nnlpp_vor(j+1,o)=iver(j)
            area_vor(j+1,o)=area(j)
         END DO
      END DO
      DEALLOCATE(nver,iver)
      DEALLOCATE(pla,vrt,area,d2)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
