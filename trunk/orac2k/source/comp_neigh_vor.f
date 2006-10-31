      SUBROUTINE comp_neigh_vor(nstart_h,nend_h,nstart_ah,nend_ah
     &     ,heavy_vor,beta,nato,ngrp,grppt,cutoff_vor,xp0,yp0,zp0
     &     ,xpg,ypg,zpg,co,iret,errmsg)

************************************************************************
*   Time-stamp: <2006-10-24 13:41:53 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed Jun 18 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      USE VORONOI_Mod, ONLY: nnlpp_vor,rsq_nnlpp_vor,nsq_nnlpp_vor
     &     ,area_vor,volume_vor,maxpla,max_neigh,compress,solvent
     &     ,dist_min,sigmas,dist_id,wat_mol,watp,wat,Res_Class
     &     ,only_water
      USE Module_Neighbors, ONLY: Neigh_Start=>Start, Neigh_Delete
     &     =>Delete, Neighbors, neigha, neighb, neighc
      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER nig_nnl,nnnlpp_vor,iret,nstart_h,nend_h,nstart_ah,nend_ah
      INTEGER nato,ngrp,grppt(2,*)
      REAL*8  xp0(*),yp0(*),zp0(*),xpg(*),ypg(*),zpg(*),co(3,3)
     &     ,cutoff_vor
      LOGICAL heavy_vor
      CHARACTER*7 beta(*)
      CHARACTER*80 errmsg

*----------------------- VARIABLES IN SCRATCH COMMON ------------------*

      INCLUDE 'parst.h'

      INTEGER index(max_neigh+1)
      INTEGER indexa(max_neigh),i_min
      REAL*8    rsq(max_neigh),rsqq(max_neigh),xmap(m1),ymap(m1),zmap(m1
     &     )

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8  rs,xgg,ygg,zgg,xpgi,ypgi,zpgi,xpi,ypi,zpi,xc,yc,zc,xg,yg
     &     ,zg,cutoff2,rmin,rsp,xmin,ymin,zmin,xb,yb,zb
      INTEGER map,i,i1,j1,j,m,n,jj,mz,o,p_nn,kk,kl
      LOGICAL ok
      TYPE(Neighbors), DIMENSION (:), POINTER :: neigh

      INCLUDE 'pbc.h'

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      neigh=>neigha
      iret=0
      cutoff2=cutoff_vor**2
      n=0
      o=0
      p_nn=0
      DO i=nstart_h,nend_h
         p_nn=p_nn+1
         m=neigh(p_nn) % no
         xpgi=xpg(i)
         ypgi=ypg(i)
         zpgi=zpg(i)
         DO jj=1,m
            j=neigh(p_nn) % nb(jj)
            xgg=xpgi-xpg(j)
            ygg=ypgi-ypg(j)
            zgg=zpgi-zpg(j)
            xmap(jj)=2.0D0*PBC(xgg)
            ymap(jj)=2.0D0*PBC(ygg)
            zmap(jj)=2.0D0*PBC(zgg)
         END DO
         IF(heavy_vor) THEN
            DO i1=grppt(1,i),grppt(2,i)
               index(1)=0
               o=o+1
               dist_min(i1)=1.0D9
               IF(beta(i1)(1:1) .NE. 'h') THEN
                  xpi=xp0(i1)
                  ypi=yp0(i1)
                  zpi=zp0(i1)
                  map=0
                  rmin=1.0D9
                  i_min=-100
                  ok=.FALSE.
                  DO jj=1,m
                     j=neigh(p_nn) % nb(jj)
                     DO j1=grppt(1,j),grppt(2,j)
                        IF(i1 .NE. j1 .AND. beta(j1)(1:1) .NE. 'h') THEN
                           xg=xpi-xp0(j1)-xmap(jj)
                           yg=ypi-yp0(j1)-ymap(jj)
                           zg=zpi-zp0(j1)-zmap(jj)
                           xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
                           yc=           co(2,2)*yg+co(2,3)*zg
                           zc=                      co(3,3)*zg
                           rs=xc*xc+yc*yc+zc*zc
                           
                           IF(solvent(i1) .AND. (.NOT. solvent(j1)))
     &                          THEN
                              rsp=SQRT(rs)-sigmas(j1)
                              IF(rmin .GT. rsp) THEN
                                 rmin=rsp
                                 i_min=j1
                                 xmin=xc
                                 ymin=yc
                                 zmin=zc
                                 ok=.TRUE.
                              END IF
                           END IF
                           IF(rs .LT. cutoff2) THEN
                              map=map+1
                              IF(map .GT. max_neigh) THEN
                                 iret=1
                                 errmsg=
     &' While counting neighbors, _MAX_NEIGH_ in config.h is '
     &//'too small. Abort'
                                 WRITE(6,*) map,max_neigh
                                 RETURN
                              END IF
                              rsq(map)=rs
                              rsqq(map)=SQRT(rs)-sigmas(j1)
                              index(map+1)=j1
                           END IF
                        END IF
                     END DO
                  END DO                  
                  index(1)=map
                  IF(ok) THEN
                     dist_min(i1)=rmin 
                     dist_id(i1)=i_min

                     rsp=SQRT(xmin**2+ymin**2+zmin**2)
                     xc=xmin/rsp
                     yc=ymin/rsp
                     zc=zmin/rsp
                     kk=wat_mol(i1)
                     IF(kk == 0) THEN
                        WRITE(*,*) ' 0 ',kk
                        STOP
                     END IF
                     IF(kk > SIZE(watp)) THEN
                        WRITE(*,*) ' > SIZE(watp) ',kk,i1
                        STOP
                     END IF
c$$$
c$$$ ----  Added for Ions
c$$$                    
                     IF(Res_Class(i_min) >= 1) THEN
                        watp(kk)%o=Res_Class(i_min)
                        watp(kk)%rsp=rmin
                        IF(kk .GT. 0) THEN
                           DO kl=1,4
                              xb=wat(kk)%v(1,kl)
                              yb=wat(kk)%v(2,kl)
                              zb=wat(kk)%v(3,kl)
                              watp(kk)%cosa(kl)=xb*xc+yb*yc+zb*zc+1.0D0
                           END DO
                        END IF
                     END IF
                  END IF
               END IF

*=======================================================================
*--- Order the neighbor atoms according to the ascending order ---------
*--- of distance from the atom i ---------------------------------------
*=======================================================================

               mz=index(1)
               IF(only_water) THEN
                  nsq_nnlpp_vor(:,o)=0
                  rsq_nnlpp_vor(:,o)=0.0D0
               END IF
               IF(mz .NE. 0) THEN
                  CALL indexx(mz,rsq,indexa)
*=======================================================================
*--- Take only the first max_neigh -------------------------------------
*=======================================================================
                  IF(only_water) THEN
                     IF(mz > max_neigh) THEN
                        WRITE(*,*) 'Increase max_neigh, please.'
     &                       ,' Program stops'
                        STOP
                     END IF
                     nsq_nnlpp_vor(1,o)=mz

                     DO jj=1,mz
                        nsq_nnlpp_vor(jj+1,o)=index(1+indexa(jj))
                        rsq_nnlpp_vor(jj+1,o)=rsqq(indexa(jj))
                     END DO
                  END IF
                  
                  IF(mz .GT. maxpla) THEN
                     DO jj=1,maxpla
                        nnlpp_vor(jj,o)=index(1+indexa(jj))
                     END DO
                  ELSE
                     DO jj=1,mz
                        nnlpp_vor(jj,o)=index(1+indexa(jj))
                     END DO
                     DO jj=mz+1,maxpla
                        nnlpp_vor(jj,o)=index(1+indexa(mz))
                     END DO
                  END IF
               ELSE
                  DO jj=1,maxpla
                     nnlpp_vor(jj,o)=0
                  END DO
               END IF
               IF(mz .NE. 0 .AND. mz .LT. maxpla) THEN
                  DO j=mz+1,maxpla
                     indexa(j)=indexa(mz)
                  END DO
               END IF
            END DO
         ELSE
            DO i1=grppt(1,i),grppt(2,i)
               o=o+1
               xpi=xp0(i1)
               ypi=yp0(i1)
               zpi=zp0(i1)
               map=0
               DO jj=1,m
                  j=neigh(p_nn) % nb(jj)
                  DO j1=grppt(1,j),grppt(2,j)
                     IF(i1 .NE. j1) THEN
                        xg=xpi-xp0(j1)-xmap(jj)
                        yg=ypi-yp0(j1)-ymap(jj)
                        zg=zpi-zp0(j1)-zmap(jj)
                        xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
                        yc=           co(2,2)*yg+co(2,3)*zg
                        zc=                      co(3,3)*zg
                        rs=xc*xc+yc*yc+zc*zc
                        IF(rs .LT. cutoff2) THEN
                           map=map+1
                           rsq(map)=rs
                           rsqq(map)=SQRT(rs)-sigmas(j1)
                           index(map+1)=j1
                        END IF
                     END IF
                  END DO
                  IF(map .GT. max_neigh) THEN
                     iret=1
                     errmsg=
     &' While counting neighbors, _MAX_NEIGH_ in config.h is '
     &//'too small. Abort'
                     RETURN
                  END IF
               END DO
               index(1)=map

*=======================================================================
*--- Order the neighbor atoms according to the ascending order ---------
*--- of distance from the atom i ---------------------------------------
*=======================================================================

               mz=index(1)
               IF(mz .NE. 0) THEN
                  IF(mz .EQ. 1) THEN
                     indexa(1)=1
                  ELSE
                     CALL indexx(mz,rsq,indexa)
                  END IF

                  IF(only_water) THEN
                     nsq_nnlpp_vor(:,o)=0
                     rsq_nnlpp_vor(:,o)=0.0D0
                     nsq_nnlpp_vor(1,o)=mz
                     DO jj=1,mz
                        nsq_nnlpp_vor(jj+1,o)=index(1+indexa(jj))
                        rsq_nnlpp_vor(jj+1,o)=rsqq(jj)
                     END DO
                  END IF

*=======================================================================
*--- Take only the first max_neigh -------------------------------------
*=======================================================================

                  IF(mz .GT. maxpla) THEN
                     DO jj=1,maxpla
                        nnlpp_vor(jj,o)=index(1+indexa(jj))
                     END DO
                  ELSE
                     DO jj=1,mz
                        nnlpp_vor(jj,o)=index(1+indexa(jj))
                     END DO
                     DO jj=mz+1,maxpla
                        nnlpp_vor(jj,o)=index(1+indexa(mz))
                     END DO
                  END IF
               ELSE
                  DO jj=1,maxpla
                     nnlpp_vor(jj,o)=0
                  END DO
               END IF
               IF(mz .NE. 0 .AND. mz .LT. maxpla) THEN
                  DO j=mz+1,maxpla
                     indexa(j)=indexa(mz)
                  END DO
               END IF
            END DO
         END IF
         n=n+m+1
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
