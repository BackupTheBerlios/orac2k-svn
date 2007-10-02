      SUBROUTINE insert_solute(kprint,co,mb,mc,xp0,yp0,zp0,xpa,ypa,zpa
     &     ,xpb,ypb,zpb,index,nbtype,radius,sigma,nbtype_slv,ntap,nmol
     &     ,nato_slv,iret,errmsg)

************************************************************************
*   Time-stamp: <2005-01-29 13:14:53 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Feb 11 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER mb,mc
      INTEGER, ALLOCATABLE :: grid(:)

      INTEGER iret,ntap,nato_slv,nmol,kprint,index(*),nbtype(*)
     &     ,nbtype_slv(*)
      REAL*8  co(3,3),xp0(*),yp0(*),zp0(*),xpa(*),ypa(*),zpa(*),xpb(*)
     &     ,ypb(*),zpb(*),sigma(*),radius
      CHARACTER*80 errmsg

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER m,l,i,ma,na,j,count,k,o,ia,ib,ic,iaa,ibb,icc,idx,idy,idz
     &     ,nmaxx,nmaxy,nmaxz,p,len,GET_IND,IE,M_get_length
      REAL*8  fact,xpi,ypi,zpi,xsp,ysp,zsp,xc,yc,zc,radi,radj,radij,aux
      REAL*8  bin_size,binx,biny,binz,dx,dy,dz,y,OFF_SET
      LOGICAL ok
      INCLUDE 'pbc.h'
      GET_IND(y)=INT(y)+(SIGN(1.0D0,y-INT(y))-1.0D0)/2+1
      OFF_SET(y)=y+1.0D0
      IE(i,j,k,l)=(l-1)*40*nmaxx*nmaxy+(k-1)*40*nmaxx+(j-1)*40+i

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      iret=0
      fact=1.0D0/(2.0D0**(1.0D0/6.0D0))
      WRITE(kprint,1)

*=======================================================================
*---- Set the bin size along the three axis. It assumes that -----------
*---- alpha, beta and gamma are sufficiently close to 90.0   -----------
*=======================================================================

      bin_size=4.0D0

      nmaxx=GET_IND(2.0D0*co(1,1)/bin_size)
      nmaxy=GET_IND(2.0D0*co(2,2)/bin_size)
      nmaxz=GET_IND(2.0D0*co(3,3)/bin_size)

      binx=2.0D0/nmaxx
      biny=2.0D0/nmaxy
      binz=2.0D0/nmaxz

      len=mb*nmaxx*nmaxy*nmaxz
      ALLOCATE(grid(len))
      DO i=1,nmaxx
         DO j=1,nmaxy
            DO k=1,nmaxz
               grid(IE(1,i,j,k))=0
            END DO
         END DO
      END DO

      DO i=1,ntap
         xpi=xp0(i)
         ypi=yp0(i)
         zpi=zp0(i)
         xpi=OFF_SET(xpi)
         ypi=OFF_SET(ypi)
         zpi=OFF_SET(zpi)

         dx=xpi/binx
         dy=ypi/biny
         dz=zpi/binz

         idx=GET_IND(dx)
         idy=GET_IND(dy)
         idz=GET_IND(dz)

         IF(idx .LT. 1) idx=idx+nmaxx
         IF(idy .LT. 1) idy=idy+nmaxy
         IF(idz .LT. 1) idz=idz+nmaxz

         IF(idx .GT. nmaxx) idx=idx-nmaxx
         IF(idy .GT. nmaxy) idy=idy-nmaxy
         IF(idz .GT. nmaxz) idz=idz-nmaxz

         grid(IE(1,idx,idy,idz))=grid(IE(1,idx,idy,idz))+1
         count=grid(IE(1,idx,idy,idz))+1
         IF(count .GT. mb) THEN
            iret=1
            errmsg=
     &        ' First dimensions of grid insufficient. Hardcoded' / /
     &        ' to grid(40,*,*,*).'
            RETURN
         END IF
         grid(IE(count,idx,idy,idz))=i
      END DO

      count=0
      DO m=1,nmol
         l=nato_slv*(m-1)
         p=0
         DO i=1,nato_slv
            xpi=xpa(l+i)
            ypi=ypa(l+i)
            zpi=zpa(l+i)
            xpi=OFF_SET(xpi)
            ypi=OFF_SET(ypi)
            zpi=OFF_SET(zpi)

            dx=xpi/binx
            dy=ypi/biny
            dz=zpi/binz

            idx=GET_IND(dx)
            idy=GET_IND(dy)
            idz=GET_IND(dz)
            xpi=xpa(l+i)
            ypi=ypa(l+i)
            zpi=zpa(l+i)

            DO iaa=idx-1,idx+1
               ia=iaa
               IF(iaa .LT. 1)     ia=iaa+nmaxx
               IF(iaa .GT. nmaxx) ia=iaa-nmaxx
               DO ibb=idy-1,idy+1
                  ib=ibb
                  IF(ibb .LT. 1)     ib=ibb+nmaxy
                  IF(ibb .GT. nmaxy) ib=ibb-nmaxy
                  DO icc=idz-1,idz+1
                     ic=icc
                     IF(icc .LT. 1)     ic=icc+nmaxz
                     IF(icc .GT. nmaxz) ic=icc-nmaxz
                     ok=.TRUE.
                     radi=sigma(nbtype_slv(i))*fact
                     o=grid(IE(1,ia,ib,ic))
                     DO k=1,o
                        p=p+1
                        j=grid(IE(1+k,ia,ib,ic))
                        radj=sigma(nbtype(j))*fact
                        xsp=xpi-xp0(j)
                        ysp=ypi-yp0(j)
                        zsp=zpi-zp0(j)
                        xsp=xsp-2.0D0*PBC(xsp)
                        ysp=ysp-2.0D0*PBC(ysp)
                        zsp=zsp-2.0D0*PBC(zsp)
                        xc=co(1,1)*xsp+co(1,2)*ysp+co(1,3)*zsp
                        yc=            co(2,2)*ysp+co(2,3)*zsp
                        zc=                        co(3,3)*zsp
                        xpb(p)=xc**2+yc**2+zc**2
                        radij=(radi+radj)*radius
                        radij=radij*radij
                        IF(xpb(p) .LT. radij) THEN
                           ok=.FALSE.
                           GOTO 100
                        END IF
                     END DO
                  END DO
               END DO
            END DO
         END DO
100      CONTINUE
         
         aux=1.0D10
         ypb(m)=aux
         IF(ok) THEN
            count=count+1
            DO j=1,p
               IF(aux .GT. xpb(j)) aux=xpb(j)
            END DO
            IF(p .EQ. 0) aux=300.0D0
            ypb(m)=aux
         END IF
      END DO

*=======================================================================
*--- Order the solvent molecules according to the ascending order ------
*--- of distance from the solute atoms ---------------------------------
*=======================================================================

      CALL indexx(nmol,ypb,index)

      WRITE(kprint,2) nmol-count,count
      
      nmol=count
      DO m=1,nmol
         ma=(index(m)-1)*nato_slv
         na=(m-1)*nato_slv
         DO j=1,nato_slv
            xpb(na+j)=xpa(ma+j)
            ypb(na+j)=ypa(ma+j)
            zpb(na+j)=zpa(ma+j)
         END DO
      END DO

      DO m=1,nmol*nato_slv
         xpa(m)=xpb(m)
         ypa(m)=ypb(m)
         zpa(m)=zpb(m)
      END DO
      DEALLOCATE(grid)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

1     FORMAT(/ /
     &     '             =============================================='
     &/ /  '             Inserting solute molecules            ------->'
     &/ /)

2     FORMAT(
     &     '             ',i6,' molecules have been removed',i6
     &     ,' remain '/ /)

      RETURN
      END
