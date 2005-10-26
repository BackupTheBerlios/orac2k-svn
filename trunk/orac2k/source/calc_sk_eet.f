      SUBROUTINE calc_sk_eet(xp0,yp0,zp0,whe,nstart,nend,rkcut,delsk
     &     ,sofk,nsofk,psofk)

************************************************************************
*                                                                      *
*                                                                      *
*                                                                      *
*---- Last update 01/07/93 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi CECAM, Orsay France 1993               *
*                                                                      *
*     EXTERNALS XERROR.                                                *
*                                                                      *
************************************************************************

C======================= DECLARATIONS ==================================

      IMPLICIT none

C----------------------- ARGUMENTS -------------------------------------

      INTEGER nstart,nend,psofk
      REAL*8  rkcut,delsk
      REAL*8  xp0(*),yp0(*),zp0(*),sofk(psofk),whe(*),nsofk(psofk)

C----------------------- PARAMETERS ------------------------------------

C-------------------- VARIABLES IN COMMONS -----------------------------
      
      INCLUDE 'unit.h'
      INCLUDE 'parst.h'
      REAL*8 ri0(m1),cthi0(m1),phi0(m1),ind(m1)
      COMMON /rag1/ ri0,cthi,phi0,ind

C-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,l,m,n,k,ii,jj
      REAL*8  rk,ik,gmn_r,gmn_i,ri,xc,yc,zc,cthi,phi,kr,plgndr,sj,sy,sjp
     &     ,temp,syp,yrnm,yinm,weight,eps,gmn,delski,sqpi,jn,xpi,ypi,zpi
      INTEGER kbox,nmax,nkcut,nweight,map
      INTEGER*8 n_a
      LOGICAL near0

C================ EXECUTABLE STATEMENTS ================================

      delski=1.0D0/delsk
      nkcut=rkcut/delsk+1
      nmax=rkcut+2.5D0
      nweight=0
      map=0
      DO i=nstart,nend
         IF(.NOT. near0(whe(i))) THEN
            map=map+1
            ind(map)=i
         END IF
         nweight=nweight+IDINT(whe(i))
      END DO
      weight=1.0D0/DBLE(nweight)
      sqpi=DSQRT(pi)
      
      DO ii=1,map
         i=ind(ii)
         xc=xp0(i)
         yc=yp0(i)
         zc=zp0(i)
         ri=DSQRT(xc*xc+yc*yc+zc*zc)
         ri0(ii)=ri
         cthi0(ii)=zc/ri
         phi0(ii)=DATAN(yc/xc)
      END DO

      DO k=1,nkcut
         rk=DBLE(k)*delsk
         ik=0.0D0
         nmax=2.0D0*pi*rk+4
         DO n=0,nmax
            n_a=n
            DO m=0,n
               gmn_r=0.0D0
               gmn_i=0.0D0
               DO i=1,map
                  ri=ri0(i)
                  cthi=cthi0(i)
                  phi=phi0(i)
                  kr=ri*rk*pi
                  CALL sphbes(n,kr,sj,sy,sjp,syp)
                  temp=plgndr(n,m,cthi)*sj
                  yrnm=DCOS(m*phi)*temp
                  yinm=DSIN(m*phi)*temp
                  gmn_r=gmn_r+yrnm
                  gmn_i=gmn_i+yinm
               END DO
               gmn=gmn_i**2+gmn_r**2
               eps=2.0D0
               IF(m .EQ. 0) eps=1.0D0
               ik=ik+0.5D0*eps*gmn
            END DO
         END DO
         kbox=MIN0(IDINT(delski*rk+0.5D0),psofk)
         sofk(kbox)=sofk(kbox)+ik*weight
         nsofk(kbox)=nsofk(kbox)+1.0D0
         WRITE(6,'(2e17.7)') rk,ik*weight
         WRITE(78,'(2e17.7)') rk,ik*weight
         ik=0.0D0
         DO ii=1,map
            i=ind(ii)
            xpi=xp0(i)
            ypi=yp0(i)
            zpi=zp0(i)
            DO jj=ii+1,map
               j=ind(jj)
               xc=xp0(i)-xp0(j)
               yc=yp0(i)-yp0(j)
               zc=zp0(i)-zp0(j)
               ri=DSQRT(xc*xc+yc*yc+zc*zc)
               kr=rk*ri*pi
               ik=ik+DSIN(kr)/kr
            END DO
         END DO
         WRITE(6,'(2e17.7)') rk,ik*weight
         WRITE(79,'(2e17.7)') rk,ik*weight*2.0D0
      END DO


C================ END EXECUTABLE STATEMENTS ============================

      STOP
      END
