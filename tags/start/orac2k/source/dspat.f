      SUBROUTINE dspat(co,oc,xp0,yp0,zp0,asitp,bsitp,nasitp,nbsitp,
     x           nrigg,na,iret,errmsg)

************************************************************************
*                                                                      *
*     Displace the secondary sites of rigid fragments to the           *
*     geometrically appropriate position. For instance, in a           *
*     planar fragment the secondari sites are moved to the plane       *
*     if they happened not to be exactly on it.                        *
*                                                                      *
*                                                                      *
*---- Last update 11/13/92 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi CECAM Universite Paris-Sud 1992        *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER nrigg,na,iret
      REAL*8  co(3,3),oc(3,3),xp0(*),yp0(*),zp0(*)
      INTEGER asitp(4,*),bsitp(na,*),nasitp(*),nbsitp(*)
      CHARACTER*80 errmsg

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,k,m,n,nbsat,nscat,nato,oa,np,ia,ib,ic
      REAL*8  eps,aux1,dummy,rcond,det,xi,yi,zi,xbi,ybi,zbi,xai,yai,zai
      REAL*8  ra,cb,norm,sgns,ver(3),xbd,ybd,zbd,xpi,ypi,zpi,pi,tol,
     x        angle
      DATA tol/23.0D0/

*-------------------- WORKING ARRAYS -----------------------------------

      INCLUDE 'parst.h'
      PARAMETER (oa = m13+4)
      REAL*8  bx(oa),by(oa),bz(oa),abx(oa),aby(oa),abz(oa),abr(oa),
     x        cx(oa),cy(oa),cz(oa)
      COMMON /rag1/ bx,by,bz,abx,aby,abz,abr,cx,cy,cz
      INCLUDE 'pbc.h'

*==================== EXECUTABLE STATEMENTS ============================


!=======================================================================
!------ Calculate the projection matrix for secondary sites ------------
!=======================================================================

      pi=4.0D0*DATAN(1.0D0)
      DO np=1,nrigg
          nscat=nbsitp(np)
          nbsat=nasitp(np)
          nato=nscat+nbsat

!=======================================================================
!------ Check on number of primary sites if 4 go to the next fragment --
!=======================================================================

          IF(nbsat .GT. 3) GOTO 100

          DO ia=1,nbsat
              i=asitp(ia,np)
              xpi=xp0(i)
              ypi=yp0(i)
              zpi=zp0(i)
              bx(ia)=co(1,1)*xpi+co(1,2)*ypi+co(1,3)*zpi
              by(ia)=co(2,1)*xpi+co(2,2)*ypi+co(2,3)*zpi
              bz(ia)=co(3,1)*xpi+co(3,2)*ypi+co(3,3)*zpi
          END DO
          DO ia=1,nscat
              i=bsitp(ia,np)
              ib=ia+nbsat
              xpi=xp0(i)
              ypi=yp0(i)
              zpi=zp0(i)
              bx(ib)=co(1,1)*xpi+co(1,2)*ypi+co(1,3)*zpi
              by(ib)=co(2,1)*xpi+co(2,2)*ypi+co(2,3)*zpi
              bz(ib)=co(3,1)*xpi+co(3,2)*ypi+co(3,3)*zpi
          END DO
          xi=bx(1)
          yi=by(1)
          zi=bz(1)
          n=nbsat-1
          ic=asitp(1,np)
          DO i=1,n
              ia=asitp(i+1,np)
              xpi=2.0D0*PBC(xp0(ia)-xp0(ic))
              ypi=2.0D0*PBC(yp0(ia)-yp0(ic))
              zpi=2.0D0*PBC(zp0(ia)-zp0(ic))
              cx(i)=co(1,1)*xpi+co(1,2)*ypi+co(1,3)*zpi
              cy(i)=co(2,1)*xpi+co(2,2)*ypi+co(2,3)*zpi
              cz(i)=co(3,1)*xpi+co(3,2)*ypi+co(3,3)*zpi
              abx(i)=bx(i+1)-xi-cx(i)
              aby(i)=by(i+1)-yi-cy(i)
              abz(i)=bz(i+1)-zi-cz(i)
              abr(i)=abx(i)**2+aby(i)**2+abz(i)**2
          END DO
          DO i=1,nscat
              ib=bsitp(i,np)
              xpi=2.0D0*PBC(xp0(ib)-xp0(ic))
              ypi=2.0D0*PBC(yp0(ib)-yp0(ic))
              zpi=2.0D0*PBC(zp0(ib)-zp0(ic))
              cx(i+n)=co(1,1)*xpi+co(1,2)*ypi+co(1,3)*zpi
              cy(i+n)=co(2,1)*xpi+co(2,2)*ypi+co(2,3)*zpi
              cz(i+n)=co(3,1)*xpi+co(3,2)*ypi+co(3,3)*zpi
              abx(i+n)=bx(nbsat+i)-xi-cx(i+n)
              aby(i+n)=by(nbsat+i)-yi-cy(i+n)
              abz(i+n)=bz(nbsat+i)-zi-cz(i+n)
              abr(i+n)=abx(i+n)**2+aby(i+n)**2+abz(i+n)**2
          END DO
          IF(n .EQ. 1) THEN
              DO i=1,nscat
                  cb=abx(1)*abx(n+i)+aby(1)*aby(n+i)+abz(1)*abz(n+i)
                  angle=DACOS(cb/DSQRT(abr(1)*abr(n+i)))*180.0D0/pi
                  IF(angle .GT. tol) THEN
                      iret=1
                      errmsg='In DSPAT: One of the fragments is'
     x //' far from being linear. Cannot adjust. Abort.'
                      RETURN
                  END IF
                  bx(nbsat+i)=abx(1)*cb/abr(1)
                  by(nbsat+i)=aby(1)*cb/abr(1)
                  bz(nbsat+i)=abz(1)*cb/abr(1)
                  bx(nbsat+i)=bx(nbsat+i)+xi+cx(i+n)
                  by(nbsat+i)=by(nbsat+i)+yi+cy(i+n)
                  bz(nbsat+i)=bz(nbsat+i)+zi+cz(i+n)
              END DO
              DO ia=1,nscat
                  i=bsitp(ia,np)
                  xbi=bx(nbsat+ia)
                  ybi=by(nbsat+ia)
                  zbi=bz(nbsat+ia)
                  xbd=2.0D0*PBC(xp0(i))
                  ybd=2.0D0*PBC(yp0(i))
                  zbd=2.0D0*PBC(zp0(i))
                  xp0(i)=oc(1,1)*xbi+oc(1,2)*ybi+oc(1,3)*zbi
                  yp0(i)=oc(2,1)*xbi+oc(2,2)*ybi+oc(2,3)*zbi
                  zp0(i)=oc(3,1)*xbi+oc(3,2)*ybi+oc(3,3)*zbi
                  xp0(i)=xp0(i)+xbd
                  yp0(i)=yp0(i)+ybd
                  zp0(i)=zp0(i)+zbd
              END DO
          ELSE
              ver(1)=aby(1)*abz(2)-abz(1)*aby(2)
              ver(2)=abz(1)*abx(2)-abx(1)*abz(2)
              ver(3)=abx(1)*aby(2)-aby(1)*abx(2)
              norm=0.0D0
              DO i=1,3
                  norm=norm+ver(i)**2
              END DO
              norm=DSQRT(norm)
              ver(1)=ver(1)/norm
              ver(2)=ver(2)/norm
              ver(3)=ver(3)/norm
              DO i=1,nscat
                  cb=ver(1)*abx(n+i)+ver(2)*aby(n+i)+ver(3)*abz(n+i)
                  angle=DACOS(cb)*180.0D0/pi
                  IF(angle .GT. 90.0D0+tol .OR. 
     x                            angle .LT. 90.0D0-tol) THEN
                      WRITE(6,*) angle,bsitp(i,np),np
                      WRITE(6,*)
     x                 DSQRT(abr(1)),DSQRT(abr(2)),DSQRT(abr(3))
                      iret=1
                      errmsg='In DSPAT: One of the fragments is'
     x //' far from being planar. Cannot adjust. Abort.'
                      RETURN
                  END IF
                  bx(nbsat+i)=abx(n+i)-ver(1)*cb
                  by(nbsat+i)=aby(n+i)-ver(2)*cb
                  bz(nbsat+i)=abz(n+i)-ver(3)*cb
                  bx(nbsat+i)=bx(nbsat+i)+xi+cx(i+n)
                  by(nbsat+i)=by(nbsat+i)+yi+cy(i+n)
                  bz(nbsat+i)=bz(nbsat+i)+zi+cz(i+n)
              END DO
              DO ia=1,nscat
                  i=bsitp(ia,np)
                  xbi=bx(nbsat+ia)
                  ybi=by(nbsat+ia)
                  zbi=bz(nbsat+ia)
                  xp0(i)=oc(1,1)*xbi+oc(1,2)*ybi+oc(1,3)*zbi
                  yp0(i)=oc(2,1)*xbi+oc(2,2)*ybi+oc(2,3)*zbi
                  zp0(i)=oc(3,1)*xbi+oc(3,2)*ybi+oc(3,3)*zbi
              END DO
          END IF
100       CONTINUE
      END DO

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
