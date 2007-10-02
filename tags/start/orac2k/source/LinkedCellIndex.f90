SUBROUTINE LinkeCellIndex(ncx,ncy,ncz,IndX,IndY,IndZ,ctoff,co)
!!$***********************************************************************
!!$   Time-stamp: <01/04/12 19:11:42 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Apr 12 2001 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*

  IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

  INTEGER :: ncx,ncy,ncz
  INTEGER, DIMENSION (:), POINTER :: IndX,IndY,IndZ
  REAL(8)  :: ctoff,co(3,3)

!!$----------------------- VARIABLES IN COMMON --------------------------*

!!$------------------------- LOCAL VARIABLES ----------------------------*

  REAL(8) :: dist_ijk,dx,dy,dz,sqcut,rmin
  INTEGER :: i,j,k,n,istart,jstart,kstart,imax,jmax,kmax,nxmax,nymax & 
       ,nzmax,warnx,warny,warnz,nind,ncxA,ncyA,nczA

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

!!$c----------------------------
!!$c initialization
!!$c----------------------------
  sqcut=ctoff**2
  dx=2.0D0/ncx
  dy=2.0D0/ncy
  dz=2.0D0/ncz
  imax=0
  jmax=0
  kmax=0

!!$c----------------------------
!!$c calculation of the half-sphere of indeces.  First way around:
!!$c----------------------------

  ncxA=1.5*(ctoff/co(1,1))/dx
  ncya=1.5*(ctoff/co(2,2))/dy
  nczA=1.5*(ctoff/co(3,3))/dz
  
  nind=1
  istart=1-ncxA
  DO i=istart,ncxA-1
     jstart=1-ncyA
     DO j=jstart,ncyA-1
        kstart=1-nczA
        DO k=kstart,nczA-1
           rmin=dist_ijk(i,j,k,dx,dy,dz,co)
           IF(rmin.lt.sqcut) THEN
              nind=nind+1
              IF(i.eq.0 .and. j.eq.0 .and. k.eq.0) nind=nind-1
           END IF
        END DO
     END DO
  END DO

  ALLOCATE(IndX(nind),IndY(nind),IndZ(nind))

!!$c----------------------------
!!$c calculation of the half-sphere of indeces.  Do it seriously now.
!!$c----------------------------

  nind=1
  IndX(1)=0
  IndY(1)=0
  IndZ(1)=0
  istart=1-ncxA
  DO i=istart,ncxA-1
     jstart=1-ncyA
     DO j=jstart,ncyA-1
        kstart=1-nczA
        DO k=kstart,nczA-1
           rmin=dist_ijk(i,j,k,dx,dy,dz,co)
           IF(rmin.lt.sqcut) THEN
              nind=nind+1
              IndX(nind)=i
              IndY(nind)=j
              IndZ(nind)=k
              IF(i.eq.0 .and. j.eq.0 .and. k.eq.0) nind=nind-1
           END IF
        END DO
     END DO
  END DO

!!$c----------------------------
!!$c pour eviter que, lors du
!!$c comptage des paires d'atomes
!!$c une cellule n'apparaisse
!!$c deux fois, il faut ajouter
!!$c le test suivant :
!!$c si ncx est pair,
!!$c il faut imax <= ncx/2
!!$c si ncx est impair,
!!$c il faut imax <= (ncx+1)/2
!!$c (idem pour jmax et kmax)
!!$c----------------------------

  nxmax=(ncx+1)/2
  nymax=(ncy+1)/2
  nzmax=(ncz+1)/2
  warnx=0
  warny=0
  warnz=0
  IF(imax.ge.nxmax) warnx=1
  IF(jmax.ge.nymax) warny=1
  IF(kmax.ge.nzmax) warnz=1
  IF(warnx.eq.1 .or. warny.eq.1 .or. warnz.eq.1) THEN
     print*,"des cellules risquent d'etre comptees deux fois"
     print*,"diminuez le cutoff ou :"
     IF(warnx.eq.1) print*,"augmentez ncx"
     IF(warny.eq.1) print*,"augmentez ncy"
     IF(warnz.eq.1) print*,"augmentez ncz"
     STOP
  END IF

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE LinkeCellIndex
