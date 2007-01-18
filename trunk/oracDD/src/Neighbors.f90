MODULE Neighbors

!!$***********************************************************************
!!$   Time-stamp: <2007-01-18 17:27:58 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Jan 18 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program oracDD ----*

  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f
  USE Node
  USE Cell
  IMPLICIT none
  PRIVATE
  PUBLIC Neighbors_, Neighbors__Atoms
  TYPE :: Neighbors__Ind
     INTEGER :: i,j,k
  END TYPE Neighbors__Ind
     
  INTEGER, SAVE :: ncx,ncy,ncz
  TYPE(Neighbors__Ind), ALLOCATABLE, SAVE :: ind_xyz(:)
  INTEGER, ALLOCATABLE, SAVE :: Chainp(:),Cellpi(:),Cellpj(:),Cellpk(:),Headp(:)
CONTAINS
  FUNCTION Neighbors_(rcut,nx,ny,nz) RESULT(out)
    LOGICAL :: out
    INTEGER :: nx,ny,nz
    REAL(8) :: rcut
    INTEGER :: vect0(3)
    INTEGER, POINTER :: vect(:)
    REAL(8) :: sqcut,dx,dy,dz,rmin
    INTEGER :: imax,jmax,kmax,i,j,k,istart,jstart,kstart,nind,warnx,warny, warnz&
         &,nxmax, nymax,nzmax

    IF(ALLOCATED(ind_xyz)) DEALLOCATE(ind_xyz)
    
    out=.TRUE.
    ncx=nx; ncy=ny; ncz=nz    

    sqcut=rcut**2

    dx=2.d0/ncx
    dy=2.d0/ncy
    dz=2.d0/ncz
    imax=0
    jmax=0
    kmax=0
    nind=1

    vect0=(/0, 0, 0/) 
    IF(.NOT. Node_()) STOP


    CALL Node__Push(vect0)   
    istart=1-ncx
    DO i=istart,ncx-1
       jstart=1-ncy
       DO j=jstart,ncy-1
          kstart=1-ncz
          DO k=kstart,ncz-1
             rmin=dist_ijk(i,j,k,dx,dy,dz)
             IF(rmin < sqcut) then
                vect0=(/i, j, k/) 
                IF(imax < abs(i)) imax=abs(i)
                IF(jmax < abs(j)) jmax=abs(j)
                IF(kmax < abs(k)) kmax=abs(k)
                IF(.NOT. (i == 0 .AND. j == 0 .AND. k == 0)) THEN
                   nind=nind+1
                   CALL Node__Push(vect0)
                END IF
             END IF
          END DO
       END DO
    END DO

    nind=Node__Size()
    ALLOCATE(ind_xyz(nind))
    nind=0
    DO WHILE(Node__Pop(vect))
       nind=nind+1
       ind_xyz(nind) % i=vect(1)
       ind_xyz(nind) % j=vect(2)
       ind_xyz(nind) % k=vect(3)
    END DO

    nxmax=(ncx+1)/2
    nymax=(ncy+1)/2
    nzmax=(ncz+1)/2
    warnx=0
    warny=0
    warnz=0
    IF(imax.ge.nxmax) warnx=1
    IF(jmax.ge.nymax) warny=1
    IF(kmax.ge.nzmax) warnz=1
    IF(warnx == 1 .OR. warny == 1 .OR. warnz == 1) THEN
       errmsg_f='Neighbor cells might be counted twice: Lower the&
            & cutoff or increase the No. of cell '
       IF(warnx == 1) THEN
          errmsg_f=TRIM(errmsg_f)//'along x '
       END IF
       IF(warny == 1) THEN
          errmsg_f=TRIM(errmsg_f)//'along y '
       END IF
       IF(warnz == 1) THEN
          errmsg_f=TRIM(errmsg_f)//'along z '
       END IF
       CALL Add_Errors(-1,errmsg_f)
       out=.TRUE. 
    END IF
  CONTAINS
    FUNCTION dist_ijk(ni,nj,nk,dx,dy,dz) RESULT(out)
      REAL(8) :: out
      REAL(8) ::  dx,dy,dz,co(3,3)
      INTEGER ::  ni,nj,nk

      REAL(8) ::  d,dmin,dt
      REAL(8) ::  lx,ly,lz
      REAL(8) ::  mx,my,mz
      REAL(8) ::  dmx,dmy,dmz
      REAL(8) ::  msq,dmsq,lambda,s

      INTEGER ::  nv(8,3)
      INTEGER ::  i,j,imin,jmin
      INTEGER ::  ndtmax


      ndtmax=500

      nv(:,1)=(/0, 0, 0, 0, 1, 1, 1, 1/)
      nv(:,2)=(/0, 0, 1, 1, 0, 0, 1, 1/)
      nv(:,3)=(/0, 1, 0, 1, 0, 1, 0, 1/)


!!$
!!$--- Minimum distance between corners of the cells (0, 0, 0) and 
!!$--- (ni, nj, nk)
!!$

      dmin=1.0D8

      do i=1,8
         do j=1,8
            lx=(ni+nv(j,1)-nv(i,1))*dx
            ly=(nj+nv(j,2)-nv(i,2))*dy
            lz=(nk+nv(j,3)-nv(i,3))*dz

            mx=co(1,1)*lx+co(1,2)*ly+co(1,3)*lz
            my=co(2,1)*lx+co(2,2)*ly+co(2,3)*lz
            mz=co(3,1)*lx+co(3,2)*ly+co(3,3)*lz

            d=mx*mx+my*my+mz*mz

            if(d.lt.dmin) then
               dmin=d
               imin=i
               jmin=j
            endif

         end do
      end do
!!$
!!$--- Check if the minimal distance is not on the edge of the cube 
!!$

      lx=(ni+nv(jmin,1)-nv(imin,1))*dx
      ly=(nj+nv(jmin,2)-nv(imin,2))*dy
      lz=(nk+nv(jmin,3)-nv(imin,3))*dz

      mx=co(1,1)*lx+co(1,2)*ly+co(1,3)*lz
      my=co(2,1)*lx+co(2,2)*ly+co(2,3)*lz
      mz=co(3,1)*lx+co(3,2)*ly+co(3,3)*lz

      msq=mx*mx+my*my+mz*mz
!!$
!!$--- Loop on the three edges
!!$

      do i=1,3
         dt=sign(1.,0.5-nv(imin,i))*dx
         
         dmx=co(1,i)*dt
         dmy=co(2,i)*dt
         dmz=co(3,i)*dt
         
         s=mx*dmx+my*dmy+mz*dmz
         dmsq=dmx*dmx+dmy*dmy+dmz*dmz
         lambda=-s/dmsq
         
         if((lambda.gt.0.) .and. (lambda.lt.1.)) then
            d=msq-s*s/dmsq
            if(d.lt.dmin) dmin=d
         endif
      enddo

      do i=1,3
         dt=sign(1.,0.5-nv(jmin,i))*dx
         
         dmx=co(1,i)*dt
         dmy=co(2,i)*dt
         dmz=co(3,i)*dt
         
         s=mx*dmx+my*dmy+mz*dmz
         dmsq=dmx*dmx+dmy*dmy+dmz*dmz
         lambda=-s/dmsq
         
         if((lambda.gt.0.) .and. (lambda.lt.1.)) then
            d=msq-s*s/dmsq
            if(d.lt.dmin) dmin=d
         endif
      enddo

      out=dmin
    END FUNCTION dist_ijk
  END FUNCTION Neighbors_
  FUNCTION Neighbors__Atoms(x,y,z) RESULT(out)
    LOGICAL :: out
    REAL(8) :: x(:),y(:),z(:)
    REAL(8) :: x1,y1,z1,dx,dy,dz
    INTEGER :: n,nx,ny,nz,natp,numcell
    
    out=.TRUE.
    IF(.NOT. Neighbors__Valid()) THEN
       out=.FALSE.
       errmsg_f='Must construct cell indeces before atomic indeces can be obtained'
       CALL Add_Errors(-1,errmsg_f)
       RETURN
    END IF
    IF(.NOT. ALLOCATED(headp)) THEN
       ALLOCATE(headp(SIZE(Ind_xyz)))
       natp=SIZE(x)
       ALLOCATE(chainp(natp),cellpi(natp),cellpj(natp),cellpk(natp))
    END IF
    headp=0

!!$=======================================================================
!!$     Compute chain list for system
!!$=======================================================================

    dx=2.d0/ncx
    dy=2.d0/ncy
    dz=2.d0/ncz

    DO n=1,natp
       x1=x(n)/dx
       y1=y(n)/dy
       z1=z(n)/dz
       nx=INT(x1)+(SIGN(1.D0,x1-INT(x1))-1.)/2
       ny=INT(y1)+(SIGN(1.D0,y1-INT(y1))-1.)/2
       nz=INT(z1)+(sign(1.d0,z1-int(z1))-1.)/2
       nx=MOD(MOD(nx,ncx)+ncx,ncx)
       ny=MOD(MOD(ny,ncy)+ncy,ncy)
       nz=MOD(MOD(nz,ncz)+ncz,ncz)
       cellpi(n)=nx
       cellpj(n)=ny
       cellpk(n)=nz
       numcell=nz+ncz*(ny+ncy*nx)+1
       chainp(n)=headp(numcell)
       headp(numcell)=n
    END DO
    

  END FUNCTION Neighbors__Atoms
  FUNCTION Neighbors__Valid() RESULT(out)
    LOGICAL :: out
    out=ALLOCATED(ind_xyz)
  END FUNCTION Neighbors__Valid
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE Neighbors
