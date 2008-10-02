!!$/---------------------------------------------------------------------\
!!$   Copyright  © 2006-2007 Massimo Marchi <Massimo.Marchi at cea.fr>   |
!!$                                                                      |
!!$    This software is a computer program named oracDD whose            |
!!$    purpose is to simulate and model complex molecular systems.       |
!!$    The code is written in fortran 95 compliant with Technical        |
!!$    Report TR 15581, and uses MPI-1 routines for parallel             |
!!$    coding.                                                           |
!!$                                                                      |
!!$    This software is governed by the CeCILL license under             |
!!$    French law and abiding by the rules of distribution of            |
!!$    free software.  You can  use, modify and/ or redistribute         |
!!$    the software under the terms of the CeCILL icense as              |
!!$    circulated by CEA, CNRS and INRIA at the following URL            |
!!$    "http://www.cecill.info".                                         |
!!$                                                                      |
!!$    As a counterpart to the access to the source code and rights      |
!!$    to copy, modify and redistribute granted by the license,          |
!!$    users are provided only with a limited warranty and the           |
!!$    software's author, the holder of the economic rights, and         |
!!$    the successive licensors have only limited liability.             |
!!$                                                                      |
!!$    The fact that you are presently reading this means that you       |
!!$    have had knowledge of the CeCILL license and that you accept      |
!!$    its terms.                                                        |
!!$                                                                      |
!!$    You should have received a copy of the CeCILL license along       |
!!$    with this program; if not, you can collect copies on the URL's    |
!!$    "http://www.cecill.info/licences/Licence_CeCILL_V2-en.html"       |
!!$    "http://www.cecill.info/licences/Licence_CeCILL_V2-fr.html"       |
!!$                                                                      |
!!$----------------------------------------------------------------------/
MODULE Neighbors

!!$***********************************************************************
!!$   Time-stamp: <2007-01-19 18:57:55 marchi>                           *
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
  PUBLIC Neighbors_, Neighbors__Particles, ind_xyz, Neighbors__Ind&
       &, Neighbors__Chain, chain_xyz, Head_xyz

  TYPE :: Neighbors__Ind
     INTEGER :: i,j,k
  END TYPE Neighbors__Ind     
  TYPE :: Neighbors__Chain
     INTEGER :: i,j,k
     INTEGER :: p
  END TYPE Neighbors__Chain
  TYPE(Neighbors__Ind), ALLOCATABLE, SAVE :: Ind_xyz(:)
  TYPE(Neighbors__Chain), ALLOCATABLE, SAVE :: Chain_xyz(:)
  INTEGER, ALLOCATABLE, SAVE :: Head_xyz(:)
  INTEGER, SAVE :: ncx,ncy,ncz
CONTAINS
!!$
!!$--- Constructor for Ind_xyz
!!$
  FUNCTION Neighbors_(rcut,nx,ny,nz) RESULT(out)
    LOGICAL :: out
    INTEGER :: nx,ny,nz
    REAL(8) :: rcut

    INTEGER :: vect0(3)
    INTEGER, POINTER :: vect(:)
    REAL(8) :: sqcut,dx,dy,dz,rmin
    INTEGER :: imax,jmax,kmax,i,j,k,istart,jstart,kstart,warnx,warny, warnz&
         &,nxmax, nymax,nzmax,nind,jend,kend

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

    vect0=(/0, 0, 0/) 
    IF(.NOT. Node_()) STOP

    CALL Node__Push(vect0)   

    istart=0
    DO i=istart,ncx-1
       jend=ncy-1
       IF(i == 0) jend=0
       DO j=1-ncy,jend
          kend=ncz-1
          IF(i == 0 .AND. j == 0) kend=0
          DO k=1-ncz,kend
             rmin=dist_ijk(i,j,k,dx,dy,dz)
             IF(rmin < sqcut) then
                vect0=(/i, j, k/) 
                IF(imax < abs(i)) imax=abs(i)
                IF(jmax < abs(j)) jmax=abs(j)
                IF(kmax < abs(k)) kmax=abs(k)
                IF(.NOT. (i == 0 .AND. j == 0 .AND. k == 0)) THEN
                   CALL Node__Push(vect0)
                END IF
             END IF
          END DO
       END DO
    END DO
    nind=Node__Size()
    WRITE(*,*) 'nind ',nind
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
          errmsg_f=TRIM(errmsg_f)//' along x '
       END IF
       IF(warny == 1) THEN
          errmsg_f=TRIM(errmsg_f)//' along y '
       END IF
       IF(warnz == 1) THEN
          errmsg_f=TRIM(errmsg_f)//' along z '
       END IF
       CALL Add_Errors(-1,errmsg_f)
       out=.TRUE. 
    END IF
  CONTAINS
    FUNCTION dist_ijk(ni,nj,nk,dx,dy,dz) RESULT(out)
      REAL(8) :: out
      REAL(8) ::  dx,dy,dz
      INTEGER ::  ni,nj,nk

      REAL(8) ::  d,dmin,dt
      REAL(8) ::  lx,ly,lz
      REAL(8) ::  mx,my,mz
      REAL(8) ::  dmx,dmy,dmz
      REAL(8) ::  msq,dmsq,lambda,s

      INTEGER, PARAMETER :: nv(8,3)=RESHAPE((/&
           & 0, 1, 0, 0, 1, 1, 0, 1 &
           &,0, 0, 1, 0, 1, 0, 1, 1 &
           &,0, 0, 0, 1, 0, 1, 1, 1 &
           &/),(/8, 3/))
      INTEGER ::  i,j,imin,jmin
      INTEGER, SAVE ::  ndtmax=0

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
!!$
!!$--- Constructor for Chain_xyz and Head_xyz
!!$
  FUNCTION Neighbors__Particles(x,y,z) RESULT(out)
    LOGICAL :: out
    REAL(8) :: x(:),y(:),z(:)
    REAL(8) :: x1,y1,z1,dx,dy,dz
    INTEGER :: n,nx,ny,nz,natp,numcell,l
    
    out=.TRUE.
    IF(.NOT. Neighbors__Valid()) THEN
       out=.FALSE.
       errmsg_f='Must construct cell indeces before atomic indeces can be obtained'
       CALL Add_Errors(-1,errmsg_f)
       RETURN
    END IF
    IF(.NOT. ALLOCATED(Head_xyz)) THEN
       ALLOCATE(Head_xyz(ncx*ncy*ncz))
       natp=SIZE(x)
       ALLOCATE(Chain_xyz(natp))
    END IF
    Head_xyz=0
    Chain_xyz (:) % p = 0

!!$=======================================================================
!!$     Compute chain list for the system
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
       Chain_xyz (n) % i=nx
       Chain_xyz (n) % j=ny
       Chain_xyz (n) % k=nz
       numcell=nz+ncz*(ny+ncy*nx)+1
       Chain_xyz (n) % p=Head_xyz(numcell)
       Head_xyz(numcell)=n
       
    END DO
  END FUNCTION Neighbors__Particles
  SUBROUTINE Neighbors__Delete
    DEALLOCATE(Ind_xyz,Chain_xyz,Head_xyz)
  END SUBROUTINE Neighbors__Delete
  FUNCTION Neighbors__Valid() RESULT(out)
    LOGICAL :: out
    out=ALLOCATED(ind_xyz)
  END FUNCTION Neighbors__Valid
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE Neighbors
