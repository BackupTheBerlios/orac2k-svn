!!$/---------------------------------------------------------------------\
!!$   Copyright  � 2006-2007 Massimo Marchi <Massimo.Marchi at cea.fr>   |
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
MODULE NeighCells
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Jan 25 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*
  
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f
  USE Node
  USE Cell
  USE PI_Communicate
  IMPLICIT none
  PRIVATE
  PUBLIC NeighCells_
  TYPE :: NeighCell__Ind
     INTEGER :: i,j,k
  END TYPE NeighCell__Ind
  TYPE :: NeighCell__Chain
     INTEGER :: i,j,k
     INTEGER :: p
  END TYPE NeighCell__Chain
  TYPE(NeighCell__Ind), ALLOCATABLE, SAVE :: Ind_xyz(:)
  TYPE(NeighCell__Chain), ALLOCATABLE, SAVE :: Chain_xyz(:)
  INTEGER, ALLOCATABLE, SAVE :: Head_xyz(:)
  INTEGER, SAVE :: ncx,ncy,ncz
  INTEGER :: nproc,npx,npy,npz
CONTAINS
!!$
!!$--- Constructor
!!$
  FUNCTION NeighCells_(rcut) RESULT(out)
    LOGICAL :: out
    INTEGER :: nx,ny,nz
    REAL(8) :: rcut
    out=.TRUE.
    IF(.NOT. ALLOCATED(PE)) THEN
       errmsg_f='Need to set up communications before NeighCell can be called'
       CALL Add_Errors(-1,errmsg_f)
       out=.FALSE.
       RETURN
    END IF
    CALL PI__GetParameters(nproc,npx,npy,npz)

    nx=NINT(a/rcut)
    ny=NINT(b/rcut)
    nz=NINT(c/rcut)
    IF(DBLE(MOD(nx,npx))/DBLE(npx) <= 0.5D0) THEN
       ncx=nx-MOD(nx,npx)
    ELSE
       ncx=nx+MOD(npx-MOD(nx,npx),npx)
    END IF

    IF(DBLE(MOD(ny,npy))/DBLE(npy) <= 0.5D0) THEN
       ncy=ny-MOD(ny,npy)
    ELSE
       ncy=ny+MOD(npy-MOD(ny,npy),npy)
    END IF

    IF(DBLE(MOD(nz,npz))/DBLE(npz) <= 0.5D0) THEN
       ncz=nz-MOD(nz,npz)
    ELSE
       ncz=nz+MOD(npz-MOD(nz,npz),npz)
    END IF
  END FUNCTION NeighCells_

  FUNCTION NeighCells__Index(rcut) RESULT(out)
    LOGICAL :: out
    REAL(8) :: rcut
    INTEGER :: vect0(3)
    INTEGER, POINTER :: vect(:)
    REAL(8) :: sqcut,dx,dy,dz,rmin
    INTEGER :: imax,jmax,kmax,i,j,k,istart,jstart,kstart,warnx,warny, warnz&
         &,nxmax, nymax,nzmax,nind

    out=.TRUE.

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
  END FUNCTION NeighCells__Index

END MODULE NeighCells
