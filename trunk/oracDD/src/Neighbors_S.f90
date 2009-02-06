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
MODULE Neighbors_S

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
  USE Groups
  USE PI_
  IMPLICIT none
  PRIVATE
  PUBLIC Neighbors_S_, Neighbors_S__Particles, Neighbors_S__Ind&
       &, Neighbors_S__Chain, Neighbors_S__, clst, Neighbors_S_Check&
       &,Neighbors_S__nc, nc, Neighbors_S__Delete, Chain_xyz,&
       & Head_xyz,clst_pme

  TYPE :: Neighbors_S__Ind
     INTEGER :: i,j,k
  END TYPE Neighbors_S__Ind     
  TYPE :: Neighbors_S__Chain
     INTEGER :: i,j,k
     INTEGER :: p
  END TYPE Neighbors_S__Chain
  TYPE :: Neighbors_S__nc
     INTEGER :: x,y,z
     REAL(8) :: rcut
  END TYPE Neighbors_S__nc
  TYPE :: Neighbors_S__
     TYPE(Neighbors_S__Ind), ALLOCATABLE :: Ind_xyz(:)
  END TYPE Neighbors_S__
  TYPE(Neighbors_S__Chain), ALLOCATABLE,SAVE :: Chain_xyz(:)
  INTEGER, ALLOCATABLE,SAVE :: Head_xyz(:)

  TYPE(Neighbors_S__), SAVE, TARGET :: clst(3)
  TYPE(Neighbors_S__), SAVE, TARGET :: clst_pme(3)


  INTEGER, SAVE :: calls=0
  INTEGER, SAVE :: ncx,ncy,ncz
  TYPE(Neighbors_S__nc), SAVE :: nc(3) 
CONTAINS
!!$
!!$--- Constructor for Ind_xyz
!!$
  FUNCTION Neighbors_S_(i_n,rcut,nx,ny,nz) RESULT(out)
    LOGICAL :: out
    INTEGER :: nx,ny,nz,i_n
    REAL(8) :: rcut

    INTEGER :: vect0(3)
    INTEGER, POINTER :: vect(:)=>NULL()
    REAL(8) :: sqcut,dx,dy,dz,rmin
    INTEGER :: imax,jmax,kmax,i,j,k,istart,warnx,warny, warnz&
         &,nxmax, nymax,nzmax,nind,ncount0,kend,jend
    Integer :: iend,jstart,kstart,iv(3)

    calls=calls+1
    i_n=calls
    out=.TRUE.
    IF(calls > SIZE(nc)) THEN
       WRITE(errmsg_f,'(''Neighbor list can only handle '',i2,'' cutoffs'')')
       CALL Add_Errors(-1,errmsg_f)
       out=.FALSE. 
       RETURN
    END IF
    ncx=nx; ncy=ny; ncz=nz
    nc(calls) % x=nx; nc(calls) % y=ny; nc(calls) % z=nz; nc(calls) % rcut=rcut
    IF(ALLOCATED(clst(calls) % ind_xyz)) DEALLOCATE(clst(calls) % ind_xyz)
    IF(ALLOCATED(clst_pme(calls) % ind_xyz)) DEALLOCATE(clst_pme(calls) % ind_xyz)

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

    iend=ncx-1
    jstart=1-ncy
    kstart=1-ncz

    iv(1)=PI_npx
    iv(2)=PI_npy
    iv(3)=PI_npz
    IF(iv(1) == 1 .AND. iv(2) /= 1) THEN
       iend=0
    ELSE IF(iv(1) == 1 .AND. iv(2) == 1) THEN
       iend=0; jstart=0;
    END IF

    istart=0
    ncount0=0
    DO i=istart,iend
       jend=ncy-1
       IF(i == 0) jend=0
       DO j=jstart,jend
          kend=ncz-1
          IF(i == 0 .AND. j == 0) kend=0
          DO k=kstart,kend
             ncount0=ncount0+1
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
    ALLOCATE(clst(calls) % ind_xyz(nind))

    IF(ASSOCIATED(vect)) DEALLOCATE(vect)
    nind=0
    DO WHILE(Node__Pop(vect))
       nind=nind+1
       clst(calls) % ind_xyz(nind) % i=vect(1)
       clst(calls) % ind_xyz(nind) % j=vect(2)
       clst(calls) % ind_xyz(nind) % k=vect(3)
       If(Pi_node_cart == 0) Write(65,*) 2.0D0*Dble(vect(1))&
            &/Dble(ncx),2.0D0*Dble(vect(2))/Dble(ncy),2.0D0&
            &*Dble(vect(3))/Dble(ncz) 
    END DO

    imax=0
    jmax=0
    kmax=0
    vect0=(/0, 0, 0/) 
    IF(.NOT. Node_()) STOP
    CALL Node__Push(vect0)
    istart=1-ncx
    jstart=1-ncy
    ncount0=0
    DO i=1-ncx,ncx-1
       DO j=1-ncy,ncy-1
          DO k=1-ncz,ncz-1
             ncount0=ncount0+1
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
    ALLOCATE(clst_pme(calls) % ind_xyz(nind))
    nind=0
    IF(ASSOCIATED(vect)) DEALLOCATE(vect)
    DO WHILE(Node__Pop(vect))
       nind=nind+1
       clst_pme(calls) % ind_xyz(nind) % i=vect(1)
       clst_pme(calls) % ind_xyz(nind) % j=vect(2)
       clst_pme(calls) % ind_xyz(nind) % k=vect(3)
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
          WRITE(*,*) 'rcut =',rcut
          WRITE(*,*) kmax,nzmax,ncz
       END IF
       CALL Add_Errors(-1,errmsg_f)
       out=.FALSE. 
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
  END FUNCTION Neighbors_S_
  FUNCTION Neighbors_S_Check(i_n,rcut,nx,ny,nz) RESULT(out)
    INTEGER :: out(3)
    INTEGER :: nx,ny,nz,i_n
    REAL(8) :: rcut

    INTEGER :: mfx,mfy,mfz
    REAL(8) :: sqcut,dx,dy,dz,rmin
    INTEGER :: imax,jmax,kmax,i,j,k,istart,jstart,kstart,warnx,warny, warnz&
         &,nxmax, nymax,nzmax,nind,kend

    mfx=0
    mfy=0
    mfz=0
    ncx=nx; ncy=ny; ncz=nz

    sqcut=rcut**2

    dx=2.d0/ncx
    dy=2.d0/ncy
    dz=2.d0/ncz
    imax=0
    jmax=0
    kmax=0

    istart=0
    DO i=istart,ncx-1
       jstart=1-ncy
       IF(i == 0) jstart=0
       DO j=jstart,ncy-1
          kend=ncz-1
          IF(i == 0 .AND. j == 0) kend=0
          DO k=1-ncz,kend
             rmin=dist_ijk(i,j,k,dx,dy,dz)
             IF(rmin < sqcut) then
                IF(imax < abs(i)) imax=abs(i)
                IF(jmax < abs(j)) jmax=abs(j)
                IF(kmax < abs(k)) kmax=abs(k)
             END IF
          END DO
       END DO
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
       IF(warnx == 1) THEN
          mfx=1
       END IF
       IF(warny == 1) THEN
          mfy=1
       END IF
       IF(warnz == 1) THEN
          mfz=1
       END IF
    END IF
    out=(/mfx, mfy, mfz/)
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
  END FUNCTION Neighbors_S_Check
!!$
!!$--- Constructor for Chain_xyz and Head_xyz
!!$
  FUNCTION Neighbors_S__Particles(i_n) RESULT(out)
    LOGICAL :: out
    INTEGER :: i_n
    REAL(8) :: x1,y1,z1,dx,dy,dz
    INTEGER :: n,nx,ny,nz,numcell,l,ntap,ngrp,AtSt,AtEn,mm,count0

    ngrp=SIZE(Groupa)
    out=.TRUE.
    IF(.NOT. Neighbors_S__Valid(i_n)) THEN
       out=.FALSE.
       errmsg_f='1 Must construct cell indeces before atomic indeces can be obtained'
       CALL Add_Errors(-1,errmsg_f)
       RETURN
    END IF
    ncx = nc(i_n) % x; ncy = nc(i_n) % y; ncz = nc(i_n) % z

    IF(ALLOCATED(Head_xyz)) THEN
       DEALLOCATE(Head_xyz,Chain_xyz)
    END IF
    ALLOCATE(Head_xyz(ncx*ncy*ncz))
    ALLOCATE(Chain_xyz(ngrp))

    Head_xyz=0
    Chain_xyz (:) % p = 0

!!$=======================================================================
!!$     Compute chain list for the system
!!$=======================================================================

    dx=2.d0/ncx
    dy=2.d0/ncy
    dz=2.d0/ncz
    count0=0
    DO n=1,ngrp
       IF(Groupa(n) % Knwn <= 0) CYCLE
       count0=count0+(Groupa(n) % AtEn - Groupa(n) % AtSt+1)
       x1=Groupa(n) % xa/dx
       y1=Groupa(n) % ya/dy
       z1=Groupa(n) % za/dz
       nx=INT(x1)+(SIGN(1.D0,x1-INT(x1))-1.0D0)/2
       ny=INT(y1)+(SIGN(1.D0,y1-INT(y1))-1.0D0)/2
       nz=INT(z1)+(sign(1.d0,z1-int(z1))-1.0D0)/2
       nx=MOD(MOD(nx,ncx)+ncx,ncx)
       ny=MOD(MOD(ny,ncy)+ncy,ncy)
       nz=MOD(MOD(nz,ncz)+ncz,ncz)
       numcell=nz+ncz*(ny+ncy*nx)+1
       Chain_xyz (n) % i=nx
       Chain_xyz (n) % j=ny
       Chain_xyz (n) % k=nz
       Chain_xyz (n) % p=Head_xyz(numcell)
       Head_xyz(numcell)=n
    END DO
  END FUNCTION Neighbors_S__Particles

  SUBROUTINE Neighbors_S__Delete
    INTEGER :: i_n
    DO i_n=1,SIZE(clst)
       DEALLOCATE(clst(i_n) % Ind_xyz)
       DEALLOCATE(clst_pme(i_n) % Ind_xyz)
    END DO
  END SUBROUTINE Neighbors_S__Delete
  FUNCTION Neighbors_S__Valid(i_n) RESULT(out)
    INTEGER :: i_n
    LOGICAL :: out
    out=(ALLOCATED(clst(i_n) % ind_xyz) .AND. ALLOCATED(clst_pme(i_n) % ind_xyz))
  END FUNCTION Neighbors_S__Valid
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE Neighbors_S
