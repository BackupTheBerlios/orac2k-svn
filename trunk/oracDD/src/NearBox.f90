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
MODULE NearBox
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

  USE PI_
  USE Constants
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f
  USE Node
  USE Cell
  USE Geometry
  IMPLICIT none
  PRIVATE
  PUBLIC NearBox_,Neighc_, Neighc_pme_, Ind_Large, NearBox__&
       &,NearBox__Neigh, NearBox__Map, NearBox__MapLarge&
       &,NearBox__Planes,Planes,Ind_Small,NearBox__Param

  INTEGER, PARAMETER :: CellMax_=3

  TYPE :: NearBox__Map
     INTEGER :: nx,ny,nz
  END TYPE NearBox__Map

  TYPE :: NearBox__MapLarge
     TYPE(NearBox__Map), ALLOCATABLE :: pt(:)
  END TYPE NearBox__MapLarge
  
  TYPE :: NearBox__Neigh
     INTEGER, ALLOCATABLE :: c(:,:)
  END TYPE NearBox__Neigh
  TYPE :: NearBox__
     TYPE(NearBox__Neigh), ALLOCATABLE :: Nei(:)
  END type NearBox__
  TYPE :: NearBox__Planes
     REAL(8), DIMENSION(4) :: vxy0,vxz0,vyz0,vxy1,vxz1,vyz1
  END type NearBox__Planes

  TYPE(NearBox__), TARGET, SAVE :: Neighc_(CellMax_)
  TYPE(NearBox__), TARGET, SAVE :: Neighc_pme_(CellMax_)
  TYPE(NearBox__Planes), ALLOCATABLE, SAVE :: Planes(:)
  TYPE(NearBox__MapLarge), ALLOCATABLE, SAVE :: Ind_Large(:,:,:)
  INTEGER, ALLOCATABLE, SAVE :: Ind_Small(:,:,:)
  INTEGER, SAVE :: ncx,ncy,ncz
CONTAINS
!!$
!!$--- Constructor
!!$
  FUNCTION NearBox_(rcut0,rcut,nprocs,npax,npay,npaz,i_cut) RESULT(out)
    INTEGER :: npax,npay,npaz,i_cut,nprocs
    LOGICAL :: out
    INTEGER :: nx,ny,nz,mx,my,mz,i,j,k,n,np,count1,o,iv,jv,kv
    REAL(8) :: rcut,rcut0

    CHARACTER(len=max_char) :: labs0
    LOGICAL, ALLOCATABLE :: mask(:,:,:)
    INTEGER :: vec0(3)
    INTEGER, POINTER :: vec(:)=>NULL()
    INTEGER :: npx,npy,npz,k1,k2,k3,il_Max,jl_Max,kl_Max,il_Min&
         &,jl_Min,kl_Min,mp,m,ox,oy,oz,iter_inst
    INTEGER, SAVE ::  iter_max=20
    TYPE :: Neighbors__Ind
       INTEGER :: i,j,k
    END TYPE Neighbors__Ind
    TYPE(Neighbors__Ind), ALLOCATABLE :: Ind_xyz(:)
    REAL(8), SAVE :: dx,dy,dz,ddx,ddy,ddz    
    INTEGER, SAVE :: Calls=0
    TYPE(NearBox__Neigh), POINTER :: Nei(:)
    
    out=.TRUE.
    npx=npax; npy=npay; npz=npaz

!!$
!!$--- If it the first time allocate space for root
!!$

    IF(i_cut > CellMax_) THEN
       errmsg_f='Number of Cutoffs exceed 3. Abort'
       CALL Add_Errors(-1,errmsg_f)
       RETURN
    END IF

    CALL Init
    ALLOCATE(Neighc_(i_cut) % Nei(nprocs))
    ALLOCATE(Neighc_pme_(i_cut) % Nei(nprocs))
    ALLOCATE(mask(ncx,ncy,ncz))

    Nei=>Neighc_(i_cut) % Nei
    CALL CellZero
    CALL AllCells(nei)
    Nei=>Neighc_pme_(i_cut) % Nei
    CALL CellZero_pme
    CALL AllCells(nei)
    CALL PlanesBoxes
  CONTAINS
    INCLUDE 'NearBox__CellZero.f90'
    INCLUDE 'NearBox__CellZero_pme.f90'
    SUBROUTINE AllCells(nei)
      TYPE(NearBox__Neigh) :: Nei(:)
      INTEGER :: i,j,k,mx,my,mz,m0x,m0y,m0z,m1x,m1y,m1z
      REAL(8) :: x1,y1,z1,sqcut,rmin,xc,yc,zc
      INTEGER :: vect0(3),count1,numcell,n0
      INTEGER, POINTER :: vect(:)

      n0=SIZE(Ind_xyz)
      mask=.FALSE.
      DO i=1,npx
         x1=(i-1)*ddx/dx
         DO j=1,npy
            y1=(j-1)*ddy/dy
            DO k=1,npz
               z1=(k-1)*ddz/dz
               numcell=(i-1)*npy*npz+(j-1)*npz+k

               mx=INT(x1)+(SIGN(1.D0,x1-INT(x1))-1.)/2
               my=INT(y1)+(SIGN(1.D0,y1-INT(y1))-1.)/2
               mz=INT(z1)+(sign(1.d0,z1-int(z1))-1.)/2
               m1x=MOD(MOD(mx,ncx)+ncx,ncx)
               m1y=MOD(MOD(my,ncy)+ncy,ncy)
               m1z=MOD(MOD(mz,ncz)+ncz,ncz)

               ALLOCATE(Nei(numcell) % c (3,n0))
               DO m=1,SIZE(Ind_xyz)
                  iv=Ind_xyz(m) % i
                  jv=Ind_xyz(m) % j
                  kv=Ind_xyz(m) % k
                  nx=mod(mod(m1x+iv,ncx)+ncx,ncx)+1
                  ny=mod(mod(m1y+jv,ncy)+ncy,ncy)+1
                  nz=mod(mod(m1z+kv,ncz)+ncz,ncz)+1
                  vect0(1)=nx; vect0(2)=ny; vect0(3)=nz
                  Nei(numcell) % c(:,m) = vect0
               END DO
            END DO
         END DO
      END DO
!!$      IF(PI_Node_Cart == 0) THEN
!!$         DO i=1,npx
!!$            DO j=1,npy
!!$               DO k=1,npz
!!$                  numcell=(i-1)*npy*npz+(j-1)*npz+k
!!$                  n0=SIZE(Nei(numcell) % c,2)
!!$                  WRITE(100+numcell-1,'(''ijk '',3i8)') i,j,k
!!$                  DO n=1,n0
!!$                     WRITE(100+numcell-1,'(3i8)') Nei(numcell) % c(:,n)
!!$                  END DO
!!$               END DO
!!$            END DO
!!$         END DO
!!$      END IF
!!$      STOP
    END SUBROUTINE AllCells
    SUBROUTINE PlanesBoxes
      INTEGER :: i,j,k,numcell
      REAL(8) :: x1,x2,y1,y2,z1,z2,v1(3),v2(3),v3(3),r1(3),r2(3)&
           &,r3(3),vxy1(4),vxz1(4),vyz1(4),vxy2(4),vxz2(4),vyz2(4)
      
      ALLOCATE(Planes(nprocs))
      DO i=1,npx
         x1=(i-1)*ddx
         x2=i*ddx
         DO j=1,npy
            y1=(j-1)*ddy
            y2=j*ddy
            DO k=1,npz
               z1=(k-1)*ddz
               z2=k*ddz
               numcell=(i-1)*npy*npz+(j-1)*npz+k
!!$
!!$--- Get xy plane low side
!!$
               v1=(/x1,y1,z1/)
               v2=(/x2,y1,z1/)
               v3=(/x1,y2,z1/)
               r1=Convert(v1)
               r2=Convert(v2)
               r3=Convert(v3)
               vxy1=Equation_Plane(r1,r2,r3)

!!$
!!$--- Get xz plane low side
!!$
               v1=(/x1,y1,z1/)
               v2=(/x2,y1,z1/)
               v3=(/x1,y1,z2/)
               r1=Convert(v1)
               r2=Convert(v2)
               r3=Convert(v3)
               vxz1=Equation_Plane(r1,r2,r3)

!!$
!!$--- Get yz plane low side
!!$
               v1=(/x1,y1,z1/)
               v2=(/x1,y2,z1/)
               v3=(/x1,y1,z2/)
               r1=Convert(v1)
               r2=Convert(v2)
               r3=Convert(v3)
               vyz1=Equation_Plane(r1,r2,r3)

!!$
!!$--- Get xy plane upper side
!!$
               v1=(/x1,y1,z2/)
               v2=(/x2,y1,z2/)
               v3=(/x1,y2,z2/)
               r1=Convert(v1)
               r2=Convert(v2)
               r3=Convert(v3)
               vxy2=Equation_Plane(r1,r2,r3)

!!$
!!$--- Get xz plane upper side
!!$
               v1=(/x1,y2,z1/)
               v2=(/x2,y2,z1/)
               v3=(/x1,y2,z2/)
               r1=Convert(v1)
               r2=Convert(v2)
               r3=Convert(v3)
               vxz2=Equation_Plane(r1,r2,r3)

!!$
!!$--- Get yz plane upper side
!!$
               v1=(/x2,y1,z1/)
               v2=(/x2,y2,z1/)
               v3=(/x2,y1,z2/)
               r1=Convert(v1)
               r2=Convert(v2)
               r3=Convert(v3)
               vyz2=Equation_Plane(r1,r2,r3)


               Planes(numcell) % vxy0 = vxy1
               Planes(numcell) % vxz0 = vxz1
               Planes(numcell) % vyz0 = vyz1

               Planes(numcell) % vxy1 = vxy2
               Planes(numcell) % vxz1 = vxz2
               Planes(numcell) % vyz1 = vyz2
            END DO
         END DO
      END DO
    END SUBROUTINE PlanesBoxes
    SUBROUTINE Init
      INTEGER :: i,j,k,mx,my,mz,m0x,m0y,m0z,m1x,m1y,m1z
      REAL(8) :: x1,y1,z1,sqcut
      INTEGER :: vect0(3),n0,numcell,nvalues(3)
      INTEGER, POINTER :: vect(:)
      TYPE :: counter3D
         INTEGER :: n0
      END type counter3D
      TYPE(Counter3D), ALLOCATABLE :: dummy(:,:,:)
      LOGICAL :: okkk

      okkk=.FALSE.
!!$       nx=NINT(a/rcut0)
!!$       ny=NINT(b/rcut0)
!!$       nz=NINT(c/rcut0)
      nx=24
      ny=24
      nz=24
      iter_inst=0
      DO WHILE(.NOT. okkk)
         iter_inst=iter_inst+1
         IF(nx < npx .OR. ny < npy .OR. nz < npz) THEN
            IF(nx < npx) THEN
               WRITE(labs0,'(i3,i3)') nx,npx
               errmsg_f='Number of neighbor list cells along X is smaller than&
                    & the number of processors: There are '//labs0(1:3)//' neighbor &
                    &list cells and '//labs0(4:6)//' No. processors on X'
            ELSE IF(ny < npy) THEN
               WRITE(labs0,'(i3,i3)') ny,npy
               errmsg_f='Number of neighbor list cells along Y is smaller than&
                    & the number of processors: There are '//labs0(1:3)//' neighbor &
                    &list cells and '//labs0(4:6)//' No. processors on Y'
            ELSE IF(nz < npz) THEN
               WRITE(labs0,'(i3,i3)') nz,npz
               errmsg_f='Number of neighbor list cells along Z is smaller than&
                    & the number of processors: There are '//labs0(1:3)//' neighbor &
                    &list cells and '//labs0(4:6)//' No. processors on Z'
            END IF
            CALL Add_Errors(-1,errmsg_f)
            out=.FALSE.
            RETURN
         END IF
         
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
         
         nvalues=Check(i_cut,rcut,ncx,ncy,ncz)
         IF(nvalues(1) == 0 .AND. nvalues(2) == 0 .AND. nvalues(3) == 0) THEN
            okkk=.TRUE.
         ELSE
            nx=nx+nvalues(1)
            ny=ny+nvalues(2)
            nz=nz+nvalues(3)
         END IF
!!$
!!$--- Abort after iter_Max unsuccessfull iterations
!!$
         IF(iter_inst > iter_Max) okkk=.TRUE.
      END DO

      dx=2.d0/DBLE(ncx)
      dy=2.d0/DBLE(ncy)
      dz=2.d0/DBLE(ncz)
      ddx=2.d0/DBLE(npx)
      ddy=2.d0/DBLE(npy)
      ddz=2.d0/DBLE(npz)


      ALLOCATE(dummy(npx,npy,npz))
      ALLOCATE(Ind_Large(npx,npy,npz))
      ALLOCATE(Ind_Small(ncx,ncy,ncz))
      dummy(:,:,:) % n0=0

      DO i=1,ncx
         x1=(i-1)*dx/ddx
         DO j=1,ncy
            y1=(j-1)*dy/ddy
            DO k=1,ncz
               z1=(k-1)*dz/ddz
               mx=INT(x1)+(SIGN(1.0D0,x1-INT(x1))-1.0D0)/2
               my=INT(y1)+(SIGN(1.0D0,y1-INT(y1))-1.0D0)/2
               mz=INT(z1)+(sign(1.0D0,z1-int(z1))-1.0D0)/2
               mx=MOD(MOD(mx,npx)+npx,npx)+1
               my=MOD(MOD(my,npy)+npy,npy)+1
               mz=MOD(MOD(mz,npz)+npz,npz)+1
               dummy(mx,my,mz) % n0=dummy(mx,my,mz) % n0+1
            END DO
         END DO
      END DO
      DO mx=1,npx
         DO my=1,npy
            DO mz=1,npz
               n0=dummy(mx,my,mz) % n0
               ALLOCATE(Ind_Large(mx,my,mz) % pt(n0))
            END DO
         END DO
      END DO
         
      dummy(:,:,:) % n0=0
      DO i=1,ncx
         x1=(i-1)*dx/ddx
         DO j=1,ncy
            y1=(j-1)*dy/ddy
            DO k=1,ncz
               z1=(k-1)*dz/ddz
               mx=INT(x1)+(SIGN(1.0D0,x1-INT(x1))-1.0D0)/2
               my=INT(y1)+(SIGN(1.0D0,y1-INT(y1))-1.0D0)/2
               mz=INT(z1)+(sign(1.0D0,z1-int(z1))-1.0D0)/2
               mx=MOD(MOD(mx,npx)+npx,npx)+1
               my=MOD(MOD(my,npy)+npy,npy)+1
               mz=MOD(MOD(mz,npz)+npz,npz)+1

               numcell=(mx-1)*npy*npz+(my-1)*npz+mz

               dummy(mx,my,mz) % n0=dummy(mx,my,mz) % n0+1
               n0=dummy(mx,my,mz) % n0
               Ind_Large(mx,my,mz) % pt(n0) % nx = i
               Ind_Large(mx,my,mz) % pt(n0) % ny = j
               Ind_Large(mx,my,mz) % pt(n0) % nz = k
               Ind_Small(i,j,k) = numcell
            END DO
         END DO
      END DO
    END SUBROUTINE Init
  END FUNCTION NearBox_
  FUNCTION Check(i_n,rcut,nx,ny,nz) RESULT(out)
    INTEGER :: out(3)
    INTEGER :: nx,ny,nz,i_n
    REAL(8) :: rcut

    INTEGER :: mfx,mfy,mfz,ncx,ncy,ncz
    REAL(8) :: sqcut,dx,dy,dz,rmin
    INTEGER :: imax,jmax,kmax,i,j,k,istart,jstart,kstart,warnx,warny, warnz&
         &,nxmax, nymax,nzmax,kend

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
  END FUNCTION Check
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
  FUNCTION Convert(v1) RESULT(out)
    REAL(8) :: v1(3),out(3)
    REAL(8) :: r1(3)
!!$    r1(1)=co(1,1)*v1(1)+co(1,2)*v1(2)+co(1,3)*v1(3)
!!$    r1(2)=co(2,1)*v1(1)+co(2,2)*v1(2)+co(2,3)*v1(3)
!!$    r1(3)=co(3,1)*v1(1)+co(3,2)*v1(2)+co(3,3)*v1(3)
    r1(1)=v1(1)
    r1(2)=v1(2)
    r1(3)=v1(3)
    out=r1
  END FUNCTION Convert
  SUBROUTINE NearBox__Param(nx,ny,nz)
    INTEGER :: nx,ny,nz
    nx=ncx; ny=ncy; nz=ncz
  END SUBROUTINE NearBox__Param
END MODULE NearBox
