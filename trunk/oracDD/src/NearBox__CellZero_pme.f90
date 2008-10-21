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
SUBROUTINE CellZero_pme
  INTEGER :: i,j,k,mx,my,mz,m0x,m0y,m0z,m1x,m1y,m1z,nind,np,n,nx&
       &,ny,nz,mxlast,mylast,mzlast,count0
  REAL(8) :: x1,y1,z1,sqcut,x0,y0,z0,rmin
  INTEGER :: vect0(3)
  INTEGER, POINTER :: vect(:)
  
  sqcut=(rcut)**2
  
!!$
!!$--- Edges of cell zero in unit of small cells
!!$
  
  x0=0.0D0 ; y0=0.0D0; z0=0.0D0
  x1=ddx/dx ; y1=ddy/dy ; z1=ddz/dz
  
  m0x=0
  m0y=0
  m0z=0
  
  mx=INT(x1)+(SIGN(1.D0,x1-INT(x1))-1.)/2
  my=INT(y1)+(SIGN(1.D0,y1-INT(y1))-1.)/2
  mz=INT(z1)+(sign(1.d0,z1-int(z1))-1.)/2
  m1x=MOD(MOD(mx,ncx)+ncx,ncx)
  m1y=MOD(MOD(my,ncy)+ncy,ncy)
  m1z=MOD(MOD(mz,ncz)+ncz,ncz)
  
  mxlast=m1x
  mylast=m1y
  mzlast=m1z
  IF(mx == ncx) mxlast=ncx-1
  IF(my == ncy) mylast=ncy-1
  IF(mz == ncz) mzlast=ncz-1
  
  mask=.FALSE.
  DO n=1,SIZE(Ind_Large(1,1,1) % pt)
     i=Ind_Large(1,1,1) % pt(n) % nx
     j=Ind_Large(1,1,1) % pt(n) % ny
     k=Ind_Large(1,1,1) % pt(n) % nz
     mask(i,j,k)=.TRUE.
  END DO

!!$      
!!$--- Cells x > x1
!!$
  
  IF(.NOT. Node_()) STOP
  count0=0
  IF(m1x /= 0) THEN
     mx=m1x
     DO my=m0y,mylast
        DO mz=m0z,mzlast         
           
           DO i=0,ncx-1
              DO j=1-ncy,ncy-1
                 DO k=1-ncz,ncz-1
                    rmin=dist_ijk(i,j,k,dx,dy,dz)
                    count0=count0+1
                    IF(rmin < sqcut) THEN
                       nx=mod(mod(mx+i,ncx)+ncx,ncx)
                       ny=mod(mod(my+j,ncy)+ncy,ncy)
                       nz=mod(mod(mz+k,ncz)+ncz,ncz)
                       IF(.NOT. mask(nx+1,ny+1,nz+1)) THEN
                          vect0=(/nx, ny, nz/) 
                          CALL Node__Push(vect0)
                          mask(nx+1,ny+1,nz+1)=.TRUE.
                       END IF
                    END IF
                 END DO
              END DO
           END DO
           
        END DO
     END DO
     mx=m0x
     DO my=m0y,mylast
        DO mz=m0z,mzlast         
           
           DO i=1-ncx,-1
              DO j=1-ncy,ncy-1
                 DO k=1-ncz,ncz-1
                    rmin=dist_ijk(i,j,k,dx,dy,dz)
                    count0=count0+1
                    IF(rmin < sqcut) THEN
                       nx=mod(mod(mx+i,ncx)+ncx,ncx)
                       ny=mod(mod(my+j,ncy)+ncy,ncy)
                       nz=mod(mod(mz+k,ncz)+ncz,ncz)
                       IF(.NOT. mask(nx+1,ny+1,nz+1)) THEN
                          vect0=(/nx, ny, nz/) 
                          CALL Node__Push(vect0)
                          mask(nx+1,ny+1,nz+1)=.TRUE.
                       END IF
                    END IF
                 END DO
              END DO
           END DO
           
        END DO
     END DO
  END IF

!!$      
!!$--- Cells y < y0
!!$

  IF(m1y /= 0) THEN
     my=m0y
     DO mx=m0x,mxlast
        DO mz=m0z,mzlast
           DO i=m0x,mxlast
              DO j=1-ncy,-1
                 DO k=1-ncz,ncz-1
                    rmin=dist_ijk(i,j,k,dx,dy,dz)
                    count0=count0+1
                    IF(rmin < sqcut) THEN
                       nx=mod(mod(mx+i,ncx)+ncx,ncx)
                       ny=mod(mod(my+j,ncy)+ncy,ncy)
                       nz=mod(mod(mz+k,ncz)+ncz,ncz)
                       IF(.NOT. mask(nx+1,ny+1,nz+1)) THEN
                          vect0=(/nx, ny, nz/) 
                          CALL Node__Push(vect0)
                          mask(nx+1,ny+1,nz+1)=.TRUE.
                       END IF
                    END IF
                 END DO
              END DO
           END DO
           
        END DO
     END DO
     my=m1y
     DO mx=m0x,mxlast
        DO mz=m0z,mzlast
           
           DO i=m0x,mxlast
              DO j=0,ncy-1
                 DO k=1-ncz,ncz-1
                    count0=count0+1
                    rmin=dist_ijk(i,j,k,dx,dy,dz)
                    IF(rmin < sqcut) THEN
                       nx=mod(mod(mx+i,ncx)+ncx,ncx)
                       ny=mod(mod(my+j,ncy)+ncy,ncy)
                       nz=mod(mod(mz+k,ncz)+ncz,ncz)
                       IF(.NOT. mask(nx+1,ny+1,nz+1)) THEN
                          vect0=(/nx, ny, nz/) 
                          CALL Node__Push(vect0)
                          mask(nx+1,ny+1,nz+1)=.TRUE.
                       END IF
                    END IF
                 END DO
              END DO
           END DO
           
           
        END DO
     END DO
  END IF
!!$      
!!$--- Cells z < z0
!!$

  IF(m1z /= 0) THEN
     mz=m0z
     DO mx=m0x,mxlast
        DO my=m0y,mylast
           
           DO i=m0x,mxlast
              DO j=m0y,mylast
                 DO k=1-ncz,-1
                    rmin=dist_ijk(i,j,k,dx,dy,dz)
                    count0=count0+1
                    IF(rmin < sqcut) THEN
                       nx=mod(mod(mx+i,ncx)+ncx,ncx)
                       ny=mod(mod(my+j,ncy)+ncy,ncy)
                       nz=mod(mod(mz+k,ncz)+ncz,ncz)
                       IF(.NOT. mask(nx+1,ny+1,nz+1)) THEN
                          vect0=(/nx, ny, nz/) 
                          CALL Node__Push(vect0)
                          mask(nx+1,ny+1,nz+1)=.TRUE.
                       END IF
                    END IF
                 END DO
              END DO
           END DO
        END DO
     END DO
     mz=m1z
     DO mx=m0x,mxlast
        DO my=m0y,mylast
           
           DO i=m0x,mxlast
              DO j=m0y,mylast
                 DO k=0,ncz-1
                    rmin=dist_ijk(i,j,k,dx,dy,dz)
                    count0=count0+1
                    IF(rmin < sqcut) THEN
                       nx=mod(mod(mx+i,ncx)+ncx,ncx)
                       ny=mod(mod(my+j,ncy)+ncy,ncy)
                       nz=mod(mod(mz+k,ncz)+ncz,ncz)
                       IF(.NOT. mask(nx+1,ny+1,nz+1)) THEN
                          vect0=(/nx, ny, nz/) 
                          CALL Node__Push(vect0)
                          mask(nx+1,ny+1,nz+1)=.TRUE.
                       END IF
                    END IF
                 END DO
              END DO
           END DO
        END DO
     END DO
  END IF
  
  nind=Node__Size()
  IF(ALLOCATED(ind_xyz)) DEALLOCATE(ind_xyz)
  ALLOCATE(ind_xyz(nind))
  nind=0
  DO WHILE(Node__Pop(vect))
     nind=nind+1
     ind_xyz(nind) % i=vect(1)
     ind_xyz(nind) % j=vect(2)
     ind_xyz(nind) % k=vect(3)
  END DO
END SUBROUTINE CellZero_pme
