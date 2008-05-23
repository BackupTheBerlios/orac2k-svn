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
MODULE PI_Neighbors

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
  USE PI_
  IMPLICIT none
  PRIVATE
  PUBLIC PI_Neighbors_,Head_xyz,Chain_xyz,Chain

  TYPE :: Chain
     INTEGER :: i,j,k
     INTEGER :: p
  END TYPE Chain

  TYPE(Chain), ALLOCATABLE, SAVE :: Chain_xyz(:)
  INTEGER, ALLOCATABLE, SAVE :: Head_xyz(:)
  INTEGER, SAVE :: ncx,ncy,ncz
CONTAINS
!!$
!!$--- Constructor for PI_Neighbors_
!!$
  FUNCTION PI_Neighbors_(x,y,z) RESULT(out)
    LOGICAL :: out
    REAL(8) :: x(:),y(:),z(:)
    REAL(8) :: x1,y1,z1,dx,dy,dz
    INTEGER :: n,nx,ny,nz,natp,numcell,l
    INTEGER, SAVE :: calls = 0 

    IF(calls == 0) THEN
       ncx=PI_npx
       ncy=PI_npy
       ncz=PI_npz
       ALLOCATE(Head_xyz(ncx*ncy*ncz))
       natp=SIZE(x)
       ALLOCATE(Chain_xyz(natp))
    END IF

    out=.TRUE.

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
  END FUNCTION PI_Neighbors_

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE PI_Neighbors
