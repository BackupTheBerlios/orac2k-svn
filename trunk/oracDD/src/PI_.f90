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
MODULE PI_
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Apr  5 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program ORAC ----*

#ifdef HAVE_MPI
  USE MPI
#endif
  USE Print_Defs
  IMPLICIT none
  PRIVATE
  PUBLIC :: PI__, PI__Finalize, PI_Comm,PI_Node,PI_nprocs,PI_npx&
       &,PI_npy,PI_npz, PI_Comm_Cart,PI__Node, PI__Ranks,&
       & PI__Setup_Cart, PI_Node_Cart,&
       & PI__Setup_SndRcv,ierr,status,PI_Comm_Z,PI_Node_Z,PI_Nprocs_Z&
       &, PI__Ranks_Z, PI_Node_FFTW, PI_Comm_FFTW, PI__Ranks_FFTW

  INTEGER, SAVE :: PI_Comm,PI_node=0,PI_nprocs,PI_npx=1,PI_npy=1&
       &,PI_npz=1, PI_Comm_Cart, PI_Node_Cart,PI_Comm_Z,PI_Node_Z&
       &,PI_Nprocs_Z,PI_group, PI_Node_FFTW, PI_Comm_FFTW,&
       & PI_Group_Cart,PI_Group_Z,PI_Group_FFTW
  TYPE :: PI_shift
     INTEGER :: p_s,p_d,m_s,m_d
  END type PI_shift
  TYPE(PI_Shift), SAVE :: PI_Nei_x,PI_Nei_y,PI_Nei_z
  TYPE :: PI__Node
     INTEGER :: n, nx, ny, nz
  END TYPE PI__Node
  TYPE(PI__Node), ALLOCATABLE, SAVE :: PI__Ranks(:),PI__Ranks_Z(:),PI__Ranks_FFTW(:)
  TYPE(PI_Shift), SAVE :: PI_Nei (3)


  INTEGER :: ierr
  INTEGER :: Status(MPI_Status_Size)
CONTAINS
  SUBROUTINE PI__
    INTEGER :: ierr
#ifdef HAVE_MPI
    CALL MPI_Init( ierr )
    PI_comm=MPI_COMM_WORLD
    call MPI_COMM_RANK( PI_comm, PI_node, ierr )
    call MPI_COMM_SIZE( PI_comm, PI_nprocs, ierr )
    IF(PI_node /=0 ) THEN
       OPEN(unit=kprint,file="/dev/null")
    END IF
#endif
  END SUBROUTINE PI__
  SUBROUTINE PI__Nodes
#ifdef HAVE_MPI
    
#endif
  END SUBROUTINE PI__Nodes
  SUBROUTINE PI__Finalize
    INTEGER :: ierr

#ifdef HAVE_MPI
    CALL MPI_Finalize( ierr )
#else
    STOP
#endif

  END SUBROUTINE PI__Finalize
  SUBROUTINE PI__Setup_Cart

    INTEGER :: dims(3),ierr,coord(3),coord2d(2)
    INTEGER :: m,n,mx,my,mz,npx,npy,npz,color,key
    LOGICAL :: periods(3)=(/.TRUE., .TRUE., .TRUE. /)
    LOGICAL :: reorder=.TRUE.
    LOGICAL :: remains(3)=(/.TRUE., .TRUE., .FALSE. /)
    INTEGER, ALLOCATABLE :: coords(:,:),ranks(:),ranks2(:)
    INTEGER, ALLOCATABLE :: coord2ds(:,:)
    INTEGER, SAVE :: PI_Comm_1,PI_Node_1,PI_Group_1
#ifdef HAVE_MPI

    dims(1)=PI_npx
    dims(2)=PI_npy
    dims(3)=PI_npz

    ALLOCATE(coords(3,PI_Nprocs))
    ALLOCATE(PI__Ranks(PI_Nprocs))
    ALLOCATE(ranks(PI_Nprocs))

    CALL MPI_CART_CREATE(PI_Comm, 3, dims, periods, reorder, PI_Comm_Cart, ierr)
    CALL MPI_COMM_RANK(PI_Comm_Cart, PI_Node_Cart, ierr)
    CALL MPI_CART_COORDS(PI_Comm_Cart, PI_Node_Cart, 3, coord, ierr)

    CALL MPI_ALLGATHER(coord,3,MPI_INTEGER4,coords,3,MPI_INTEGER4&
         &,PI_Comm_Cart, ierr)

    npx=PI_npx
    npy=PI_npy
    npz=PI_npz

    DO m=1,PI_nprocs
       n=m-1
       mx=coords(1,m)+1
       my=coords(2,m)+1
       mz=coords(3,m)+1
       PI__Ranks (m) % n = (mx-1)*npy*npz+(my-1)*npz+mz
       PI__Ranks (m) % nx = mx-1
       PI__Ranks (m) % ny = my-1
       PI__Ranks (m) % nz = mz-1
       ranks((mz-1)*npx*npy+(my-1)*npx+mx)=n
    END DO

    CALL MPI_COMM_GROUP(PI_Comm_Cart, PI_Group_Cart, ierr)
    CALL MPI_GROUP_INCL(PI_Group_Cart,PI_nprocs,ranks,PI_Group_FFTW,&
         & ierr)
    CALL MPI_COMM_CREATE(PI_Comm_Cart,PI_Group_FFTW,PI_Comm_FFTW, ierr)
    CALL MPI_COMM_RANK(PI_Comm_FFTW, PI_Node_FFTW, ierr)
    ALLOCATE(PI__Ranks_FFTW(PI_Nprocs))
    DO m=1,PI_nprocs
       n=m-1
       PI__Ranks_FFTW (m) % n = PI__Ranks(ranks(m)+1) % n
       PI__Ranks_FFTW (m) % nx = PI__Ranks(ranks(m)+1) % nx
       PI__Ranks_FFTW (m) % ny = PI__Ranks(ranks(m)+1) % ny
       PI__Ranks_FFTW (m) % nz = PI__Ranks(ranks(m)+1) % nz
    END DO
    mx=PI__Ranks_FFTW (PI_Node_FFTW+1) % nx+1
    my=PI__Ranks_FFTW (PI_Node_FFTW+1) % ny+1
    key=(my-1)*npx+mx
    color=PI__Ranks_FFTW (PI_Node_FFTW+1) % nz

    CALL MPI_COMM_SPLIT(PI_Comm_FFTW,color,key,PI_Comm_Z,ierr)
    CALL MPI_COMM_RANK(PI_Comm_Z, PI_Node_Z, ierr)
    CALL MPI_COMM_SIZE( PI_comm_Z, PI_nprocs_Z, ierr )

    coord2d(1)=mx-1
    coord2d(2)=my-1
    ALLOCATE(coord2ds(2,PI_Nprocs_Z))
    CALL MPI_ALLGATHER(coord2d,2,MPI_INTEGER4,coord2ds,2,MPI_INTEGER4&
         &,PI_Comm_Z, ierr)

    ALLOCATE(PI__Ranks_Z(PI_Nprocs_Z))
    DO m=1,PI_nprocs_Z
       n=m-1
       mx=coord2ds(1,m)+1
       my=coord2ds(2,m)+1
       PI__Ranks_Z (m) % n = (my-1)*npx+mx
       PI__Ranks_Z (m) % nx = mx-1
       PI__Ranks_Z (m) % ny = my-1
    END DO
#endif
  END SUBROUTINE PI__Setup_Cart
  SUBROUTINE PI__Setup_SndRcv
    INTEGER :: ierr

#ifdef HAVE_MPI
    CALL MPI_CART_SHIFT(PI_Comm_Cart,0,1,PI_Nei(1) % p_s,PI_Nei(1) %&
         & p_d,ierr) 
    CALL MPI_CART_SHIFT(PI_Comm_Cart,0,-1,PI_Nei(1) % m_s,PI_Nei(1) %&
         & m_d,ierr)
    CALL MPI_CART_SHIFT(PI_Comm_Cart,1,1,PI_Nei(2) % p_s,PI_Nei(2) %&
         & p_d,ierr) 
    CALL MPI_CART_SHIFT(PI_Comm_Cart,1,-1,PI_Nei(2) % m_s,PI_Nei(2) %&
         & m_d,ierr)
    CALL MPI_CART_SHIFT(PI_Comm_Cart,2,1,PI_Nei(3) % p_s,PI_Nei(3) %&
         & p_d,ierr) 
    CALL MPI_CART_SHIFT(PI_Comm_Cart,2,-1,PI_Nei(3) % m_s,PI_Nei(3) %&
         & m_d,ierr)
#endif

  END SUBROUTINE PI__Setup_SndRcv
END MODULE PI_
