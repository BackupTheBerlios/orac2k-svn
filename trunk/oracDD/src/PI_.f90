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

  USE Print_Defs
  IMPLICIT none
#ifdef HAVE_MPI
  include 'mpif.h'
#endif
  PRIVATE
  PUBLIC :: PI__, PI__Finalize, PI_Comm,PI_Node,PI_nprocs,PI_npx&
       &,PI_npy,PI_npz, PI_Comm_Cart,PI__Node, PI__Ranks,&
       & PI__Get_Ranks, PI_Node_Cart
  INTEGER, SAVE :: PI_Comm,PI_node=0,PI_nprocs,PI_npx=1,PI_npy=1&
       &,PI_npz=1, PI_Comm_Cart, PI_Node_Cart
  TYPE :: PI__Node
     INTEGER :: n, nx, ny, nz
  END TYPE PI__Node
  TYPE(PI__Node), ALLOCATABLE :: PI__Ranks(:)
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
#endif

  END SUBROUTINE PI__Finalize
  SUBROUTINE PI__Get_Ranks

    INTEGER :: dims(3),ierr,coord(3)
    INTEGER :: m,n,mx,my,mz,npx,npy,npz
    LOGICAL :: periods(3)=(/.TRUE., .TRUE., .TRUE. /)
    LOGICAL :: reorder=.TRUE.
    INTEGER, POINTER :: coords(:,:)
    
#ifdef HAVE_MPI

    dims(1)=PI_npx
    dims(2)=PI_npy
    dims(3)=PI_npz


    ALLOCATE(coords(3,PI_Nprocs))
    ALLOCATE(PI__Ranks(PI_Nprocs))

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
    END DO
#endif
  END SUBROUTINE PI__Get_Ranks
END MODULE PI_
