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
MODULE PI_Communicate
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Jan 26 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*
  
#ifdef HAVE_MPI
  USE mpi
#endif
  USE UNITS
  USE PI_
  IMPLICIT none
  PRIVATE
  PUBLIC 
  REAL(8), ALLOCATABLE, SAVE :: Buffer(:)
CONTAINS
  SUBROUTINE PI__Communicate(i_n)
    INTEGER :: i_n
    
    INTEGER :: n_p

    n_pi=PI__Ranks(PI_Node_Cart+1) % n 
    CALL Gather_Send

    DO n=1,PI_Nprocs
       nn=PI__Ranks(n+1) % n 
       IF(SIZE(Pe_nb(i_n) % tbl (n_pi, nn) % exc) /= 0) THEN
          mm=SIZE(Pe_nb(i_n) % tbl (n_pi, nn) % exc)
          DO m=1,mm
             (/ox, oy, oz/) = Pe_nb(i_n) % tbl (n_pi, nn) % exc(:,m)
             CALL Gather_Send(ox,oy,oz)
          END DO
       END IF
    END DO
  CONTAINS
    SUBROUTINE Gather_Send(ox,oy,oz)
      INTEGER, OPTIONAL :: ox,oy,oz
      INTEGER, SAVE :: ntotal
      
      IF(.NOT. PRESENT(ox)) THEN
         ntotal=0
      END IF
      

    END SUBROUTINE Gather_Send
  END SUBROUTINE PI__Communicate
  
  SUBROUTINE PI__Main
    INTEGER :: n
    n=SIZE(Atoms)*3
    CALL PI__Buffer_A(n)
  END SUBROUTINE PI__Main
  SUBROUTINE PI__Buffer_A(n)
    INTEGER :: n
    INTEGER :: ierr

    ierr=0
#ifdef HAVE_MPI
    ALLOCATE(Buffer(n))
    CALL MPI_BUFFER_ATTACH(Buffer, n*8, ierr)
#endif
  END SUBROUTINE PI__Buffer_A
  SUBROUTINE PI__Buffer_D
    INTEGER :: n
    INTEGER :: ierr

    ierr=0
#ifdef HAVE_MPI
    
    n=SIZE(Buffer)*8
    CALL MPI_BUFFER_DETACH(Buffer, n, ierr)
    DEALLOCATE(Buffer)
#endif
  END SUBROUTINE PI__Buffer_D
END MODULE PI_Communicate
