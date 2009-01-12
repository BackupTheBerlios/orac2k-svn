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
MODULE PI_Collectives
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Dec 11 2008 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*
  

#ifdef HAVE_MPI
  USE mpi
#endif
  USE PI_
  USE IndBox, ONLY: IndBox_a_p
  USE Print_Defs
  USE CONSTANTS, ONLY: max_data
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f, errmsg_w
  IMPLICIT none
  PRIVATE
  PUBLIC PI_Gather_, PI_AllGather_,PI_Write_, PI_ErrSignal_
CONTAINS
  SUBROUTINE PI_Gather_(xc,yc,zc)
    REAL(8) :: xc(:),yc(:),zc(:)
    REAL(8), ALLOCATABLE :: x(:),y(:),z(:)
    INTEGER, ALLOCATABLE :: iBuff(:),locals(:),displ(:)
    INTEGER :: natom,nlocal,n,nn,pt

#ifdef HAVE_MPI
    natom=SIZE(xc)

    ALLOCATE(x(natom),y(natom),z(natom),iBuff(natom))
    ALLOCATE(locals(PI_Nprocs),displ(PI_Nprocs))

    IF(PI_Nprocs > 1) THEN
       nlocal=SIZE(IndBox_a_p)
       CALL MPI_ALLGATHER(nlocal,1,MPI_INTEGER4,locals,1&
            &,MPI_INTEGER4,PI_Comm_Cart,ierr)
       displ(1)=0
       DO n=2,PI_Nprocs
          displ(n)=displ(n-1)+locals(n-1)
       END DO
       
       pt=displ(PI_Node_Cart+1)
       DO nn=1,SIZE(IndBox_a_p)
          n=IndBox_a_p(nn)
          iBuff(pt+nn)=n
          x(pt+nn)=xc(n)
          y(pt+nn)=yc(n)
          z(pt+nn)=zc(n)
       END DO

       CALL MPI_GATHERV(x(pt+1),nlocal,MPI_REAL8,x,locals&
            &,displ,MPI_REAL8,0,PI_Comm_Cart,ierr)
       CALL MPI_GATHERV(y(pt+1),nlocal,MPI_REAL8,y,locals&
            &,displ,MPI_REAL8,0,PI_Comm_Cart,ierr)
       CALL MPI_GATHERV(z(pt+1),nlocal,MPI_REAL8,z,locals&
            &,displ,MPI_REAL8,0,PI_Comm_Cart,ierr)
       CALL MPI_GATHERV(iBuff(pt+1),nlocal,MPI_INTEGER4,iBuff,locals&
            &,displ,MPI_INTEGER4,0,PI_Comm_Cart,ierr)       

       IF(PI_Node_Cart == 0) THEN
          xc(iBuff(:))=x
          yc(iBuff(:))=y
          zc(iBuff(:))=z
       END IF
    END IF
#else
    RETURN
#endif    
  END SUBROUTINE PI_Gather_

  SUBROUTINE PI_AllGather_(xc,yc,zc)
    REAL(8) :: xc(:),yc(:),zc(:)
    REAL(8), ALLOCATABLE :: x(:),y(:),z(:)
    INTEGER, ALLOCATABLE :: iBuff(:),locals(:),displ(:)
    INTEGER :: natom,nlocal,n,nn,pt

#ifdef HAVE_MPI
    natom=SIZE(xc)

    ALLOCATE(x(natom),y(natom),z(natom),iBuff(natom))
    ALLOCATE(locals(PI_Nprocs),displ(PI_Nprocs))

    IF(PI_Nprocs > 1) THEN
       nlocal=SIZE(IndBox_a_p)
       CALL MPI_ALLGATHER(nlocal,1,MPI_INTEGER4,locals,1&
            &,MPI_INTEGER4,PI_Comm_Cart,ierr)
       displ(1)=0
       DO n=2,PI_Nprocs
          displ(n)=displ(n-1)+locals(n-1)
       END DO
       
       pt=displ(PI_Node_Cart+1)
       DO nn=1,SIZE(IndBox_a_p)
          n=IndBox_a_p(nn)
          iBuff(pt+nn)=n
          x(pt+nn)=xc(n)
          y(pt+nn)=yc(n)
          z(pt+nn)=zc(n)
       END DO

       CALL MPI_ALLGATHERV(x(pt+1),nlocal,MPI_REAL8,x,locals&
            &,displ,MPI_REAL8,PI_Comm_Cart,ierr)
       CALL MPI_ALLGATHERV(y(pt+1),nlocal,MPI_REAL8,y,locals&
            &,displ,MPI_REAL8,PI_Comm_Cart,ierr)
       CALL MPI_ALLGATHERV(z(pt+1),nlocal,MPI_REAL8,z,locals&
            &,displ,MPI_REAL8,PI_Comm_Cart,ierr)
       CALL MPI_ALLGATHERV(iBuff(pt+1),nlocal,MPI_INTEGER4,iBuff,locals&
            &,displ,MPI_INTEGER4,PI_Comm_Cart,ierr)

       xc(iBuff(:))=x
       yc(iBuff(:))=y
       zc(iBuff(:))=z
    END IF
#else
    RETURN
#endif    
  END SUBROUTINE PI_AllGather_
  SUBROUTINE PI_Write_(unit,xc,yc,zc,iv)
    INTEGER, OPTIONAL :: unit
    INTEGER :: iv(:)
    REAL(8) :: xc(:),yc(:),zc(:)
    REAL(8), ALLOCATABLE :: x(:),y(:),z(:)
    INTEGER :: natom,n,U_print
    u_Print=kprint
    IF(PRESENT(unit)) u_Print=unit

#ifdef HAVE_MPI
    IF(Pi_Nprocs > 1) THEN
       natom=SIZE(xc)
       ALLOCATE(x(natom),y(natom),z(natom))
       x=xc;y=yc;z=zc
       CALL PI_Gather_(x,y,z)
       IF(PI_Node == 0) WRITE(u_Print,'(i8,3e17.7)') (iv(n),x(iv(n)),y(iv(n)),z(iv(n)),n=1,SIZE(iv))
    ELSE
       IF(PI_Node == 0) WRITE(u_Print,'(i8,3e17.7)') (iv(n),xc(iv(n)),yc(iv(n)),zc(iv(n)),n=1,SIZE(iv))
    END IF
#else
    WRITE(u_Print,'(i8,3e17.7)') (iv(n),xc(iv(n)),yc(iv(n)),zc(iv(n)),n=1,SIZE(iv))
#endif    
  END SUBROUTINE PI_Write_
  SUBROUTINE PI_ErrSignal_(errmsg_f,out)
    CHARACTER(len=max_data) :: errmsg_f
    LOGICAL :: out 
    LOGICAL, ALLOCATABLE :: couts(:)
    CHARACTER(len=max_data), ALLOCATABLE :: errmsg(:)
    CHARACTER(len=max_data) :: string
    INTEGER :: cout,n

#ifdef HAVE_MPI
    IF(PI_Nprocs == 1) RETURN
    cout=out
    ALLOCATE(couts(PI_Nprocs),errmsg(PI_Nprocs))

    CALL MPI_ALLGATHER(cout,1,MPI_LOGICAL,couts,1,MPI_LOGICAL,PI_Comm_Cart,ierr)
    CALL MPI_ALLGATHER(errmsg_f,max_data,MPI_CHARACTER,errmsg,max_data,MPI_CHARACTER,PI_Comm_Cart,ierr)
    DO n=1,PI_Nprocs
       IF(.NOT. couts(n)) THEN
          WRITE(string,'(''From Node No. '',i3,'' : '')') n-1
          errmsg_f=TRIM(string)//TRIM(errmsg(n))
          out=couts(n)
          CALL Add_Errors(-1,errmsg_f)
          RETURN
       END IF
    END DO
#endif

  END SUBROUTINE PI_ErrSignal_

END MODULE PI_Collectives
