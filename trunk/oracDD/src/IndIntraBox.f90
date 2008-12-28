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
MODULE IndIntraBox
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Tue Aug  5 2008 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*

#include "config.h"

#ifdef HAVE_MPI
  USE mpi
#endif
  USE PI_
  USE PrmUtilities
  USE SystemTpg
  USE SystemPrm
  USE PI_IntraMaps
  USE Groups
  USE Atom
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f, errmsg_w
  IMPLICIT none
  PRIVATE
  PUBLIC IndIntraBox_n0_, IntraParam, Param_Bonds,Param_Angles,Param_Imph&
       &,Param_Dihed,Param_Constr,Indx_Bonds,Indx_Angles,Indx_Imph&
       &,Indx_Dihed,Indx_Int14,Indx_Constr, IndIntraBox_n1_

  INTEGER :: natom_local
  INTEGER, ALLOCATABLE, SAVE :: indBox_a_p(:)
  LOGICAL, ALLOCATABLE, SAVE :: oks(:),okt(:)
  TYPE :: IntraParam
     REAL(8), ALLOCATABLE :: Pot(:)
  END type IntraParam
  TYPE(IntraParam), ALLOCATABLE, SAVE :: Param_Bonds(:),Param_Angles(:)&
       &,Param_Imph(:),Param_Dihed(:),Param_Constr(:)
  INTEGER, ALLOCATABLE, SAVE :: Indx_Bonds(:,:),Indx_Angles(:,:)&
       &,Indx_Imph(:,:),Indx_Dihed(:,:),Indx_Int14(:,:),Indx_Constr(:,:)
  INTEGER, SAVE :: dummys=0
CONTAINS
  FUNCTION IndIntraBox_n0_() RESULT(out)
    LOGICAL :: out
    INTEGER :: n,m
    INTEGER, ALLOCATABLE :: via(:)
    out=.FALSE.

    IF(.NOT. IndIntraBox_()) RETURN
 
    CALL IntraPot( Prm % Bonds, Tpg % Bonds, Indx_Bonds, Param_Bonds)
    CALL IntraPot( Prm % Constr, Tpg % Bonds, Indx_Constr, Param_Constr)
    CALL IntraPot( Prm % Angles, Tpg % Angles, Indx_Angles, Param_Angles)
    CALL IntraPot( Prm % Imph, Tpg % Imph, Indx_Imph, Param_Imph)

    out=.TRUE.
  END FUNCTION IndIntraBox_n0_
  FUNCTION IndIntraBox_n1_() RESULT(out)
    LOGICAL :: out
    INTEGER :: n,m
    out=.FALSE.

    IF(.NOT. IndIntraBox_()) RETURN
 
    CALL IntraPot(Prm % Dihed, Tpg % Dihed, Indx_Dihed, Param_Dihed)
    CALL IntraPot14(Tpg % Int14, Indx_Int14)

    out=.TRUE.
  END FUNCTION IndIntraBox_n1_

  SUBROUTINE IntraPot(Prm, Tpg, Indx, Param)
    TYPE(SystemPrm__Chain) :: Prm(:)
    INTEGER :: Tpg(:,:)
    INTEGER, ALLOCATABLE :: Indx(:,:)
    TYPE(IntraParam), ALLOCATABLE :: Param(:)
    
    INTEGER, ALLOCATABLE :: via(:),ind_o(:)
    INTEGER :: n_Tpg,m_Prm,count0,count1,n,nn,ng,n1,n2,nnn


    m_Prm=SIZE(Prm)
    n_Tpg=SIZE(Tpg,1)

    ALLOCATE(via(n_Tpg),ind_o(m_Prm))

    ind_o=0
    count0=0
    DO nn=1,m_Prm
       n=Prm (nn) % pt
       IF(n < 0) CYCLE
       via=Tpg(:,n)
       n1=COUNT(oks(via(:)))
       IF(n1 == 0) CYCLE
       n2=COUNT(okt(via(:)))
       IF(n1+n2 == n_Tpg) THEN
          count0=count0+1
          ind_o(count0)=nn
       END IF
    END DO
    IF(ALLOCATED(Indx)) DEALLOCATE(Indx)
    IF(ALLOCATED(Param)) DEALLOCATE(Param)

    ALLOCATE(Indx(n_Tpg,count0),Param(count0))

    DO nnn=1,count0
       nn=ind_o(nnn)
       n=Prm (nn) % pt
       via=Tpg (:,n)
       ng=SIZE(Prm (nn) % g)
       ALLOCATE(Param(nnn) % pot(ng))
       Indx(:,nnn)=via
       Param(nnn) % pot=Prm(nn) % g
    END DO
    CALL MPI_ALLREDUCE(count0,count1,1,MPI_INTEGER4,MPI_SUM,PI_Comm_Cart,ierr)
  END SUBROUTINE IntraPot
  SUBROUTINE IntraPot14(Tpg, Indx)
    INTEGER :: Tpg(:,:)
    INTEGER, ALLOCATABLE :: Indx(:,:)
    
    INTEGER, ALLOCATABLE :: via(:),ind_o(:)
    INTEGER :: n_Tpg,m_Tpg,count0,n,nn,ng,n1,n2
    
    
    n_Tpg=SIZE(Tpg,1)
    m_Tpg=SIZE(Tpg,2)

    ALLOCATE(via(n_Tpg),ind_o(m_Tpg))
    count0=0
    DO n=1,m_Tpg
       via=Tpg(:,n)
       n1=COUNT(oks(via(:)))
       IF(n1 == 0) CYCLE
       n2=COUNT(okt(via(:)))
       IF(n1+n2 == n_Tpg) THEN
          count0=count0+1
          ind_o(count0)=n
       END IF
    END DO
    IF(ALLOCATED(Indx)) DEALLOCATE(Indx)
    ALLOCATE(Indx(n_Tpg,count0))

    DO nn=1,count0
       n=ind_o(nn)
       via=Tpg (:,n)
       Indx(:,nn)=via
    END DO    
  END SUBROUTINE IntraPot14

  FUNCTION IndIntraBox_() RESULT(out)
    LOGICAL :: out
    INTEGER :: n,m,count_a_p,mm
    
      
    count_a_p=0
    count_a_p=COUNT(Groupa(Atoms(:) % Grp_No) % knwn == 1)

    IF(ALLOCATED(IndBox_a_p)) DEALLOCATE(IndBox_a_p)
    IF(ALLOCATED(oks)) DEALLOCATE(oks)
    IF(ALLOCATED(okt)) DEALLOCATE(okt)

    ALLOCATE(IndBox_a_p(count_a_p),oks(SIZE(Atoms)),okt(SIZE(Atoms)))
    
    count_a_p=0 
    oks=.FALSE.
    okt=.FALSE.
    DO n=1,SIZE(Atoms)
       mm=Atoms(n) % Grp_No
       m=Groupa(mm) % knwn 
       IF(m == 1) THEN
          count_a_p=count_a_p+1
          IndBox_a_p(count_a_p)=n
       END IF
       IF(m == 1) oks(n)=.TRUE.
       IF(m == 2) okt(n)=.TRUE.
    END DO
    out=count_a_p /= 0
    IF(.NOT. out) THEN
       errmsg_f='No Primary Atoms for Intramolecular interaction found in the unit box'
       CALL Add_Errors(-1,errmsg_f)
    END IF
  END FUNCTION IndIntraBox_
END MODULE IndIntraBox
