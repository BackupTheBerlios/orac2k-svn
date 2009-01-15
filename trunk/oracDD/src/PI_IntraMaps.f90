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
MODULE PI_IntraMaps
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Nov 21 2008 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*

#ifdef HAVE_MPI
  USE mpi
#endif
  USE PI_
  USE SystemTpg
  USE Groups
  USE Atom
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, error_args,&
       & errmsg_f
  USE IndBox, ONLY: IndBox_g_t,IndBox_g_p
  IMPLICIT none
  PRIVATE
  PUBLIC IntraMaps_n0_, IntraMaps_n1_, Map_n0, Map_n1
  INTEGER, ALLOCATABLE, TARGET, SAVE :: Map_n0(:),Map_n1(:)
  INTEGER, SAVE :: Init_Calls=0
  
CONTAINS
  SUBROUTINE IntraMaps_n0_
    INTEGER :: n,n0,AtSt,AtEn,q,count0,nn
    INTEGER, ALLOCATABLE :: ind_o(:)
    LOGICAL :: ok
    
    IF(.NOT. ALLOCATED(IndBox_g_p)) THEN
       errmsg_f='Internal Error. PI_IntraMaps must be called after IndBox'
       CALL Add_Errors(-1,errmsg_f)
       CALL Print_Errors()
       STOP
    END IF

    n0=SIZE(IndBox_g_p)
    ALLOCATE(ind_o(n0))
    ind_o=0
    count0=0
    DO nn=1,n0
       n=IndBox_g_p(nn)
       AtSt=Groupa(n) % AtSt
       AtEn=Groupa(n) % AtEn
       ok=.FALSE.
       DO q=AtSt,AtEn
          IF(Choose_Atom(q, Tpg % Bonds, Atoms_Tpg(q) % Bonds)) THEN
             ok=.TRUE.
             EXIT
          END IF
          IF(Choose_Atom(q, Tpg % Angles, Atoms_Tpg(q) % Angles)) THEN
             ok=.TRUE.
             EXIT
          END IF
          IF(Choose_Atom(q, Tpg % Imph, Atoms_Tpg(q) % Imph)) THEN
             ok=.TRUE.
             EXIT
          END IF
       END DO
       IF(ok) THEN
          count0=count0+1
          ind_o(count0)=n
       END IF
    END DO

    IF(ALLOCATED(Map_n0)) DEALLOCATE(Map_n0)
    ALLOCATE(Map_n0(count0))
    Map_n0=ind_o(1:count0)
  END SUBROUTINE IntraMaps_n0_
  SUBROUTINE IntraMaps_n1_
    INTEGER :: n,n0,AtSt,AtEn,q,count0,nn
    INTEGER, ALLOCATABLE :: ind_o(:)
    LOGICAL :: ok

    IF(.NOT. ALLOCATED(IndBox_g_p)) THEN
       errmsg_f='Internal Error. PI_IntraMaps must be called after IndBox'
       CALL Add_Errors(-1,errmsg_f)
       CALL Print_Errors()
       STOP
    END IF

    n0=SIZE(IndBox_g_p)

    ALLOCATE(ind_o(n0))
    ind_o=0
    count0=0
    DO nn=1,n0
       n=IndBox_g_p(nn)
       AtSt=Groupa(n) % AtSt
       AtEn=Groupa(n) % AtEn
       ok=.FALSE.
       DO q=AtSt,AtEn
          IF(Choose_Atom(q, Tpg % Dihed, Atoms_Tpg(q) % Dihed)) THEN
             ok=.TRUE.
             EXIT
          END IF
          IF(Choose_Atom(q, Tpg % Int14, Atoms_Tpg(q) % Int14)) THEN
             ok=.TRUE.
             EXIT
          END IF
       END DO
       IF(ok) THEN
          count0=count0+1
          ind_o(count0)=n
       END IF
    END DO
    IF(ALLOCATED(Map_n1)) DEALLOCATE(Map_n1)
    ALLOCATE(Map_n1(count0))
    Map_n1=ind_o(1:count0)
  END SUBROUTINE IntraMaps_n1_
  FUNCTION Choose_Atom(q, Tpga, Tpgb) RESULT(out)
    LOGICAL :: out
    INTEGER :: q,Tpga(:,:), Tpgb(:)
    
    INTEGER :: n,p
    LOGICAL :: ok
    INTEGER :: r,t,s,ss,ssk
    

    n=SIZE(Tpgb)
    p=SIZE(Tpga,1)
    ok=.FALSE.
    
    DO r=1,n
       t=Tpgb(r)
       DO s=1,p
          ss=Tpga(s,t)
          IF(ss == q) CYCLE
          ssk=Atoms(ss) % knwn
          IF(ssk == 0) THEN
             ok=.TRUE.
             EXIT
          END IF
       END DO
       IF(ok) EXIT
    END DO
    out=ok
  END FUNCTION Choose_Atom
END MODULE PI_IntraMaps
