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
MODULE CHARMM_Intra
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Nov 28 2008 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*

#ifdef HAVE_MPI
  USE mpi
#endif
  USE Errors,ONLY: Add_errors=>Add, Print_Errors, errmsg_f
  USE Units
  USE IndIntraBox, ONLY: Param_Bonds,Indx_Bonds,Param_Angles&
       &,Indx_Angles,Param_Dihed,Indx_Dihed,Param_Imph,Indx_Imph&
       &,Indx_Int14
  USE Potential
  USE Atom
  USE Groups
  USE PI_
  USE LennardJones, ONLY: LennardJones__Par
  IMPLICIT none
  PRIVATE
  PUBLIC Bonds,Angles,Dihed,Imph,Int14,Bonds_,Angles_,Dihed_,Imph_&
       &,Int14coul_,Int14conf_ 
  REAL(8), SAVE :: Ubond_slv,Ubond_slt,Ubend_slv,Ubend_slt,Uitors_slv&
       &,Uitors_Slt,Uptors_Slv,Uptors_Slt,Uint14coul_Slv&
       &,Uint14coul_Slt,Uint14conf_Slv,Uint14conf_Slt
  INTEGER, SAVE :: Calls=0
  REAL(8), SAVE :: Conv_Fact,Degree_To_Rad=pi/180.0_8
CONTAINS
  INCLUDE "CHARMM_Bonds.f90" 
  INCLUDE "CHARMM_Angles.f90" 
  INCLUDE "CHARMM_Imph.f90" 
  INCLUDE "CHARMM_Dihed.f90" 
  INCLUDE "CHARMM_Int14.f90"
  SUBROUTINE Bonds_(u_Slv,u_Slt)
    REAL(8) :: u_Slv,u_Slt
#ifdef HAVE_MPI
    CALL MPI_ALLREDUCE(Ubond_slv,u_slv,1,MPI_REAL8,MPI_SUM,PI_Comm_Cart,ierr)
    CALL MPI_ALLREDUCE(Ubond_slt,u_slt,1,MPI_REAL8,MPI_SUM,PI_Comm_Cart,ierr)
#else
    u_slv=Ubond_slv
    u_slt=Ubond_slt
#endif
  END SUBROUTINE Bonds_
  SUBROUTINE Angles_(u_Slv,u_Slt)
    REAL(8) :: u_Slv,u_Slt
#ifdef HAVE_MPI
    CALL MPI_ALLREDUCE(Ubend_slv,u_slv,1,MPI_REAL8,MPI_SUM,PI_Comm_Cart,ierr)
    CALL MPI_ALLREDUCE(Ubend_slt,u_slt,1,MPI_REAL8,MPI_SUM,PI_Comm_Cart,ierr)
#else
    u_slv=Ubend_slv
    u_slt=Ubend_slt
#endif
  END SUBROUTINE Angles_

  SUBROUTINE Imph_(u_Slv,u_Slt)
    REAL(8) :: u_Slv,u_Slt
#ifdef HAVE_MPI
    CALL MPI_ALLREDUCE(Uitors_slv,u_slv,1,MPI_REAL8,MPI_SUM,PI_Comm_Cart,ierr)
    CALL MPI_ALLREDUCE(Uitors_slt,u_slt,1,MPI_REAL8,MPI_SUM,PI_Comm_Cart,ierr)
#else
    u_slv=Uitors_slv
    u_slt=Uitors_slt
#endif
  END SUBROUTINE Imph_
  SUBROUTINE Dihed_(u_Slv,u_Slt)
    REAL(8) :: u_Slv,u_Slt
#ifdef HAVE_MPI
    CALL MPI_ALLREDUCE(Uptors_slv,u_slv,1,MPI_REAL8,MPI_SUM,PI_Comm_Cart,ierr)
    CALL MPI_ALLREDUCE(Uptors_slt,u_slt,1,MPI_REAL8,MPI_SUM,PI_Comm_Cart,ierr)
#else
    u_slv=Uptors_slv
    u_slt=Uptors_slt
#endif
  END SUBROUTINE Dihed_
  SUBROUTINE Int14coul_(u_Slv,u_Slt)
    REAL(8) :: u_Slv,u_Slt
#ifdef HAVE_MPI
    CALL MPI_ALLREDUCE(Uint14coul_slv,u_slv,1,MPI_REAL8,MPI_SUM,PI_Comm_Cart,ierr)
    CALL MPI_ALLREDUCE(Uint14coul_slt,u_slt,1,MPI_REAL8,MPI_SUM,PI_Comm_Cart,ierr)
#else
    u_slv=Uint14coul_slv
    u_slt=Uint14coul_slt
#endif
  END SUBROUTINE Int14coul_
  SUBROUTINE Int14conf_(u_Slv,u_Slt)
    REAL(8) :: u_Slv,u_Slt
#ifdef HAVE_MPI
    CALL MPI_ALLREDUCE(Uint14conf_slv,u_slv,1,MPI_REAL8,MPI_SUM,PI_Comm_Cart,ierr)
    CALL MPI_ALLREDUCE(Uint14conf_slt,u_slt,1,MPI_REAL8,MPI_SUM,PI_Comm_Cart,ierr)
#else
    u_slv=Uint14conf_slv
    u_slt=Uint14conf_slt
#endif
  END SUBROUTINE Int14conf_

END MODULE CHARMM_Intra
