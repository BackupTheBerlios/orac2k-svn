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
MODULE Intra
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

#define _INIT_EXCHANGE_     0
#define _EXCHANGE_ONLY_     1
  USE PI_
  USE CHARMM_Intra, ONLY: CHARMM_Bonds=>Bonds,CHARMM_Angles=>Angles&
       &,CHARMM_Dihed=>Dihed, CHARMM_Imph=>Imph, CHARMM_Int14=>Int14&
       &,CHARMM_Bonds_=>Bonds_,CHARMM_Angles_=>Angles_,CHARMM_Dihed_&
       &=>Dihed_, CHARMM_Imph_=>Imph_, CHARMM_Int14coul_=>Int14coul_ &
       &, CHARMM_Int14conf_=>Int14conf_  
  USE Keyword, ONLY: charmm
  USE Forces, ONLY: fp_n0,fp_n1
  USE IndIntraBox, ONLY: Param_Bonds,Indx_Bonds,Param_Angles&
       &,Indx_Angles,Param_Dihed,Indx_Dihed,Param_Imph,Indx_Imph&
       &,Indx_Int14
  USE IndBox
  USE Atom
  USE Groups
  USE PI_Communicate, ONLY: PI__FoldIntra
  IMPLICIT none
  PRIVATE
  PUBLIC Intra_n0_, Intra_n1_
  REAL(8), SAVE :: ubond_Slv,ubond_Slt,ubend_Slv,ubend_Slt,uitors_Slv&
       &,uitors_Slt,uptors_Slv,uptors_Slt,ucoulint14_Slv,&
       &ucoulint14_Slt,uconfint14_Slv,uconfint14_Slt
CONTAINS
  SUBROUTINE Intra_n0_(init)
    INTEGER :: init
    INTEGER :: n,nn
    fp_n0(:) % x  = 0.0_8
    fp_n0(:) % y  = 0.0_8
    fp_n0(:) % z  = 0.0_8
    IF(charmm) THEN
       CALL Intra_(CHARMM_Bonds,fp_n0 % x,fp_n0 % y,fp_n0 % z)

       CALL Intra_(CHARMM_Angles,fp_n0 % x,fp_n0 % y,fp_n0 % z)
       CALL Intra_(CHARMM_Imph,fp_n0 % x,fp_n0 % y,fp_n0 % z)

       CALL Energy_(CHARMM_Bonds_,ubond_Slv,ubond_Slt)

       CALL Energy_(CHARMM_Angles_,ubend_Slv,ubend_Slt)
       CALL Energy_(CHARMM_Imph_,uitors_Slv,uitors_Slt)

       CALL PI__FoldIntra(fp_n0,1,init)

       IF(PI_Node_Cart == 0) WRITE(*,*) 'u1 ',ubond_Slv,ubond_Slt
       IF(PI_Node_Cart == 0) WRITE(*,*) 'u2 ',ubend_Slv,ubend_Slt
       IF(PI_Node_Cart == 0) WRITE(*,*) 'u3 ',uitors_Slv,uitors_Slt
!!$       DO n=1,SIZE(Atoms)
!!$          IF(Groupa(Atoms(n) % Grp_No) % knwn /= 1) CYCLE
!!$          WRITE(100+PI_Node_Cart,'(i8,4x,3e14.5)') n,fp_n0 (n) %  x,&
!!$               & fp_n0 (n) %  y, fp_n0 (n) %  z
!!$       END DO
    END IF
  END SUBROUTINE Intra_n0_
  SUBROUTINE Intra_n1_(init)
    INTEGER :: init
    INTEGER :: n,nn
    fp_n1(:) % x  = 0.0_8
    fp_n1(:) % y  = 0.0_8
    fp_n1(:) % z  = 0.0_8
    IF(charmm) THEN
       CALL Intra_(CHARMM_Dihed,fp_n1 % x,fp_n1 % y,fp_n1 % z)
       CALL Intra_(CHARMM_Int14,fp_n1 % x,fp_n1 % y,fp_n1 % z)
       CALL Energy_(CHARMM_Dihed_,uptors_Slv,uptors_Slt)
       CALL Energy_(CHARMM_Int14coul_,ucoulint14_Slv,ucoulint14_Slt)
       CALL Energy_(CHARMM_Int14conf_,uconfint14_Slv,uconfint14_Slt)

       CALL PI__FoldIntra(fp_n1,2,init)

       IF(PI_Node_Cart == 0) WRITE(*,*) 'u4 ',uptors_Slv,uptors_Slt
       IF(PI_Node_Cart == 0) WRITE(*,*) 'u5 ',uconfint14_Slv,ucoulint14_Slv
       IF(PI_Node_Cart == 0) WRITE(*,*) 'u6 ',uconfint14_Slt,ucoulint14_Slt
!!$       DO nn=1,SIZE(IndBox_a_p)
!!$          n=IndBox_a_t(IndBox_a_p(nn))
!!$          WRITE(200+PI_Node_Cart,'(i8,4x,3e14.5)') n, fp_n1 (n) %  x,&
!!$               & fp_n1 (n) %  y, fp_n1 (n) %  z
!!$       END DO
    END IF
  END SUBROUTINE Intra_n1_

  SUBROUTINE Intra_(Routine,fpx,fpy,fpz)
    REAL(8) :: fpx(:),fpy(:),fpz(:)
    INTERFACE
       SUBROUTINE Routine(fpx,fpy,fpz)
         REAL(8) :: fpx(:),fpy(:),fpz(:)
       END SUBROUTINE Routine
    END INTERFACE
    CALL Routine(fpx,fpy,fpz)
  END SUBROUTINE Intra_
  SUBROUTINE Energy_(Routine,u_slv,u_slt)
    REAL(8) :: u_slv,u_slt
    INTERFACE
       SUBROUTINE Routine(u_Slv,u_Slt)
         REAL(8) :: u_Slv,u_Slt
       END SUBROUTINE Routine
    END INTERFACE
    CALL Routine(u_slv,u_slt)
  END SUBROUTINE Energy_
END MODULE Intra
