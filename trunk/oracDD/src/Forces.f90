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
MODULE Forces
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

  USE Constants, ONLY: max_pars,max_data, max_char
  USE Errors, ONLY: Add_Errors=>Add, error_other, error_unr, error_args&
       &, errmsg_f, Print_Errors
  IMPLICIT none
  PRIVATE
  PUBLIC :: Force,fp_n0,fp_n1,fp_m,fp_l,fp_h,fp_ew,fs_n0,fs_n1,fs_m&
       &,fs_l,fs_h,fs_ew,Init,Radii,Cutoff,Memory, Pick
  TYPE :: Force
     REAL(8) :: x,y,z
  END type Force
  TYPE(Force), ALLOCATABLE, SAVE, TARGET :: fp_n0(:),fp_n1(:),fp_m(:),fp_l(:),fp_h(:),fp_ew(:)
  TYPE(Force), ALLOCATABLE, SAVE, TARGET :: fs_n0(:),fs_n1(:),fs_m(:),fs_l(:),fs_h(:),fs_ew(:)
  TYPE :: cutoff
     REAL(8) :: inn,out,shell,update
  END type cutoff
  TYPE(Cutoff), ALLOCATABLE, SAVE :: Radii(:)
CONTAINS
  SUBROUTINE Memory(natom)
    INTEGER :: natom
    ALLOCATE(fp_n0(natom),fp_n1(natom),fp_m(natom),fp_l(natom),fp_h(natom),fp_ew(natom))
  END SUBROUTINE Memory
  SUBROUTINE Init(NShell,Rcut_Table)
    INTEGER :: NShell
    REAL(8) :: Rcut_Table(3,3)
    REAL(8), DIMENSION(3) :: rcut_i,rcut_s,rcut_o,rcut_u 
    REAL(8) :: rneigh=1.5D0
    INTEGER :: n,m
    CHARACTER(len=max_char) :: lab0
    INTEGER, ALLOCATABLE :: indx(:)
    REAL(8), ALLOCATABLE :: vect(:),vect1(:)
    
    rcut_i=(/4.7D0, 7.3D0, 9.7D0/)
    rcut_s=(/0.3D0,0.3D0,0.3D0/)
    rcut_u=(/0.5D0,0.5D0,1.5D0/)

    SELECT CASE(NShell)
    CASE(1)
       rcut_i(1)=Rcut_Table(1,1)
       rcut_u(1)=rneigh

       IF(Rcut_Table(1,2) /= -1.0D0) THEN
          rcut_s(1)=Rcut_Table(1,2)
       ELSE
          rcut_s(1)=0.0D0
       END IF
       IF(Rcut_Table(1,3) /= -1.0D0) THEN
          rcut_u(1)=Rcut_Table(1,3)
       ELSE
          rcut_u(1)=0.0D0
       END IF
    CASE(2)
       rcut_i(1)=Rcut_Table(1,1)
       IF(Rcut_Table(1,2) /= -1.0D0) THEN
          rcut_s(1)=Rcut_Table(1,2)
       ELSE
          rcut_s(1)=0.0D0
       END IF
       IF(Rcut_Table(1,3) /= -1.0D0) THEN
          rcut_u(1)=Rcut_Table(1,3)
       ELSE
          rcut_u(1)=0.0D0
       END IF

       rcut_i(2)=Rcut_Table(2,1)
       rcut_u(2)=rneigh
       IF(Rcut_Table(2,2) /= -1.0D0) THEN
          rcut_s(2)=Rcut_Table(2,2)
       ELSE
          rcut_s(2)=0.0D0
       END IF
       IF(Rcut_Table(2,3) /= -1.0D0) THEN
          rcut_u(2)=Rcut_Table(2,3)
       ELSE
          rcut_u(2)=0.0D0
       END IF
    CASE(3)
       rcut_i(1)=Rcut_Table(1,1)
       IF(Rcut_Table(1,2) /= -1.0D0) THEN
          rcut_s(1)=Rcut_Table(1,2)
       ELSE
          rcut_s(1)=0.0D0
       END IF
       IF(Rcut_Table(1,3) /= -1.0D0) THEN
          rcut_u(1)=Rcut_Table(1,3)
       ELSE
          rcut_u(1)=0.0D0
       END IF

       rcut_i(2)=Rcut_Table(2,1)
       IF(Rcut_Table(2,2) /= -1.0D0) THEN
          rcut_s(2)=Rcut_Table(2,2)
       ELSE
          rcut_s(2)=0.0D0
       END IF
       IF(Rcut_Table(2,3) /= -1.0D0) THEN
          rcut_u(2)=Rcut_Table(2,3)
       ELSE
          rcut_u(2)=0.0D0
       END IF

       rcut_i(3)=Rcut_Table(3,1)
       IF(Rcut_Table(3,2) /= -1.0D0) THEN
          rcut_s(3)=Rcut_Table(3,2)
       ELSE
          rcut_s(3)=0.0D0
       END IF
       IF(Rcut_Table(3,3) /= -1.0D0) THEN
          rcut_u(3)=Rcut_Table(3,3)
       ELSE
          rcut_u(3)=0.0D0
       END IF
    CASE DEFAULT
       WRITE(lab0,'(i2)') NShell
       errmsg_f=' Can run only with three shells, whereas '&
            &//TRIM(lab0)//' were requested. '
       CALL Add_Errors(-1,errmsg_f)
       RETURN
    END SELECT

    ALLOCATE(Radii(NShell))

    IF(NShell > 1) THEN
       ALLOCATE(vect(NShell),indx(NShell))
       
       vect(:)=rcut_i(1:NShell)
       
       CALL INDEXX(NShell,vect,indx)
       DO n=1,NShell
          m=indx(n)
          Radii(n) % inn = rcut_i(m)
          Radii(n) % shell = rcut_s(m)
          Radii(n) % update = rcut_u(m)
          Radii(n) % out = rcut_i(m)+rcut_s(m)
       END DO
    ELSE
       DO n=1,NShell
          Radii(n) % inn = rcut_i(n)
          Radii(n) % shell = rcut_s(n)
          Radii(n) % update = rcut_u(n)
          Radii(n) % out = rcut_i(n)+rcut_s(n)
       END DO
    END IF
  END SUBROUTINE Init
  FUNCTION Pick(level) RESULT(out)
    INTEGER :: level
    TYPE(Force), POINTER :: out(:)
    SELECT CASE(level)
    CASE(_N0_)
       out=>fp_n0
    CASE(_N1_)
       out=>fp_n1
    CASE(_M_)
       out=>fp_m
    CASE(_L_)
       out=>fp_l
    CASE(_H_)
       out=>fp_h
    CASE DEFAULT
       errmsg_f='Only five level RESPA is allowed'
       CALL Add_Errors(-1,errmsg_f)
       CALL Print_Errors()
       STOP
    END SELECT
  END FUNCTION Pick
END MODULE Forces
