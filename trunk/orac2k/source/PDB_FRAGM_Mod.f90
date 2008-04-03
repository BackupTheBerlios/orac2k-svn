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
MODULE PDB_FRAGM_Mod
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Mar 20 2008 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program orac2k ----*

  USE INPUT_Mod, ONLY: Read_String, Parser, err_open,err_end,err_unr&
       &,err_fnf,err_args
  IMPLICIT none
  PRIVATE
  PUBLIC Initialize, Read_it, FRAGM_, FRAGM__, npdb
  CHARACTER(7), ALLOCATABLE :: beta(:),prsymb(:)
  INTEGER, SAVE :: kpdb,npdb
  CHARACTER(80), SAVE :: filename
  LOGICAL :: FRAGM__=.FALSE.
CONTAINS
  SUBROUTINE Initialize(beta_a,prsymb_a,res,ntap)
    IMPLICIT NONE
    INTEGER :: res(*),ntap
    CHARACTER(7) :: beta_a(*)
    CHARACTER(8) :: prsymb_a(*)
    INTEGER :: i
    
    ALLOCATE(beta(ntap),prsymb(ntap))
    beta=beta_a(1:ntap)
    DO i=1,ntap
       prsymb(i)=prsymb_a(res(i))
    END DO
  END SUBROUTINE Initialize
  SUBROUTINE FRAGM_(fstep,ind1,ind2,mres,grppt,xp0,yp0,zp0)
    IMPLICIT NONE 
    REAL(8) :: xp0(*),yp0(*),zp0(*),fstep
    INTEGER :: ind1(:),ind2(:),mres(:,:),grppt(:,:)
    CHARACTER(4) :: bet,bet2
    CHARACTER(3) :: ResName
    CHARACTER(1) :: Chain
    INTEGER :: n,nn,i1,i2,i,l,k
    REAL(8) :: xb,yb,zb,dr

    dr=0.0D0
    WRITE(kpdb,2) fstep
    DO n=1,ind1(1+1)-1
       i1=n
       DO i2=mres(1,i1),mres(2,i1)
          DO i=grppt(1,i2),grppt(2,i2)
             xb=xp0(i)
             yb=yp0(i)
             zb=zp0(i)
             bet(1:5)=beta(i)(1:5)
             CALL low_up(bet,5)
             DO l=1,5
                bet2(l:l)=bet(l:l)
             END DO
             ResName(1:3)=prsymb(i)(1:3)
             CALL low_up(ResName,3)
             Chain=' '
             WRITE(kpdb,1) 'ATOM ',i,bet2,ResName,i1,xb,yb,zb,dr,DFLOAT(k)
          END DO
       END DO
    END DO
    DO n=1,SIZE(ind2)
       nn=ind2(n)
       i1=ind1(1+nn)
       DO i2=mres(1,i1),mres(2,i1)
          DO i=grppt(1,i2),grppt(2,i2)
             xb=xp0(i)
             yb=yp0(i)
             zb=zp0(i)
             bet(1:5)=beta(i)(1:5)
             CALL low_up(bet,5)
             DO l=1,5
                bet2(l:l)=bet(l:l)
             END DO
             ResName(1:3)=prsymb(i)(1:3)
             CALL low_up(ResName,3)
             Chain=' '
             WRITE(kpdb,1) 'ATOM ',i,bet2,ResName,i1,xb,yb,zb,dr,DFLOAT(k)
          END DO
       END DO
    END DO
    WRITE(kpdb,'(a)')'TER  '
1   FORMAT(a5,i6,1x,a4,x,a3,1x,i5,4x,3f8.3,f8.4,f4.1)
2   FORMAT('REMARK   1 Configuration at time step ',f11.2,'      ',&
         &'                 ')
  END SUBROUTINE FRAGM_
  INCLUDE 'PDB_FRAGM__Read.f90'
END MODULE PDB_FRAGM_Mod
