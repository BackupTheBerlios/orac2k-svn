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
MODULE LA_Routines
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Tue Dec  9 2008 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*

  USE Errors,ONLY: Add_errors=>Add, Print_Errors, errmsg_f
  IMPLICIT none
  PRIVATE
  PUBLIC Matinv_
CONTAINS
  FUNCTION Matinv_(Mat,i_Mat) RESULT(out)
    LOGICAL :: out
    REAL(8) :: Mat(:,:),i_Mat(:,:)
    INTEGER :: n_dim,m_dim,info
    REAL(8), ALLOCATABLE :: Matinv(:,:),Work(:)
    INTEGER, ALLOCATABLE :: ipiv(:)


    out=.TRUE.
    n_dim=SIZE(Mat,1)
    m_dim=SIZE(Mat,2)
    IF(n_dim /= m_Dim) THEN
       errmsg_f='Cannot invert a rectangular matrix '
       CALL Add_Errors(-1,errmsg_f)
       RETURN
    END IF

    
    ALLOCATE(Matinv(n_Dim, m_Dim),ipiv(n_Dim),Work(n_Dim*20))

    CALL DCOPY(n_Dim*m_Dim,Mat,1,Matinv,1)
    CALL DGETRF(n_Dim,m_Dim,Matinv,n_Dim,ipiv,Info)
    IF(Info /= 0) THEN
       errmsg_f='DGETRF exited with info not equal zero '
       CALL Add_Errors(-1,errmsg_f)
       RETURN
    END IF
    CALL DGETRI(n_Dim,MATINV,n_Dim,Ipiv,Work,n_Dim*20,Info)
    IF(Info /= 0) THEN
       errmsg_f='DGETRI exited with info not equal zero '
       CALL Add_Errors(-1,errmsg_f)
       RETURN
    END IF
    i_Mat=Matinv
  END FUNCTION Matinv_

END MODULE LA_Routines
