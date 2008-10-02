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

  IMPLICIT none
  PRIVATE
  PUBLIC :: Force,fp_n0,fp_n1,fp_m,fp_l,fp_h,fp_ew,fs_n0,fs_n1,fs_m&
       &,fs_l,fs_h,fs_ew,rcut_i,rcut_o,rcut_s,rcut_u
  TYPE :: Force
     REAL(8) :: x,y,z
  END type Force
  TYPE(Force), ALLOCATABLE, SAVE :: fp_n0(:),fp_n1(:),fp_m(:),fp_l(:),fp_h(:),fp_ew(:)
  TYPE(Force), ALLOCATABLE, SAVE :: fs_n0(:),fs_n1(:),fs_m(:),fs_l(:),fs_h(:),fs_ew(:)
  REAL(8), DIMENSION(3), SAVE :: rcut_i=(/4.7D0, 7.3D0, 9.7D0/)&
       &,rcut_s=(/0.3D0,0.3D0,0.3D0/),rcut_o=(/4.7D0+0.3D0, 7.3D0&
       &+0.3D0, 9.7D0+0.3D0/), rcut_u=(/0.5D0,0.5D0,1.5D0/)
END MODULE Forces
