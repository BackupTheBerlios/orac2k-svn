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
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Jan 23 2009 -                                     *
!!$                                                                      *
!!$***********************************************************************

#define _IN_     IndBox_a_P(:)

!!$---- This module is part of the program oracDD ----*
SUBROUTINE ATOM__KinCompute
  REAL(8), POINTER :: vx(:),vy(:),vz(:),mass(:)
  REAL(8) :: ek_Slt,ek_slv

  vx=>Atoms(:) % vx
  vy=>Atoms(:) % vy
  vz=>Atoms(:) % vz
  mass=>Atoms(:) % mass

  ek_Slt=SUM(mass(_IN_)*(vx(_IN_)**2+vy(_IN_)**2+vz(_IN_)**2)&
       &,Atoms(_IN_) % Id_Slv == 1) 
  ek_Slv=SUM(mass(_IN_)*(vx(_IN_)**2+vy(_IN_)**2+vz(_IN_)**2)&
       &,Atoms(_IN_) % Id_Slv == 2) 

  ek_Slv=0.5*ek_Slv*efact
  ek_Slv=0.5*ek_Slv*efact
  CALL EN_Kinetic_(ek_Slv,ek_Slt)

END SUBROUTINE ATOM__KinCompute
