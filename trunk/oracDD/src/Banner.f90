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
MODULE Banner
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Wed Jan 24 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*

  IMPLICIT none
  CHARACTER(len=80), DIMENSION(31) :: Banner__Page=&
       &(/&
       &'================================================================================',&
       &'=                                                                              =',&
       &'=                                                                              =',&
       &'=                                                                              =',&
       &'=                                            oooooooooo.   oooooooooo.         =',&
       &'=                                            `888''   `Y8b  `888''   `Y8b        =',&
       &'=      .ooooo.  oooo d8b  .oooo.    .ooooo.   888      888  888      888       =',&
       &'=     d88'' `88b `888""8P `P  )88b  d88'' `"Y8  888      888  888      888       =',&
       &'=     888   888  888      .oP"888  888        888      888  888      888       =',&
       &'=     888   888  888     d8(  888  888   .o8  888     d88''  888     d88''       =',&
       &'=     `Y8bod8P'' d888b    `Y888""8o `Y8bod8P'' o888bood8P''   o888bood8P''         =',&
       &'=                                                                              =',&
       &'=                              (Version 0.1)                                   =',&
       &'=                                                                              =',&
       &'=     Copyright (C) 2006-2007 Massimo Marchi <Massimo.Marchi at cea.fr>        =',&
       &'=                                                                              =',&
       &'=                   Commissariat a l''Energie Atomique                          =',&
       &'=                       Centre d''Etudes de Saclay                              =',&
       &'=                         Gif-sur-Yvette, FRANCE                               =',&
       &'=                                                                              =',&
       &'=     "A parallel computer program to simulate and model complex               =',&
       &'=      molecular systems"                                                      =',&
       &'=                                                                              =',&
       &'=     This software is governed by the CeCILL license under French law         =',&
       &'=     and abiding by the rules of distribution of free software.  You can      =',&
       &'=     use, modify and/ or redistribute the software under the terms of         =',&
       &'=     the CeCILL icense as circulated by CEA, CNRS and INRIA at the            =',&
       &'=     following URL "http://www.cecill.info".                                  =',&
       &'=                                                                              =',&
       &'=                                                                              =',&
       &'================================================================================'/)
CONTAINS
  SUBROUTINE Banner_
    WRITE(*,'(a80)') Banner__Page
  END SUBROUTINE Banner_
END MODULE Banner
