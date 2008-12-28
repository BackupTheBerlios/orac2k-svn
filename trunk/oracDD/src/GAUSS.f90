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
MODULE Gauss
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Mon Dec  8 2008 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*
  IMPLICIT None
  PUBLIC cgauss
CONTAINS
!!$==========================================================================
!!$ Program to produce exponentially correlated colored (Gaussian) noise.
!!$ based on Fox et al Physical Review A vol.38(1988)5938 and 
!!$ modification of GASDEV from Numerical Recipes for Fortran(2nd ed.pg279)
!!$
!!$ CAPE is capital E in the article by Fox et. al.
!!$ PREV is the previous value of CGAUSS used in the next iteration
!!$ L1ME2 is the main parameters causing colored noise in Fox et al
!!$       and represents (lamda*(1-E)^2). Ditto for H in that article.
!!$
!!$ routine is illustrated in Double Precision in case it is needed in this
!!$ mode, otherwise all Double Precision variables maybe changed to REAL
!!$ but the corresponding changes must be made to CGAUS0 and the calling
!!$ programs.


  FUNCTION cgauss() RESULT(out)
    REAL(8) :: out
    INTEGER,SAVE :: iset=0
    REAL(8) ::  fac,rsq,v1,v2,h
    REAL(8), SAVE :: gset,prev=0.0_8
    REAL(8), EXTERNAL :: Duni

    if (iset.eq.0) then
1      v1=2.0_8*Duni()-1.0_8
       v2=2.0_8*Duni()-1.0_8
       rsq=v1**2+v2**2
       IF(rsq >= 1.0_8 .OR. rsq == 0.0_8) GOTO 1

       fac=SQRT(-2.0_8*LOG(rsq)/rsq)
       gset=v1*fac
       h=v2*fac
       iset=1
    else
       h=gset
       iset=0
    endif

    out=h       !in integration is previously set in PARAM

  END FUNCTION cgauss
END MODULE Gauss
