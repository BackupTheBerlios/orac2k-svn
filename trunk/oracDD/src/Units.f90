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
MODULE UNITS
  REAL(8), PARAMETER :: pi=3.1415926535897931D0, twopi=2.0D0*pi, avogad=6.0225D23&
       &, boltz=1.38054D-23,gascon=8.3143D0,planck=6.6256D-34,elechg=1.602D-19&
       &,epso=8.854D-12,boxl=2.0D0,unitm=1.0D0/(avogaD*1000.0D0),unitl=1.0D-10&
       &,unitt=1.D-15,unite=unitm*(unitl/unitt)**2&
       &,unitc=4.0D0*pi*epso*unitl*unite/(elechg*elechg)&
       &,unitp=(unite/unitl**3)/1.0D6,efact=unite*avogad&
       &,hartree = 4.35981 * 1.0D-18 &
       &,lbohr = 0.52917706D0
  REAL(8), SAVE :: LJ_Fact,unitepot,unitefield
CONTAINS
  SUBROUTINE Units_
    LJ_fact=1.0D0/(2.0D0**(1.0D0/6.0D0))
    unitepot=unite/unitc**(0.5D0)/hartree
    unitefield=unitepot*lbohr 
  END SUBROUTINE Units_
END MODULE UNITS
