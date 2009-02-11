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

!!$---- This module is part of the program oracDD ----*
SUBROUTINE Write_it_(Time)
  CHARACTER(len=max_char), SAVE :: str,str0
  REAL(8) :: Time
  INTEGER :: Iflag,nlen
  REAL(8), SAVE :: fact=efact/1000.0_8
  REAL(8) :: iop

  iop=12.0D0
  CALL Get_Energies_
  WRITE(kprint,'(a)') REPEAT('-',107)
  str0=' '
  str0=pipe//ctime//TRIM(NiceWrite_R8(Time))
  str0=TRIM(str0)//ctot//TRIM(NiceWrite_R8(_fact*fact*Total % Tot))
  str0=TRIM(str0)//ctpot//TRIM(NiceWrite_R8(_fact*fact*TotPot % Tot))
  str0=TRIM(str0)//ccdir//TRIM(NiceWrite_R8(_fact*fact*SUM(Coul_Dir % Tot)))
  str0=TRIM(str0)//clj//TRIM(NiceWrite_R8(_fact*fact*SUM(Lj % Tot)))//pipe
  WRITE(kprint,'(a)') TRIM(str0)
  str0=' '
  str0=pipe//ccrec//TRIM(NiceWrite_R8(_fact*fact*Coul_rec % Tot))
  str0=TRIM(str0)//cbond//TRIM(NiceWrite_R8(_fact*fact*Bond % Tot))
  str0=TRIM(str0)//cstre//TRIM(NiceWrite_R8(_fact*fact*stretch % Tot))
  str0=TRIM(str0)//cangle//TRIM(NiceWrite_R8(_fact*fact*Angle % Tot))
  str0=TRIM(str0)//cdihed//TRIM(NiceWrite_R8(_fact*fact*Dihed % Tot))//pipe
  WRITE(kprint,'(a)') TRIM(str0)
  str0=' '
  str0=pipe//cimph//TRIM(NiceWrite_R8(_fact*fact*Imph % Tot))
  str0=TRIM(str0)//c14lj//TRIM(NiceWrite_R8(_fact*fact*Int14LJ % Tot))
  str0=TRIM(str0)//c14coul//TRIM(NiceWrite_R8(_fact*fact*Int14Coul % Tot))
  str0=TRIM(str0)//ctemp//TRIM(NiceWrite_R8(_fact*Temp % Tot))
  str0=TRIM(str0)//blank9//blank12//pipe
  WRITE(kprint,'(a)') TRIM(str0)
  WRITE(kprint,'(a)') REPEAT('-',107)

END SUBROUTINE Write_it_
FUNCTION NiceWrite_R8(number) RESULT(out)
  CHARACTER(80) :: out
  CHARACTER(80) :: fmt
  REAL(8) :: number
  REAL(8), PARAMETER :: three=10.0D0**(12-1-3)-0.001D0,two=10.0D0**(12-1&
       &-2)-0.01D0,one=10.0D0**(12-1-1)-0.1D0,min=1.00D0
  
  REAL(8), PARAMETER :: threem=10.0D0**(11-1-3)-0.001D0,twom=10.0D0&
       &**(11-1-2)-0.01D0,onem=10.0D0**(11-1-1)-0.1D0,minm=1.00D0
  
  fmt='(e12.5)'
  IF(Number >= 0.0D0) THEN
     IF(Number <= three) fmt='(f12.3)'
     IF(Number > three .AND. Number <= two) fmt='(f12.2)'
     IF(Number > two .AND. Number <= one) fmt='(f12.1)'
     IF(Number > one) fmt='(e12.6)'
     IF(Number <= min) fmt='(e12.6)'
  ELSE
     IF(ABS(Number) <= threem) fmt='(f12.3)'
     IF(ABS(Number) > threem .AND. ABS(Number) <= twom) fmt='(f12.2)'
     IF(ABS(Number) > twom .AND. ABS(Number) <= onem) fmt='(f12.1)'
     IF(ABS(Number) > onem) fmt='(e12.5)'
     IF(ABS(Number) <= minm) fmt='(e12.5)'
  END IF
  out=' '
  WRITE(out(1:12),fmt) Number
  
END FUNCTION NiceWrite_R8
