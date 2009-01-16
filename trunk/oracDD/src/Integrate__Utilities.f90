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
!!$              - Mon Jan  5 2009 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*
FUNCTION Pick_Init(na,c) RESULT(out)
  INTEGER :: na,c,out
  
  INTEGER, SAVE :: First_Call=0
  INTEGER, ALLOCATABLE, SAVE :: Mult_Shell(:)
  INTEGER :: n
  IF(First_Call == 0) THEN
     ALLOCATE(Mult_Shell(NShell))
     DO n=1,NShell
        Mult_shell(n)=Get_Mult(n)
     END DO
     First_Call=First_Call+1
  END IF
  IF(MOD(c,Mult_shell(na)) == 0) THEN
     out=0
  ELSE
     out=1
  END IF
  IF(c == 0 .AND. na /= Nshell) out=1
CONTAINS
  RECURSIVE FUNCTION Get_Mult(n) RESULT(out)
    INTEGER :: n,out
    IF(n > NShell) THEN
       out=1
       RETURN
    END IF
    SELECT CASE(n)
    CASE(_N0_)
       out=n0_*Get_Mult(n+1)
    CASE(_N1_)
       out=n1_*Get_Mult(n+1)
    CASE(_M_)
       out=m_*Get_Mult(n+1)
    CASE(_L_)
       out=l_*Get_Mult(n+1)
    CASE(_H_)
       out=h_*Get_Mult(n+1)
    END SELECT
  END FUNCTION Get_Mult
END FUNCTION Pick_Init

FUNCTION Print_Now(n,Freq) RESULT(out)
  INTEGER :: n
  REAL(8) :: Freq
  LOGICAL :: out
  
  INTEGER :: My
  REAL(8) :: Time_at_step,dt,My_Freq
  
  SELECT CASE(n)
  CASE(_N0_)
     dt=dt_n0
  CASE(_N1_)
     dt=dt_n1
  CASE(_M_)
     dt=dt_m
  CASE(_L_)
     dt=dt_l
  CASE(_H_)
     dt=dt_h
  END SELECT
  
  Time_at_Step=Time_Step()
  
  My=INT(Freq/dt)
  IF(My == 0) My=1
  My_Freq=My*dt
  
  out=MOD(Time_at_Step,My_Freq) == 0
END FUNCTION Print_Now
FUNCTION  Get_RunLength() RESULT(out)
  TYPE(Length) :: out
  out % nstep=-1
  out % Time=-1.0_8
  
  SELECT CASE(NShell)
  CASE(_N0_)
     out % nstep=INT(Run_ % Time/dt_N0)
     out % Time=out % nstep*dt_n0
  CASE(_N1_)
     out % nstep=INT(Run_ % Time/dt_n1)
     out % Time=out % nstep*dt_n1
  CASE(_M_)
     out % nstep=INT(Run_ % Time/dt_m)
     out % Time=out % nstep*dt_m
  CASE(_L_)
     out % nstep=INT(Run_ % Time/dt_l)
     out % Time=out % nstep*dt_l
  CASE(_H_)
     out % nstep=INT(Run_ % Time/dt_h)
     out % Time=out % nstep*dt_h
  END SELECT
  WRITE(kprint,*) NShell,Run_ % Time,dt_H
END FUNCTION Get_RunLength
