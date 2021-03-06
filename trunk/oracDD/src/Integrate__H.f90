!!$/---------------------------------------------------------------------\
!!$   Copyright  � 2006-2007 Massimo Marchi <Massimo.Marchi at cea.fr>   |
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

SUBROUTINE Integrate_h
  INTEGER, SAVE :: counter=0,MyCalls=0
  INTEGER, PARAMETER :: Init=1
  Real(8), Save :: Startts,Endtts,timets
  Integer :: ha,n
  Real(8) :: tfact

  DO ha=1,h_
     startts=MPI_WTIME()

     __correct_vp(dt_h,fpp_h)

     IF(.NOT. Rattle_it(dt_h,RATTLE__Correct_)) CALL Print_Errors()
     
     CALL Integrate_l
     
     counter=counter+1
     Call SaveCounter(_H_,counter)

     CALL FORCES_Zero(_H_)
     CALL Forces_(_H_)

     Call GatherLocals(fpp_h,fp_h(:)%x,fp_h(:)%y,fp_h(:)%z,IndBox_a_p)

     __correct_vp(dt_h,fpp_h)

     IF(NShell == _H_) THEN
        IF(.NOT. RATTLE__Parameters_(Atoms(:) % mass,Atoms(:) %&
             & knwn)) CALL Print_Errors()
     END IF
     MyCalls=MyCalls+1
     endtts=MPI_WTIME()
     Timets=Timets+endtts-startts
     Write(kprint,'('' H Mycalls '',i5,'' CPU Time '',2e14.6)') MyCalls,timets/Dble(MyCalls),endtts-startts
  END DO
END SUBROUTINE Integrate_h
