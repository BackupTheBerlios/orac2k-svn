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

SUBROUTINE Integrate_m
  INTEGER, SAVE :: counter=0
  INTEGER :: Ind,p,o,ma,n1
  INTEGER, PARAMETER :: init=1
  INTEGER :: MyCalls=0
  Real(8), Save :: Startts,Endtts,timets
  Integer :: n
  Real(8) :: tfact

  DO ma=1,m_

     __correct_vp(dt_m,fpp_m)

     IF(.NOT. Rattle_it(dt_m,RATTLE__Correct_)) CALL Print_Errors()
     
     CALL Integrate_n1

!!$
!!$--- From Cartesian get reduced coordinate then recompute group coords
!!$

     counter=counter+1
     Call SaveCounter(_M_,counter)

     IF(.NOT. Atom__Convert(_X_TO_XA_,IndBox_a_p(:))) CALL Print_Errors()
     IF(.NOT. Groups__Update(IndBox_g_p)) CALL Print_Errors()
     
     CALL FORCES_Zero(_M_)
     CALL Forces_(_M_)

     Call GatherLocals(fpp_m,fp_m(:)%x,fp_m(:)%y,fp_m(:)%z,IndBox_a_p)
     __correct_vp(dt_m,fpp_m)

     IF(NShell == _M_) THEN
        IF(.NOT. RATTLE__Parameters_(Atoms(:) % mass,Atoms(:) %&
             & knwn)) CALL Print_Errors()
     END IF
     Call GatherGlobal(atoms(:)%vx,atoms(:)%vy,atoms(:)%vz,v0,IndBox_a_p)
     IF(MOD(counter,Time_of_Print % nstep) == 0) THEN
        CALL ATOM__KinCompute
        CALL EN_Total_
        CALL EN_Write_it_(dt_n0*DBLE(nstep_n0))
     END IF
     
     IF(Inout__PDB % Unit /= 0 .AND. Print_Now(_M_,Inout__PDB % Freq)) THEN
        CALL Atom__PDB(Inout__PDB % Unit)
     END IF
  END DO
END SUBROUTINE Integrate_m
