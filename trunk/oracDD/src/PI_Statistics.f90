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
MODULE PI_Statistics
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Nov  7 2008 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*

#ifdef HAVE_MPI
  USE mpi
#endif
  USE PI_
  IMPLICIT none
  PRIVATE
  PUBLIC Write_It, Time_It, Sample_Exchange, Add_Calls, Calls
  TYPE :: Statistics
     REAL(8) :: KByte_S=0.0_8,KByte_R=0.0_8
     REAL(8) :: Atoms_S=0.0_8,Atoms_R=0.0_8
  END type Statistics
  TYPE :: Times
     REAL(8) :: Tot,Comms
  END type Times
  TYPE(Statistics), SAVE :: Comms
  TYPE(Times), SAVE :: Timea
  INTEGER, SAVE :: Calls=0
CONTAINS
  SUBROUTINE Add_Calls
    Calls=Calls+1
  END SUBROUTINE Add_Calls
  SUBROUTINE Sample_Exchange(NoAtm_S,NoAtm_R)
    INTEGER :: NoAtm_S, NoAtm_R

    Comms % Atoms_S=Comms % Atoms_S+DBLE(NoAtm_s)
    Comms % Atoms_R=Comms % Atoms_R+DBLE(NoAtm_r)
    Comms % KByte_S=Comms % KByte_S+DBLE(NoAtm_s*3*8)/1024.0_8
    Comms % KByte_R=Comms % KByte_R+DBLE(NoAtm_r*3*8)/1024.0_8
  END SUBROUTINE Sample_Exchange
  SUBROUTINE Time_It(startime,endtime)
    REAL(8) :: out ,startime,endtime
    timea % Comms=timea % Comms+endtime-startime
  END SUBROUTINE Time_It
  FUNCTION Write_It RESULT(out)
    LOGICAL :: out
    REAL(8) :: Kbyte_s, Kbyte_r, Atoms_s, Atoms_r
    REAL(8) :: Kbyte_s0, Kbyte_r0, Atoms_s0, Atoms_r0
    REAL(8) :: Timeb,Timeb0,Timec,timec0

    out=.TRUE.
    IF(PI_Nprocs == 1) RETURN
#ifdef HAVE_MPI
    Kbyte_s0=Comms % KByte_s/DBLE(Calls)
    Kbyte_r0=Comms % KByte_r/DBLE(Calls)
    Atoms_s0=Comms % Atoms_s/DBLE(Calls)
    Atoms_r0=Comms % Atoms_r/DBLE(Calls)
    CALL MPI_ALLREDUCE(Kbyte_s0,Kbyte_s,1,MPI_REAL8,MPI_SUM,PI_Comm_Cart,ierr)
    CALL MPI_ALLREDUCE(Kbyte_r0,Kbyte_r,1,MPI_REAL8,MPI_SUM,PI_Comm_Cart,ierr)
    CALL MPI_ALLREDUCE(Atoms_s0,Atoms_s,1,MPI_REAL8,MPI_SUM,PI_Comm_Cart,ierr)
    CALL MPI_ALLREDUCE(Atoms_r0,Atoms_r,1,MPI_REAL8,MPI_SUM,PI_Comm_Cart,ierr)
    CALL MPI_ALLREDUCE(Timea % Comms,Timeb,1,MPI_REAL8,MPI_SUM,PI_Comm_Cart,ierr)
    CALL MPI_ALLREDUCE(Timea % Tot,Timec,1,MPI_REAL8,MPI_SUM,PI_Comm_Cart,ierr)

    IF(PI_Node_Cart == 0) THEN
       Kbyte_s=Kbyte_s/DBLE(PI_Nprocs)
       Kbyte_r=Kbyte_r/DBLE(PI_Nprocs)
       Atoms_s=Atoms_s/DBLE(PI_Nprocs)
       Atoms_r=Atoms_r/DBLE(PI_Nprocs)
       Timeb=Timeb/DBLE(PI_Nprocs)
       Timec=Timec/DBLE(PI_Nprocs)
       WRITE(*,100) Atoms_s,Kbyte_s,Atoms_r,Kbyte_r
       WRITE(*,200) Timeb,Timec,Timec-Timeb
    END IF
100 FORMAT(/'=====>    Average data transfer by each CPU per full shift    <====='&
         &/'=====>    ',f12.2,' atoms (',f12.4,' KB ) sent          <====='&
         &/'=====>    ',f12.2,' atoms (',f12.4,' KB ) received      <=&
         &===='/)
200 FORMAT(/'=====>        Timing                   <====='&
          &/'=====>   Comm. Time = ',f12.5,' s        <====='/&
          &/'=====>   Tot.  Time = ',f12.5,' s        <====='/&
          &/'=====>   Rem.  Time = ',f12.5,' s        <====='/&
          &)
#endif
  END FUNCTION Write_It

END MODULE PI_Statistics
