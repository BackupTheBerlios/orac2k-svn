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
MODULE Energies
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Jan 22 2009 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*

#ifdef HAVE_MPI
  USE mpi
#endif
  USE PI_
  USE Print_defs
  USE STRPAK, ONLY: SP_Putnum
  USE CONSTANTS
  USE Units
  IMPLICIT none
  PRIVATE
  PUBLIC :: Lj_,Coul_Dir_,Coul_Rec_,Stretch_,Angle_,Imph_,Int14LJ_&
       &,Int14Coul_,Dihed_,Write_it_,Get_energies_,Kinetic_,Temp_,Banner_,Total_
#include "Energies.h"
  Type(Energy_Type), SAVE :: Total,TotPot,Lj(3),Coul_Dir(3),Coul_Rec,Bond&
       &,Stretch,Angle,Imph,Int14LJ,Int14Coul,Dihed,Kinetic,Temp
  Type(Energy_Type), SAVE :: Total_a,TotPot_a,Lj_a(3),Coul_Dir_a(3),Coul_Rec_a&
       &,Bond_a,Stretch_a,Angle_a,Imph_a,Int14LJ_a,Int14Coul_a&
       &,Dihed_a,Kinetic_a,temp_a 
  REAL(8) :: en(4)
  INTERFACE Reduce_it_
     MODULE PROCEDURE Reduce_A_it_
     MODULE PROCEDURE Reduce_S_it_
  END INTERFACE
CONTAINS
  SUBROUTINE Reduce_A_it_(Energy)
    Type(Energy_Type) :: Energy(:)
    INTEGER :: nsize,msize,n,ierr
    REAL(8), ALLOCATABLE :: eny(:)

    nsize=SIZE(Energy)
    msize=nsize*4
    ALLOCATE(eny(msize))

    DO n=1,nsize
       eny(4*(n-1)+1)=Energy(n) % Tot
       eny(4*(n-1)+2)=Energy(n) % Slv
       eny(4*(n-1)+3)=Energy(n) % Slt
       eny(4*(n-1)+4)=Energy(n) % Mix
    END DO

#ifdef HAVE_MPI
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,eny,msize,MPI_REAL8,MPI_SUM,PI_Comm,ierr)
#endif

    DO n=1,nsize
       Energy(n) % Tot=eny(4*(n-1)+1)/Dble(PI_Nprocs)
       Energy(n) % Slv=eny(4*(n-1)+2)/Dble(PI_Nprocs)
       Energy(n) % Slt=eny(4*(n-1)+3)/Dble(PI_Nprocs)
       Energy(n) % Mix=eny(4*(n-1)+4)/Dble(PI_Nprocs)
    END DO

  END SUBROUTINE Reduce_A_it_
  SUBROUTINE Reduce_S_it_(Energy)
    Type(Energy_Type) :: Energy
    INTEGER :: msize,n,ierr
    REAL(8), ALLOCATABLE :: eny(:)


    msize=4
    ALLOCATE(eny(msize))

    eny(1)=Energy % Tot
    eny(2)=Energy % Slv
    eny(3)=Energy % Slt
    eny(4)=Energy % Mix
#ifdef HAVE_MPI
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,eny,msize,MPI_REAL8,MPI_SUM,PI_Comm,ierr)
#endif

    Energy % Tot=eny(1)/Dble(PI_Nprocs)
    Energy % Slv=eny(2)/Dble(PI_Nprocs)
    Energy % Slt=eny(3)/Dble(PI_Nprocs)
    Energy % Mix=eny(4)/Dble(PI_Nprocs)
  END SUBROUTINE Reduce_S_it_
  SUBROUTINE Total_
    CALL TotPot_

    Total % tot=TotPot % tot + Kinetic % Tot
    Total % Slv=TotPot % Slv + Kinetic % Slv
    Total % Slt=TotPot % Slt + Kinetic % Slt
    Total % Mix=TotPot % Mix + Kinetic % Mix

    Total_a % tot=Total_a % tot+Total % tot
    Total_a % Slv=Total_a % Slv+Total % Slv
    Total_a % Slt=Total_a % Slt+Total % Slt
    Total_a % Mix=Total_a % Mix+Total % Mix
    Total_a % n0=Total_a % n0+1
  END SUBROUTINE Total_
  SUBROUTINE TotPot_
    Bond % Slv = Stretch % Slv + Angle % Slv + Imph % Slv + Int14LJ % Slv&
         & + Int14Coul % Slv + Dihed % Slv 
    Bond % Slt = Stretch % Slt + Angle % Slt + Imph % Slv + Int14LJ % Slt&
         & + Int14Coul % Slt + Dihed % Slt 
    Bond % Mix = Stretch % Mix + Angle % Mix + Imph % Mix + Int14LJ % Mix&
         & + Int14Coul % Mix + Dihed % Mix 
    Bond % Tot = Stretch % Tot + Angle % Tot + Imph % Tot + Int14LJ % Tot&
         & + Int14Coul % Tot + Dihed % Tot 

    TotPot % Slv=SUM(Lj % Slv) + SUM(Coul_Dir % Slv) + Coul_Rec % Slv + Bond % Slv
    TotPot % Slt=SUM(Lj % Slt) + SUM(Coul_Dir % Slt) + Coul_Rec % Slt + Bond % Slt
    TotPot % Mix=SUM(Lj % Mix) + SUM(Coul_Dir % Mix) + Coul_Rec % Mix + Bond % Mix
    TotPot % Tot=SUM(Lj % Tot) + SUM(Coul_Dir % Tot) + Coul_Rec % Tot + Bond % Tot

    TotPot_a % Slv=TotPot_a % Slv+TotPot % Slv
    TotPot_a % Slt=TotPot_a % Slt+TotPot % Slt
    TotPot_a % Mix=TotPot_a % Mix+TotPot % Mix
    TotPot_a % Tot=TotPot_a % Tot+TotPot % Tot
    TotPot_a % n0=TotPot_a % n0+1
  END SUBROUTINE TotPot_
  SUBROUTINE Coul_Rec_(Mix)
    REAL(8) :: Mix
    REAL(8) :: tot
    INTEGER :: ierr

    tot=Mix
    Coul_Rec%Tot=Tot
    Coul_Rec%Mix=Mix

  END SUBROUTINE Coul_Rec_
  SUBROUTINE Coul_Dir_(Id,Slv,Slt,Mix)
    INTEGER :: Id
    REAL(8) :: Slv,Slt,Mix
    REAL(8) :: tot
    INTEGER :: ierr

    tot=Slv+Slt+Mix
    en(1)=tot; en(2)=Slv;en(3)=Slt;en(4)=Mix

    Coul_Dir(Id)%Tot=en(1)
    Coul_Dir(Id)%Slv=en(2)
    Coul_Dir(Id)%Slt=en(3)
    Coul_Dir(Id)%Mix=en(4)

  END SUBROUTINE Coul_Dir_
  SUBROUTINE LJ_(Id,Slv,Slt,Mix)
    INTEGER :: Id
    REAL(8) :: Slv,Slt,Mix
    REAL(8) :: tot
    INTEGER :: ierr

    tot=Slv+Slt+Mix
    en(1)=tot; en(2)=Slv;en(3)=Slt;en(4)=Mix

    LJ(Id)%Tot=en(1)
    LJ(Id)%Slv=en(2)
    LJ(Id)%Slt=en(3)
    LJ(Id)%Mix=en(4)

  END SUBROUTINE LJ_
  SUBROUTINE Int14LJ_(Slv,Slt,Mix)
    REAL(8) :: Slv,Slt,Mix
    REAL(8) :: tot
    INTEGER :: ierr

    tot=Slv+Slt+Mix
    en(1)=tot; en(2)=Slv;en(3)=Slt;en(4)=Mix

    Int14LJ%Tot=en(1)
    Int14LJ%Slv=en(2)
    Int14LJ%Slt=en(3)
    Int14LJ%Mix=en(4)

  END SUBROUTINE Int14LJ_
  SUBROUTINE Int14Coul_(Slv,Slt,Mix)
    REAL(8) :: Slv,Slt,Mix
    REAL(8) :: tot
    INTEGER :: ierr

    tot=Slv+Slt+Mix
    en(1)=tot; en(2)=Slv;en(3)=Slt;en(4)=Mix

    Int14Coul%Tot=en(1)
    Int14Coul%Slv=en(2)
    Int14Coul%Slt=en(3)
    Int14Coul%Mix=en(4)

  END SUBROUTINE Int14Coul_
  SUBROUTINE Imph_(Slv,Slt,Mix)
    REAL(8) :: Slv,Slt,Mix
    REAL(8) :: tot
    INTEGER :: ierr

    tot=Slv+Slt+Mix
    en(1)=tot; en(2)=Slv;en(3)=Slt;en(4)=Mix

    Imph%Tot=en(1)
    Imph%Slv=en(2)
    Imph%Slt=en(3)
    Imph%Mix=en(4)

  END SUBROUTINE Imph_
  SUBROUTINE Dihed_(Slv,Slt,Mix)
    REAL(8) :: Slv,Slt,Mix
    REAL(8) :: tot
    INTEGER :: ierr

    tot=Slv+Slt+Mix
    en(1)=tot; en(2)=Slv;en(3)=Slt;en(4)=Mix

    Dihed%Tot=en(1)
    Dihed%Slv=en(2)
    Dihed%Slt=en(3)
    Dihed%Mix=en(4)

  END SUBROUTINE Dihed_
  SUBROUTINE Angle_(Slv,Slt,Mix)
    REAL(8) :: Slv,Slt,Mix
    REAL(8) :: tot
    INTEGER :: ierr

    tot=Slv+Slt+Mix
    en(1)=tot; en(2)=Slv;en(3)=Slt;en(4)=Mix

    Angle%Tot=en(1)
    Angle%Slv=en(2)
    Angle%Slt=en(3)
    Angle%Mix=en(4)

  END SUBROUTINE Angle_
  SUBROUTINE Stretch_(Slv,Slt,Mix)
    REAL(8) :: Slv,Slt,Mix
    REAL(8) :: tot
    INTEGER :: ierr

    tot=Slv+Slt+Mix
    en(1)=tot; en(2)=Slv;en(3)=Slt;en(4)=Mix

    Stretch%Tot=en(1)
    Stretch%Slv=en(2)
    Stretch%Slt=en(3)
    Stretch%Mix=en(4)

  END SUBROUTINE Stretch_
  SUBROUTINE Kinetic_(Slv,Slt)
    REAL(8) :: Slv,Slt
    Kinetic%tot=Slv+Slt
    Kinetic%Slv=Slv
    Kinetic%Slt=Slt

  END SUBROUTINE Kinetic_
  SUBROUTINE Temp_(Tot,Slv,Slt)
    REAL(8) :: Slv,Slt,Tot    
    Temp%tot=Tot
    Temp%Slv=Slv
    Temp%Slt=Slt
  END SUBROUTINE Temp_
  SUBROUTINE Banner_
    WRITE(kprint,100) 
100 FORMAT(&
&'-----------------------------------------------------------------------------------------'/&
&'|                        --- Output    C  O  D  E  S  ---                               |'/&
&'|                                      ....                                             |'/&
&'|      Tot  = Total            ; TPot    = Total Potential    ; LJ      =  LJ Energy  ; |'/&
&'|      CDir = Direct Coulombic ; CRec    = Reciprocal Coul    ; Stre    = Stretching  ; |'/&
&'|     Angle = Angle bending    ; Imph    = Improper Torsions  ; 14LJ    = 14 LJ       ; |'/&
&'|      14Co = 14 Coulombic     ; Dihed   = Proper Torsions    ; Bnd     = Total Bonded; |'/&
&'|      Kin  = Kinetic          ; Temp    = System Temperature ; Time    = Time in fs  ; |'/&
&'|                                                                                       |'/&
&'-----------------------------------------------------------------------------------------'/)
  END SUBROUTINE Banner_
  SUBROUTINE Get_Energies_
    CALL Reduce_it_(Stretch)
    __avg(Stretch_a,Stretch)

    CALL Reduce_it_(Angle)
    __avg(Angle_a,Angle)

    CALL Reduce_it_(Dihed)
    __avg(Dihed_a,Dihed)


    CALL Reduce_it_(Imph)
    __avg(Imph_a,Imph)

    CALL Reduce_it_(Int14LJ)
    __avg(Int14LJ_a,Int14LJ)

    CALL Reduce_it_(Int14Coul)
    __avg(Int14Coul_a,Int14Coul)

    CALL Reduce_it_(Lj)
    __avg(Lj_a,Lj)
    CALL Reduce_it_(Coul_Dir)
    __avg(Coul_Dir_a,Coul_Dir)

    CALL Reduce_it_(Coul_Rec)
    __avg(Coul_Rec_a,Coul_Rec)

    CALL Reduce_it_(Kinetic)
    __avg(Kinetic_a,Kinetic)

    CALL Reduce_it_(Temp)

    __avg(Temp_a,Temp)

    CALL Reduce_it_(Bond)
    __avg(Bond_a,Bond)

    CALL Reduce_it_(TotPot)
    __avg(TotPot_a,TotPot)

    CALL Reduce_it_(Total)
    __avg(Total_a,Total)

  END SUBROUTINE Get_Energies_
#include "ENERGIES__Write.f90"
END MODULE Energies
