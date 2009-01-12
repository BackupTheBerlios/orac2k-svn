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
MODULE Groups
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
  
  USE PrmUtilities
  USE SystemTpg
  USE SystemPrm
  USE AtomCnt
  USE IndSequence
  USE Atom
  USE PI_
  IMPLICIT none
  PRIVATE
  PUBLIC Groups_, Groups__Chain,  Groups__Base,Groupa&
       &,Groups__InitCoords,Groups__Update_Knwn, Groups__Update
  TYPE :: Groups__Base
     REAL(8) :: x,y,z      ! Coordinates orthogonal frame
     REAL(8) :: xa,ya,za   ! Coordinates reduced frame
     REAL(8) :: mass
     INTEGER ::  AtSt,AtEn ! Atom Start, Atom End
     INTEGER :: Res_No, Grp_No, Id_Res, Id_Type, Id_slv
     INTEGER :: Knwn=1
  END TYPE Groups__Base
  
  TYPE :: Groups__Chain
     INTEGER, ALLOCATABLE :: bd(:)
     INTEGER :: nbd
     REAL(8), ALLOCATABLE :: p(:)
     INTEGER :: np
  END TYPE Groups__Chain
  TYPE(Groups__Base), ALLOCATABLE, SAVE :: Groupa(:)
  INTEGER, SAVE, POINTER :: Res_Atm(:,:)=>NULL()
  INTEGER, SAVE, POINTER :: Grp_Atm(:,:)=>NULL()
CONTAINS
  FUNCTION Groups_() RESULT(out)
    LOGICAL :: out
    TYPE :: Local
       TYPE(Groups__Chain), ALLOCATABLE :: Tpg(:)
    END TYPE Local
    TYPE(local), POINTER :: Grp_Local(:)=>NULL()
    INTEGER :: n,s,m,o1,ntot
    
    out=.TRUE.
    Grp_Atm=>IndSequence__Grp()
  END FUNCTION Groups_
  FUNCTION Groups__InitCoords() RESULT(out)
    INTEGER :: ntap,ngrp,n,p1
    REAL(8) :: pmass
    REAL(8), ALLOCATABLE :: Tmass(:)
    LOGICAL :: out
    
    out=.TRUE.
    IF(.NOT. ASSOCIATED(Grp_Atm) .OR. .NOT. ALLOCATED(Atoms)) THEN
       out=.FALSE.
       RETURN
    END IF
    
    ngrp=SIZE(Grp_Atm,2)
    IF(.NOT. ALLOCATED(Groupa)) ALLOCATE(Groupa(ngrp))
    ntap=SIZE(Atoms)
    ALLOCATE(tmass(ngrp))
    tmass=0.0D0
    Groupa(:) % xa = 0.0D0
    Groupa(:) % ya = 0.0D0
    Groupa(:) % za = 0.0D0
    Groupa(:) % x = 0.0D0
    Groupa(:) % y = 0.0D0
    Groupa(:) % z = 0.0D0
    
    DO n=1,ntap
       p1=Atoms(n) % Grp_No
       pmass=Atoms(n) % pmass
       tmass(p1)=tmass(p1)+ Atoms(n) % mass
       Groupa(p1) % xa = Groupa(p1) % xa + pmass*Atoms(n) % xa
       Groupa(p1) % ya = Groupa(p1) % ya + pmass*Atoms(n) % ya
       Groupa(p1) % za = Groupa(p1) % za + pmass*Atoms(n) % za
       
       Groupa(p1) % x = Groupa(p1) % x + pmass*Atoms(n) % x
       Groupa(p1) % y = Groupa(p1) % y + pmass*Atoms(n) % y
       Groupa(p1) % z = Groupa(p1) % z + pmass*Atoms(n) % z
       Groupa(p1) % Res_No = Atoms(n) % Res_No
       Groupa(p1) % Grp_No = Atoms(n) % Grp_No
       Groupa(p1) % Id_Res = Atoms(n) % Id_Res
       Groupa(p1) % Id_Type = Atoms(n) % Id_Type
       Groupa(p1) % Id_Slv = Atoms(n) % Id_slv
    END DO
    DO n=1,ngrp
       Groupa(n) % mass = tmass(n)
    END DO
    CALL IndSequence__Update
    Grp_Atm=>IndSequence__Grp()
    DO n=1,ngrp
       Groupa(n) % AtSt = Grp_atm(1,n)
       Groupa(n) % AtEn = Grp_atm(2,n)
    END DO
  END FUNCTION Groups__InitCoords
  FUNCTION Groups__Update_Knwn() RESULT(out)
    INTEGER :: ntap,ngrp,n,AtSt,AtEn,nn
    REAL(8) :: pmass
    LOGICAL :: out
    
    
    out=.TRUE.
    IF(.NOT. ALLOCATED(Atoms)) THEN
       out=.FALSE.
       RETURN
    END IF
    
    ngrp=SIZE(groupa)
    ntap=SIZE(Atoms)
    
    DO n=1,ngrp
       IF(Groupa(n) % knwn == 2) THEN
          Groupa(n) % xa = 0.0D0
          Groupa(n) % ya = 0.0D0
          Groupa(n) % za = 0.0D0
          Groupa(n) % x = 0.0D0
          Groupa(n) % y = 0.0D0
          Groupa(n) % z = 0.0D0
       END IF
    END DO
    
    DO nn=1,ngrp
       IF(Groupa(nn) % knwn /= 2) CYCLE
       AtSt=Groupa(nn) % AtSt
       AtEn=Groupa(nn) % Aten
       DO n=AtSt,AtEn
          pmass=Atoms(n) % pmass
          Groupa(nn) % xa = Groupa(nn) % xa + pmass*Atoms(n) % xa
          Groupa(nn) % ya = Groupa(nn) % ya + pmass*Atoms(n) % ya
          Groupa(nn) % za = Groupa(nn) % za + pmass*Atoms(n) % za
          Groupa(nn) % x = Groupa(nn) % x + pmass*Atoms(n) % x
          Groupa(nn) % y = Groupa(nn) % y + pmass*Atoms(n) % y
          Groupa(nn) % z = Groupa(nn) % z + pmass*Atoms(n) % z
       END DO
    END DO
  END FUNCTION Groups__Update_Knwn
  FUNCTION Groups__Update() RESULT(out)
    INTEGER :: ntap,ngrp,n,AtSt,AtEn,nn
    REAL(8) :: pmass
    LOGICAL :: out
    
    
    out=.TRUE.
    IF(.NOT. ALLOCATED(Atoms)) THEN
       out=.FALSE.
       RETURN
    END IF
    
    ngrp=SIZE(groupa)
    ntap=SIZE(Atoms)
    
    Groupa(:) % x=0.0D0
    Groupa(:) % y=0.0D0
    Groupa(:) % z=0.0D0
    Groupa(:) % xa=0.0D0
    Groupa(:) % ya=0.0D0
    Groupa(:) % za=0.0D0

    DO nn=1,ngrp
       IF(Groupa(nn) % knwn == 0) CYCLE
       AtSt=Groupa(nn) % AtSt
       AtEn=Groupa(nn) % Aten
       DO n=AtSt,AtEn
          pmass=Atoms(n) % pmass
          Groupa(nn) % xa = Groupa(nn) % xa + pmass*Atoms(n) % xa
          Groupa(nn) % ya = Groupa(nn) % ya + pmass*Atoms(n) % ya
          Groupa(nn) % za = Groupa(nn) % za + pmass*Atoms(n) % za
          Groupa(nn) % x = Groupa(nn) % x + pmass*Atoms(n) % x
          Groupa(nn) % y = Groupa(nn) % y + pmass*Atoms(n) % y
          Groupa(nn) % z = Groupa(nn) % z + pmass*Atoms(n) % z
       END DO
    END DO
  END FUNCTION Groups__Update
END MODULE Groups
