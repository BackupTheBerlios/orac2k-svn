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
  PUBLIC Groups_, Groups__Chain, Groups__Pot, Grp_A, Groups__Base,Groupa&
       &, Groups__InitCoords,Groups__Update_Knwn
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

  TYPE :: Groups__Pot
     TYPE(Groups__Chain), ALLOCATABLE :: bonds(:),angles(:),dihed(:)&
          &,imph(:),int14(:)
  END TYPE Groups__Pot
  TYPE :: Groups__Param
     INTEGER, ALLOCATABLE :: bonds(:),angles(:),dihed(:)&
          &,imph(:),int14(:)
  END type Groups__Param
  TYPE(Groups__Pot), ALLOCATABLE, SAVE :: Grp_A(:)
  TYPE(Groups__Param), ALLOCATABLE, SAVE :: Grp_B(:)

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
    ALLOCATE(Grp_A(SIZE(Grp_Atm,2)))
    Grp_Local=>Do_Copy(Prm % bonds, Tpg % bonds)

    DO n=1,SIZE(Grp_Local)
       s=0
       IF(ALLOCATED(Grp_Local(n) % Tpg)) s=SIZE(Grp_Local(n) % Tpg)
       IF(s == 0) CYCLE
       ALLOCATE(Grp_A(n) % bonds(s))
       DO m=1,s
          o1=SIZE(Grp_Local(n) % Tpg (m) % bd)
          ALLOCATE(Grp_A(n) % bonds (m) % bd(o1))
          Grp_A(n) % bonds (m) % nbd = o1
          o1=SIZE(Grp_Local(n) % Tpg (m) % p)
          ALLOCATE(Grp_A(n) % bonds (m) % p(o1))
          Grp_A(n) % bonds (m) % np = o1
          Grp_A(n) % bonds (m) % bd = Grp_Local(n) % Tpg (m) % bd
          Grp_A(n) % bonds (m) % p  = Grp_Local(n) % Tpg (m) % p
       END DO

       DO m=1,s
          o1=Grp_A(n) % bonds(m) % nbd
       END DO
    END DO
    Grp_Local=>Do_Copy(Prm % angles, Tpg % angles)
    DO n=1,SIZE(Grp_Local)
       s=0
       IF(ALLOCATED(Grp_Local(n) % Tpg)) s=SIZE(Grp_Local(n) % Tpg)
       IF(s == 0) CYCLE
       ALLOCATE(Grp_A(n) % angles(s))
       DO m=1,s
          o1=SIZE(Grp_Local(n) % Tpg (m) % bd)
          ALLOCATE(Grp_A(n) % angles (m) % bd(o1))
          Grp_A(n) % angles (m) % nbd = o1
          o1=SIZE(Grp_Local(n) % Tpg (m) % p)
          ALLOCATE(Grp_A(n) % angles (m) % p(o1))
          Grp_A(n) % angles (m) % np = o1
          Grp_A(n) % angles (m) % bd = Grp_Local(n) % Tpg (m) % bd
          Grp_A(n) % angles (m) % p  = Grp_Local(n) % Tpg (m) % p
       END DO
    END DO
    Grp_Local=>Do_Copy(Prm % dihed, Tpg % dihed)
    DO n=1,SIZE(Grp_Local)
       s=0
       IF(ALLOCATED(Grp_Local(n) % Tpg)) s=SIZE(Grp_Local(n) % Tpg)
       IF(s == 0) CYCLE
       ALLOCATE(Grp_A(n) % dihed(s))
       DO m=1,s
          o1=SIZE(Grp_Local(n) % Tpg (m) % bd)
          ALLOCATE(Grp_A(n) % dihed (m) % bd(o1))
          Grp_A(n) % dihed (m) % nbd = o1
          o1=SIZE(Grp_Local(n) % Tpg (m) % p)
          ALLOCATE(Grp_A(n) % dihed (m) % p(o1))
          Grp_A(n) % dihed (m) % np = o1
          Grp_A(n) % dihed (m) % bd = Grp_Local(n) % Tpg (m) % bd
          Grp_A(n) % dihed (m) % p  = Grp_Local(n) % Tpg (m) % p
       END DO
    END DO
    Grp_Local=>Do_Copy(Prm % imph, Tpg % imph)
    DO n=1,SIZE(Grp_Local)
       s=0
       IF(ALLOCATED(Grp_Local(n) % Tpg)) s=SIZE(Grp_Local(n) % Tpg)
       IF(s == 0) CYCLE
       ALLOCATE(Grp_A(n) % imph(s))
       DO m=1,s
          o1=SIZE(Grp_Local(n) % Tpg (m) % bd)
          ALLOCATE(Grp_A(n) % imph (m) % bd(o1))
          Grp_A(n) % imph (m) % nbd = o1
          o1=SIZE(Grp_Local(n) % Tpg (m) % p)
          ALLOCATE(Grp_A(n) % imph (m) % p(o1))
          Grp_A(n) % imph (m) % np = o1
          Grp_A(n) % imph (m) % bd = Grp_Local(n) % Tpg (m) % bd
          Grp_A(n) % imph (m) % p  = Grp_Local(n) % Tpg (m) % p
       END DO
    END DO
    Grp_Local=>Do_Copy(Prm % bonds, Tpg % bonds, 0)
  CONTAINS
    FUNCTION Do_Copy(p_Tpg, t_Tpg, delete) RESULT(out)
      TYPE(SystemPrm__Chain) :: p_tpg(:)
      INTEGER :: t_tpg(:,:)
      INTEGER, OPTIONAL :: delete
      TYPE(Local), POINTER :: out(:)

      INTEGER, ALLOCATABLE, SAVE :: dims_a(:)
      TYPE(Local), ALLOCATABLE, SAVE, TARGET :: Grp(:)
      INTEGER :: n,p1,G_p1,o1

      out=>NULL()
      IF(PRESENT(delete)) THEN
         IF(ALLOCATED(dims_a)) DEALLOCATE(dims_a)
         IF(ALLOCATED(Grp)) DEALLOCATE(Grp)
         RETURN
      END IF

      IF(ALLOCATED(dims_a)) DEALLOCATE(dims_a)
      ALLOCATE(dims_a(SIZE(Grp_Atm,2)))
      IF(ALLOCATED(Grp)) DEALLOCATE(Grp)
      ALLOCATE(Grp(SIZE(Grp_Atm,2)))

      Dims_a=0
      DO n=1,SIZE(P_Tpg)
         o1=p_Tpg(n)  % pt
         p1=t_Tpg (1,o1)
         G_p1=AtomCnts (p1) % Grp_No
         Dims_a(G_p1) = Dims_a(G_p1)  +1
      END DO
      DO n=1,SIZE(Dims_a)
         ALLOCATE(Grp(n) % Tpg (Dims_a(n)))
      END DO

      Dims_a=0
      DO n=1,SIZE(P_Tpg)
         o1=p_Tpg(n)  % pt
         p1=t_Tpg (1,o1)
         G_p1=AtomCnts (p1) % Grp_No
         Dims_a(G_p1) = Dims_a(G_p1)  + 1
         m=Dims_a(G_p1)
         ALLOCATE(Grp(G_p1) % Tpg(m) % bd (SIZE(t_Tpg,1)))
         ALLOCATE(Grp(G_p1) % Tpg(m) % p (SIZE(p_Tpg (n) % g)))
         
         Grp(G_p1) % Tpg(m) % bd = t_Tpg(:,o1)
         Grp(G_p1) % Tpg(m) % p = p_Tpg (n) % g
         Grp(G_p1) % Tpg(m) % nbd = SIZE(t_Tpg,1)
         Grp(G_p1) % Tpg(m) % np = SIZE(p_Tpg (n) % g)
      END DO
      out=>Grp
    END FUNCTION Do_Copy
      
  END FUNCTION Groups_
  FUNCTION Groups__InitCoords() RESULT(out)
    INTEGER :: ntap,ngrp,n,p1
    REAL(8) :: mass
    REAL(8), ALLOCATABLE :: Tmass(:)
    LOGICAL :: out

    out=.TRUE.
    IF(.NOT. ALLOCATED(Grp_A) .OR. .NOT. ALLOCATED(Atoms)) THEN
       out=.FALSE.
       RETURN
    END IF

    ngrp=SIZE(Grp_A)
    ALLOCATE(Groupa(ngrp))
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
       mass=Atoms(n) % mass
       tmass(p1)=tmass(p1)+ mass
       Groupa(p1) % xa = Groupa(p1) % xa + mass*Atoms(n) % xa
       Groupa(p1) % ya = Groupa(p1) % ya + mass*Atoms(n) % ya
       Groupa(p1) % za = Groupa(p1) % za + mass*Atoms(n) % za

       Groupa(p1) % x = Groupa(p1) % x + mass*Atoms(n) % x
       Groupa(p1) % y = Groupa(p1) % y + mass*Atoms(n) % y
       Groupa(p1) % z = Groupa(p1) % z + mass*Atoms(n) % z
       Groupa(p1) % Res_No = Atoms(n) % Res_No
       Groupa(p1) % Grp_No = Atoms(n) % Grp_No
       Groupa(p1) % Id_Res = Atoms(n) % Id_Res
       Groupa(p1) % Id_Type = Atoms(n) % Id_Type
       Groupa(p1) % Id_Slv = Atoms(n) % Id_slv
    END DO
    DO n=1,ngrp
       Groupa(n)%xa=Groupa(n)%xa/tmass(n)
       Groupa(n)%ya=Groupa(n)%ya/tmass(n)
       Groupa(n)%za=Groupa(n)%za/tmass(n)
       
       Groupa(n)%x=Groupa(n)%x/tmass(n)
       Groupa(n)%y=Groupa(n)%y/tmass(n)
       Groupa(n)%z=Groupa(n)%z/tmass(n)
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
    REAL(8) :: mass
    REAL(8) :: Tmass
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
          mass=Atoms(n) % mass
          Groupa(nn) % xa = Groupa(nn) % xa + mass*Atoms(n) % xa
          Groupa(nn) % ya = Groupa(nn) % ya + mass*Atoms(n) % ya
          Groupa(nn) % za = Groupa(nn) % za + mass*Atoms(n) % za
          Groupa(nn) % x = Groupa(nn) % x + mass*Atoms(n) % x
          Groupa(nn) % y = Groupa(nn) % y + mass*Atoms(n) % y
          Groupa(nn) % z = Groupa(nn) % z + mass*Atoms(n) % z
       END DO
    END DO
    DO n=1,ngrp
       IF(Groupa(n) % knwn /= 2) CYCLE
       tmass=Groupa(n) % Mass
       Groupa(n)%xa=Groupa(n)%xa/tmass
       Groupa(n)%ya=Groupa(n)%ya/tmass
       Groupa(n)%za=Groupa(n)%za/tmass
       
       Groupa(n)%x=Groupa(n)%x/tmass
       Groupa(n)%y=Groupa(n)%y/tmass
       Groupa(n)%z=Groupa(n)%z/tmass
    END DO
  END FUNCTION Groups__Update_Knwn
END MODULE Groups
