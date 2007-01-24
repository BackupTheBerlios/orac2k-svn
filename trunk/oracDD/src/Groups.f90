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
  IMPLICIT none
  PRIVATE
  PUBLIC Groups_, Groups__Chain, Groups__Pot, Grp_A,Grp_L
  TYPE :: Groups__Base
     REAL(8) :: x,y,z      ! Coordinates
     INTEGER ::  AtSt,AtEn ! Atom Start, Atom End
     INTEGER :: Res_No, Grp_No, Id_Res, Id_Type, Id_slv
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
  TYPE(Groups__Pot), ALLOCATABLE :: Grp_L(:), Grp_A(:)
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
    ALLOCATE(Grp_A(SIZE(Grp_Atm)))

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
      ALLOCATE(dims_a(SIZE(Grp_Atm)))
      IF(ALLOCATED(Grp)) DEALLOCATE(Grp)
      ALLOCATE(Grp(SIZE(Grp_Atm)))

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

END MODULE Groups
