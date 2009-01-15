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
MODULE AddHydrogens_

!!$***********************************************************************
!!$   Time-stamp: <2007-01-14 17:36:33 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Sun Jan 14 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program oracDD ----*

!!$-----------------------------------------------------------------------\
!!$                                                                       |
!!$--- f95 wrapper for old f77 routines, to add hydrogens to a .pdb       |
!!$--- structure                                                          |
!!$                                                                       |
!!$-----------------------------------------------------------------------/


!!$======================== DECLARATIONS ================================*

  USE Errors, ONLY: Add_Errors=>Add, Print_Errors
  USE AtomCnt
  USE PDB
  USE Print_Defs
  IMPLICIT none
  PRIVATE
  PUBLIC AddHydrogens__
CONTAINS
  FUNCTION AddHydrogens__(PDB__Coords) RESULT(out)
    LOGICAL :: out
    TYPE(AtomPdb), POINTER :: PDB__Coords(:)
    INTEGER :: n,m,p,o,q,nh,nnh,nna,nPDB__Coords,iret
    REAL(8), ALLOCATABLE :: xh(:),yh(:),zh(:),x0(:),y0(:),z0(:)
    INTEGER, ALLOCATABLE :: blist(:), ind_x(:), ind_y(:)
    LOGICAL, ALLOCATABLE :: ok_Hyd(:)
    CHARACTER(len=80) :: errmsg
    LOGICAL :: ok

    out=.TRUE.
    nPDB__Coords=SIZE(PDB__Coords)
    ALLOCATE(ind_x(SIZE(AtomCnts)))
    ind_x=0
    DO n=1,nPDB__Coords
       ind_x(PDB__Coords(n) % Serial) = n
    END DO

    ALLOCATE(ok_Hyd(nPDB__Coords))
    ok_Hyd=.TRUE.
    WHERE(PDB__Coords(:) % x < 1.0D8) ok_Hyd=.FALSE.

    iret=0
    errmsg=' '
    DO n=1,nPDB__Coords
       IF(.NOT. ok_Hyd(n)) THEN
          p=PDB__Coords(n) % Serial
          ok=.FALSE.
          nh=0
          nnh=0
          DO m=1,SIZE(AtomCnts(p) % cnt)
             o=ind_x(AtomCnts(p) % cnt(m))
             IF(PDB__Coords(o) % x < 1.0D8) THEN
                nnh=nnh+1
             ELSE
                ok=.TRUE.
                nh=nh+1
             END IF
          END DO
          IF(ok) THEN
             ALLOCATE(xh(nnh+1),yh(nnh+1),zh(nnh+1))
             ALLOCATE(x0(nh),y0(nh),z0(nh))
             ALLOCATE(blist(nnh))
             ALLOCATE(ind_y(nh))
             xh=0.0D0;yh=0.0D0;zh=0.0D0 
             x0=0.0D0;y0=0.0D0;z0=0.0D0 
             xh(1)=PDB__Coords(n) % x
             yh(1)=PDB__Coords(n) % y
             zh(1)=PDB__Coords(n) % z
             nh=0
             nnh=0
             ind_y=0
             DO m=1,SIZE(AtomCnts(p) % cnt)
                o=ind_x(AtomCnts(p) % cnt(m))
                IF(PDB__Coords(o) % x < 1.0D8) THEN
                   nnh=nnh+1
                   xh(nnh+1)=PDB__Coords(o) % x
                   yh(nnh+1)=PDB__Coords(o) % y
                   zh(nnh+1)=PDB__Coords(o) % z
                   blist(nnh)=nnh+1
                ELSE
                   nh=nh+1
                   ind_y(nh)=o
                END IF
             END DO
             IF(nna == 3 .AND. nh == 1) WRITE(kprint,*) nnh,blist
             CALL AddH(1,nnh,nh,blist,xh,yh,zh,PDB__Coords(n) %AtmName&
                  &, x0, y0, z0, iret, errmsg)

             IF(iret /= 0) THEN
                CALL Add_Errors(-1,errmsg)
                out=.FALSE.
                RETURN
             END IF

             DO m=1,nh
                q=ind_y(m)
                PDB__Coords(q) % x = x0(m)
                PDB__Coords(q) % y = y0(m)
                PDB__Coords(q) % z = z0(m)
             END DO
             DEALLOCATE(xh,yh,zh,x0,y0,z0,blist,ind_y)
          END IF
       END IF
    END DO
    DEALLOCATE(ind_x, ok_Hyd)
  END FUNCTION AddHydrogens__



!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE AddHydrogens_
