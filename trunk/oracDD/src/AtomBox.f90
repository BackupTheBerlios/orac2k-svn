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
MODULE AtomBox

!!$***********************************************************************
!!$   Time-stamp: <2007-01-23 17:12:57 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Wed Jan 17 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program oracDD ----*

  USE Units
  USE Cell
  USE PDB
  USE AtomCnt
  USE SystemPrm
  USE Solvent
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, error_args, errmsg_f
  USE Print_Defs
  IMPLICIT none
  PRIVATE
  PUBLIC AtomBox_, AtomBox__BuildSlv, AtomBox__, AtomBox__ChgFrame&
       &, AtomBox__Center, AtomBox__AtomPDB
  TYPE :: AtomBox__
     INTEGER :: Serial
     REAL(8) :: sigma
     REAL(8) :: x,y,z
  END TYPE AtomBox__
  INTEGER, SAVE :: nSlv=0
CONTAINS
  SUBROUTINE AtomBox_(PDB__Coords, Coords)
    TYPE(AtomPdb), POINTER :: PDB__Coords(:)
    TYPE(AtomBox__), POINTER :: Coords(:)
    INTEGER :: n,m
    IF(ASSOCIATED(PDB__Coords)) ALLOCATE(Coords(SIZE(PDB__Coords)))
    Coords(:) % x = PDB__Coords(:) % x
    Coords(:) % y = PDB__Coords(:) % y
    Coords(:) % z = PDB__Coords(:) % z
    Coords(:) % Serial = PDB__Coords(:) % Serial
    DO n=1,SIZE(PDB__Coords)
       m=AtomCnts(PDB__Coords(n) % Serial) % Id_Type
       Coords(n) % Sigma = Prm % LJ % Par_SE (m) % g(1)
    END DO
  END SUBROUTINE AtomBox_
  FUNCTION AtomBox__BuildSlv(Coords) RESULT(out)
    TYPE(AtomBox__), POINTER :: Coords(:)
    INTEGER, INTRINSIC :: count
    LOGICAL :: out
    INTEGER :: Ic, Jc, Kc
    INTEGER :: n,i_Type,o,p,q,m,l
    REAL(8) :: cube,cube_side,a0,b0,c0,volume0,x,y,z,x_0,y_0,z_0&
         &,cube_x,cube_y,cube_z,ko(3,3),vol0,vol
    INTEGER :: i_dvol(1),nCoords,offset,n_Tot
    REAL(8), ALLOCATABLE :: Dvol(:)
    TYPE(AtomBox__), ALLOCATABLE :: Coords0(:)
    out=.TRUE.


    IF(Solvent__Param % rho > 0.0D0) THEN
       ALLOCATE(Dvol(SIZE(Solvent__Box)))
       vol=1.0D0/Solvent__Param % rho
       DO n=1,SIZE(Solvent__Box)
          cube=(vol)*DBLE(Solvent__Box(n) % nt)
          N_tot=NINT(volume/cube)
          Ic=NINT(((a**2/(b*c))*DBLE(N_Tot))**(1.0D0/3.0D0))
          Jc=NINT(Ic*b/a)
          Kc=NINT(Ic*c/a)
          vol0=volume/DBLE((Ic*Jc*Kc)*Solvent__Box(n) % nt)
          Dvol(n)=ABS(vol-vol0)
       END DO
       i_dvol=MINLOC(Dvol)
       i_Type=i_Dvol(1)
       cube=(vol)*DBLE(Solvent__Box(i_Type) % nt)
       N_tot=NINT(volume/cube)
       Ic=NINT(((a**2/(b*c))*DBLE(N_Tot))**(1.0D0/3.0D0))
       Jc=NINT(Ic*b/a)
       Kc=NINT(Ic*c/a)
       DEALLOCATE(Dvol)
    ELSE IF(COUNT(Solvent__Param % replicate /= 0) /= 0) THEN
       DO n=1,SIZE(Solvent__Box)
          IF(TRIM(Solvent__Box(n) % type) == TRIM(Solvent__Param % Cell_Type)) THEN
             i_Type=n
             EXIT
          END IF
       END DO
       Ic=Solvent__Param % Replicate(1)
       Jc=Solvent__Param % Replicate(2)
       Kc=Solvent__Param % Replicate(3)
    ELSE
       errmsg_f='Can''t build Simulation box without solvent density&
            & or replicate parameters '
       CALL Add_Errors(-1,errmsg_f)
       out=.FALSE.; RETURN
    END IF
    
    cube_x=boxl/DBLE(Ic)
    cube_y=boxl/DBLE(Jc)
    cube_z=boxl/DBLE(Kc)
    
    nCoords=SIZE(Coords)
    

    ALLOCATE(Coords0 (Solvent__Box(i_Type) % nt*Ic*Jc*Kc*SIZE(Coords)))
    WRITE(kprint,'(a,f12.3)') ' Solvent imposed molecular volume ======> '&
         &,volume/DBLE(Solvent__Box(i_Type) % nt*Ic*Jc*Kc)

    m=0
    offset=0
    DO o=1,Ic
       x_0=DBLE(o-1)
       DO p=1,Jc
          y_0=DBLE(p-1)
          DO q=1,Kc
             z_0=DBLE(q-1)

             DO n=1,Solvent__Box(i_Type) % nt
                x=(x_0+Solvent__Box(i_Type) % T(1,n))*cube_x-1.0D0
                y=(y_0+Solvent__Box(i_Type) % T(2,n))*cube_y-1.0D0
                z=(z_0+Solvent__Box(i_Type) % T(3,n))*cube_z-1.0D0
                DO l=1,nCoords
                   m=m+1
                   Coords0(m) % x = Coords(l) % x + co(1,1)*x + co(1,2)*y + co(1,3)*z 
                   Coords0(m) % y = Coords(l) % y + co(2,1)*x + co(2,2)*y + co(2,3)*z 
                   Coords0(m) % z = Coords(l) % z + co(3,1)*x + co(3,2)*y + co(3,3)*z
                   Coords0(m) % sigma = Coords(l) % sigma
                   Coords0(m) % Serial = Coords(l) % Serial + offset
                END DO
                offset=offset+nCoords
             END DO
          END DO
       END DO
    END DO

    DEALLOCATE(Coords)
    ALLOCATE(Coords(SIZE(Coords0)))
    Coords=Coords0
    DEALLOCATE(Coords0)
  END FUNCTION AtomBox__BuildSlv  
  SUBROUTINE AtomBox__ChgFrame(Dir,Coords)
    INTEGER :: Dir
    TYPE(AtomBox__) :: Coords(:)
    TYPE(AtomBox__), ALLOCATABLE :: Coords0(:)

    ALLOCATE(Coords0(SIZE(Coords)))
    Coords0=Coords
    IF(Dir > 0) THEN
       Coords(:) % x = co(1,1)*Coords0(:) % x + co(1,2)*Coords0(:) % y &
            &+ co(1,3)*Coords0(:) % z
       Coords(:) % y = co(2,1)*Coords0(:) % x + co(2,2)*Coords0(:) % y &
            &+ co(2,3)*Coords0(:) % z
       Coords(:) % z = co(3,1)*Coords0(:) % x + co(3,2)*Coords0(:) % y &
            &+ co(3,3)*Coords0(:) % z
    ELSE
       Coords(:) % x = oc(1,1)*Coords0(:) % x + oc(1,2)*Coords0(:) % y &
            &+ oc(1,3)*Coords0(:) % z
       Coords(:) % y = oc(2,1)*Coords0(:) % x + oc(2,2)*Coords0(:) % y &
            &+ oc(2,3)*Coords0(:) % z
       Coords(:) % z = oc(3,1)*Coords0(:) % x + oc(3,2)*Coords0(:) % y &
            &+ oc(3,3)*Coords0(:) % z
    END IF
    
    DEALLOCATE(Coords0)
  END SUBROUTINE AtomBox__ChgFrame
  SUBROUTINE AtomBox__Center(Coords)
    TYPE(AtomBox__), POINTER :: Coords(:)
    INTEGER :: nCoords
    REAL(8) :: xcm,ycm,zcm

    nCoords=SIZE(Coords)
    xcm=SUM(Coords(:) % x)/DBLE(nCoords)
    ycm=SUM(Coords(:) % y)/DBLE(nCoords)
    zcm=SUM(Coords(:) % z)/DBLE(nCoords)
    Coords(:) % x = Coords(:) % x - xcm
    Coords(:) % y = Coords(:) % y - ycm
    Coords(:) % z = Coords(:) % z - zcm
  END SUBROUTINE AtomBox__Center
  
  SUBROUTINE AtomBox__AtomPDB(Coords, PDB__Coords)
    TYPE(AtomBox__) :: Coords(:)
    TYPE(AtomPdb), POINTER :: PDB__Coords(:)
    INTEGER :: n,m
    ALLOCATE(PDB__Coords(SIZE(Coords)))
    PDB__Coords(:) % x = Coords(:) % x
    PDB__Coords(:) % y = Coords(:) % y
    PDB__Coords(:) % z = Coords(:) % z
    PDB__Coords(:) % Serial = Coords(:) % Serial
  END SUBROUTINE AtomBox__AtomPDB
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE AtomBox
