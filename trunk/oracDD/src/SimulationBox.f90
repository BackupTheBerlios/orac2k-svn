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
MODULE SimulationBox

!!$***********************************************************************
!!$   Time-stamp: <2007-01-12 20:39:55 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Jan 12 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program oracDD ----*


!!$======================== DECLARATIONS ================================*


  USE Setup
  USE Inout
  USE Neighbors
  USE AddHydrogens_
  USE SystemTpg
  USE Cell
  USE Solvent
  USE Solute
  USE PDB
  USE IndSequence
  USE SecondarySeq
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f
  USE SystemPrm
  USE AtomCnt
  USE AtomBox
  USE Print_Defs

  IMPLICIT none
  PRIVATE
  PUBLIC :: SimulationBox_,nunits_Slv, Atoms_InBox,SimulationBox__Read
  TYPE(Atombox__), ALLOCATABLE, SAVE, TARGET :: Atoms_InBox(:)
  TYPE(Atombox__), POINTER, SAVE :: Slv_InBox(:)=>NULL(),Slt_InBox(:)=>NULL()

!!$
!!$--- Temporary pointers Erased after call to SimulationBox_
!!$

  TYPE(AtomBox__), POINTER, SAVE :: Total(:)=>NULL()
  TYPE(AtomBox__), POINTER, SAVE :: Slv(:)=>NULL(),Slt(:)=>NULL()&
       &,Sys(:)=>NULL()
  REAL(8), PARAMETER :: Cube_length=6.0D0
  INTEGER, SAVE :: nato_Slv,nato_Slt,nunits_Slv,nunits_SlvOld
  REAL(8), SAVE :: rcut
CONTAINS
  SUBROUTINE SimulationBox_
    INTEGER :: n,m,Begins, Ends,p, nccx,nccy,nccz
    INTEGER :: nato_Slv,nato_Slt
    LOGICAL :: Exist_Slv, Exist_Slt
    TYPE(AtomPdb), POINTER :: PDB__Coords(:)=>NULL()

!!$
!!$--- Need to build a box
!!$

    Exist_Slt=ALLOCATED(Secondary(1) % Line)
    Exist_Slv=ALLOCATED(Secondary(2) % Line)

    IF(Exist_Slt) THEN
       IF(.NOT. ALLOCATED(PDB_Solute)) THEN
          errmsg_f='Solute is defined, but no file to read from '
          CALL Add_Errors(-1,errmsg_f)
          CALL Print_Errors()
       END IF
       IF(.NOT. PDB_('Solute', PDB_Solute, PDB__Coords)) CALL Print_Errors()
       IF(.NOT. AddHydrogens__(PDB__Coords)) CALL Print_Errors()

       
       IF(ASSOCIATED(PDB__Coords)) THEN
          CALL AtomBox_(PDB__Coords,Slt)
          DEALLOCATE(PDB__Coords)
       END IF

       nato_Slt=SIZE(Slt)
       IF(Solvent__Param % Build) CALL AtomBox__Center(Slt)

       IF( .NOT. Exist_Slv) THEN
          m=SIZE(Slt)
          ALLOCATE(Atoms_InBox(m))
          Atoms_InBox=Slt
          CALL AtomBox__ChgFrame(-1,Atoms_InBox)
          Slt_InBox=>atoms_InBox
          DEALLOCATE(Slt)
          RETURN
       END IF
    END IF

    IF(Exist_Slv) THEN
       IF(.NOT. ALLOCATED(PDB_Solvent)) THEN
          errmsg_f='Solvent is defined, but no file to read from '
          CALL Add_Errors(-1,errmsg_f)
          CALL Print_Errors()
       END IF
       IF(.NOT. PDB_('Solvent', PDB_Solvent, PDB__Coords)) CALL Print_Errors()
       IF(.NOT. AddHydrogens__(PDB__Coords)) CALL Print_Errors()

       IF(ASSOCIATED(PDB__Coords)) THEN
          CALL AtomBox_(PDB__Coords,Slv)
          nato_Slv=SIZE(Slv)
          DEALLOCATE(PDB__Coords)
          IF(Solvent__Param % Build) CALL AtomBox__Center(Slv)
       END IF

       IF(Solvent__Param % added /= 0 .OR. (.NOT. Solvent__Param % Build) ) THEN
          m=SIZE(Slv)
          ALLOCATE(Atoms_InBox(m))
          Atoms_InBox=Slv
          CALL AtomBox__ChgFrame(-1,Atoms_InBox)
          Slv_InBox=>atoms_InBox
          nunits_Slv=SIZE(Slv)/nato_Slv
          DEALLOCATE(Slv)
          RETURN
       END IF

       IF(.NOT. AtomBox__BuildSlv(Slv)) CALL Print_Errors()
       
       IF(Exist_Slt) THEN
          n=SIZE(Slt)
          m=SIZE(Slv)
          ALLOCATE(Total(n+m))
          Total(1:n)=Slt
          Total(n+1:)=Slv
          DEALLOCATE(Slt,Slv)
          Slt=>Total(1:n)
          Slv=>Total(n+1:)
       
          nccx=INT(a/Cube_Length)
          nccy=INT(b/Cube_Length)
          nccz=INT(c/Cube_Length)
          rcut=MAXVAL(Slv(:) % sigma)*2.0D0
          IF(.NOT. Neighbors_(rcut, nccx, nccy, nccz)) CALL Print_Errors()

          CALL AtomBox__ChgFrame(-1,Total)
          IF(.NOT. Neighbors__Particles(Total(:) % x&
               &, Total(:) % y, Total(:) % z)) CALL Print_Errors()
          CALL Insert
       ELSE
          m=SIZE(Slv)
          ALLOCATE(Atoms_InBox(m))
          Atoms_InBox=Slv
          CALL AtomBox__ChgFrame(-1,Atoms_InBox)
          Slv_InBox=>atoms_InBox
          nunits_Slv=SIZE(Slv)/nato_Slv
          DEALLOCATE(Slv)
       END IF
    END IF
  CONTAINS
    SUBROUTINE Insert
      INTEGER ::  nx,ny,nz,i,j,k,l,m,n,nv,numcell,iv,jv,kv,nmin&
           &,nn,Size_Total,count_a,o
      REAL(8) :: x1,y1,z1,x2,y2,z2,xx,yy,zz,sqcut,d
      TYPE(AtomBox__), ALLOCATABLE, SAVE :: Temp(:)
      LOGICAL, ALLOCATABLE :: ok_mol(:)
      LOGICAL :: ok
      
      nunits_Slv=SIZE(Slv)/nato_Slv
      ALLOCATE(ok_mol(nunits_Slv))
      ok_mol=.TRUE.

      DO n=1,nunits_Slv
         ok=.TRUE.
         DO m=1,nato_Slv
            nn=(n-1)*nato_Slv + m
            x1=Slv(nn) % x
            y1=Slv(nn) % y
            z1=Slv(nn) % z
            nv=0
            i=Chain_xyz(Slv (nn) % Serial) % i
            j=Chain_xyz(Slv (nn) % Serial) % j
            k=Chain_xyz(Slv (nn) % Serial) % k
            DO o=1,SIZE(Ind_xyz)
               iv=Ind_xyz(o) % i
               jv=Ind_xyz(o) % j
               kv=Ind_xyz(o) % k
               nx=mod(mod(i+iv,nccx)+nccx,nccx)
               ny=mod(mod(j+jv,nccy)+nccy,nccy)
               nz=mod(mod(k+kv,nccz)+nccz,nccz)
               numcell=nz+nccz*(ny+nccy*nx)+1
               IF(numcell > nccx*nccy*nccz) STOP
               l=Head_xyz(numcell)
               DO WHILE(l > 0)
                  IF(l > nato_Slt) THEN
                     l=Chain_xyz(l) % p
                     CYCLE
                  END IF
                  sqcut=(Solute__Exclusion*(Total(l) % sigma+Slv(nn) % sigma)*0.5D0)**2
                  x2=x1-Total(l) % x
                  y2=y1-Total(l) % y
                  z2=z1-Total(l) % z

                  x2=x2-2.0*PBC(x2)
                  y2=y2-2.0*PBC(y2)
                  z2=z2-2.0*PBC(z2)

                  xx=co(1,1)*x2+co(1,2)*y2+co(1,3)*z2
                  yy=           co(2,2)*y2+co(2,3)*z2
                  zz=                      co(3,3)*z2
                  d=xx**2+yy**2+zz**2
                  IF(d < sqcut) THEN
                     ok=.FALSE.
                     EXIT
                  ENDIF
                  l=Chain_xyz(l) % p
               END DO
               IF(.NOT. ok) EXIT
            END DO
            IF(.NOT. ok) EXIT
         END DO
         IF(.NOT. ok) THEN
            ok_mol(n) = .FALSE.
         END IF
      END DO

      nunits_SlvOld=nunits_Slv
      nunits_Slv=COUNT(ok_mol)
      Size_Total=nato_Slt+nunits_Slv*nato_Slv
      WRITE(kprint,'(a)') ' Inserting solute into solvent ====>'
      WRITE(kprint,'(a,i6,a,i6,a,i6, a)') ' Eliminated '&
           &,nunits_SlvOld-nunits_Slv,' solvent units over ', nunits_SlvOld&
           &,' remain ',nunits_Slv,' units '
      ALLOCATE(Temp(SIZE(Total)))
      Temp=Total
      
      ALLOCATE(Atoms_InBox(Size_Total))
      Atoms_InBox(1:nato_Slt)=Temp(1:nato_Slt)

      count_a=0
      DO n=1,nunits_SlvOld
         IF(.NOT. ok_mol(n)) CYCLE
         DO m=1,nato_Slv
            count_a=count_a+1
            nn=(n-1)*nato_Slv + m
            Atoms_InBox(count_a+nato_Slt) = Temp(nn+nato_Slt)
            Atoms_InBox(count_a+nato_Slt) % Serial = count_a+nato_Slt
         END DO
      END DO

      Slt_InBox=>Atoms_InBox(1:nato_Slt)
      Slv_InBox=>Atoms_InBox(nato_Slt+1:)

      DEALLOCATE(Temp,ok_mol,Total)
      Slv=>NULL()
      Slt=>NULL()
    END SUBROUTINE Insert
    FUNCTION PBC(x) RESULT(out)
      REAL(8) :: out
      REAL(8) :: x
      out=ANINT(0.5D0*x)
    END FUNCTION PBC
  END SUBROUTINE SimulationBox_
  SUBROUTINE SimulationBox__Read
    INTEGER :: m
    TYPE(AtomPdb), POINTER :: PDB__Coords(:)=>NULL()

!!$
!!$--- Do not need to build a box: Read from PDB File the entire coordinates set
!!$

    IF(ALLOCATED(Setup__PDB)) THEN
       IF(.NOT. PDB_('System', Setup__PDB, PDB__Coords)) CALL Print_Errors()
       IF(.NOT. AddHydrogens__(PDB__Coords)) CALL Print_Errors()
       CALL AtomBox_(PDB__Coords,Sys)
       DEALLOCATE(PDB__Coords)
       m=SIZE(Sys)
       ALLOCATE(Atoms_InBox(m))
       Atoms_InBox=Sys
       CALL AtomBox__ChgFrame(-1,Atoms_InBox)
       DEALLOCATE(Sys)
       RETURN
    END IF

  END SUBROUTINE SimulationBox__Read
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE SimulationBox
