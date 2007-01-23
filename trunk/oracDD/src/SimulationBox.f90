!!$/---------------------------------------------------------------------\
!!$                                                                      |
!!$  Copyright (C) 2006-2007 Massimo Marchi <Massimo.Marchi@cea.fr>      |
!!$                                                                      |
!!$      This program is free software;  you  can  redistribute  it      |
!!$      and/or modify it under the terms of the GNU General Public      |
!!$      License version 2 as published  by  the  Free  Software         |
!!$      Foundation;                                                     |
!!$                                                                      |
!!$      This program is distributed in the hope that  it  will  be      |
!!$      useful, but WITHOUT ANY WARRANTY; without even the implied      |
!!$      warranty of MERCHANTABILITY or FITNESS  FOR  A  PARTICULAR      |
!!$      PURPOSE.   See  the  GNU  General  Public License for more      |
!!$      details.                                                        |
!!$                                                                      |
!!$      You should have received a copy of the GNU General  Public      |
!!$      License along with this program; if not, write to the Free      |
!!$      Software Foundation, Inc., 59  Temple  Place,  Suite  330,      |
!!$      Boston, MA  02111-1307  USA                                     |
!!$                                                                      |
!!$\---------------------------------------------------------------------/
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
  IMPLICIT none
  PRIVATE
  PUBLIC :: SimulationBox_,nunits_Slv
  TYPE(Atombox__), ALLOCATABLE, SAVE, TARGET :: Atoms_InBox(:)
  TYPE(Atombox__), POINTER, SAVE :: Slv_InBox(:)=>NULL(),Slt_InBox(:)=>NULL()

!!$
!!$--- Temporary pointers Erased after call to SimulationBox_
!!$

  TYPE(AtomBox__), POINTER, SAVE :: Total(:)=>NULL()
  TYPE(AtomBox__), POINTER, SAVE :: Slv(:)=>NULL(),Slt(:)=>NULL()
  REAL(8), PARAMETER :: Cube_length=6.0D0
  INTEGER, SAVE :: nato_Slv,nato_Slt,nunits_Slv,nunits_SlvI
  REAL(8), SAVE :: rcut
CONTAINS
  SUBROUTINE SimulationBox_
    INTEGER :: n,m,Begins, Ends,p, nccx,nccy,nccz
    INTEGER :: nato_Slv,nato_Slt
    LOGICAL :: Exist_Slv, Exist_Slt
    TYPE(AtomPdb), POINTER :: PDB__Coords(:)=>NULL()

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
          CALL PDB__Write(PDB__Coords)
          IF(ALLOCATED(Secondary(2) % Line)) THEN
             CALL AtomBox_(PDB__Coords,Slt)
          END IF
          DEALLOCATE(PDB__Coords)
       END IF
       nato_Slt=SIZE(Slt)
       IF(Solvent__Param % Build) CALL AtomBox__Center(Slt)
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
          CALL PDB__Write(PDB__Coords)
          CALL AtomBox_(PDB__Coords,Slv)
          nato_Slv=SIZE(Slv)
          DEALLOCATE(PDB__Coords)
          IF(Solvent__Param % Build) CALL AtomBox__Center(Slv)
       END IF

       IF(Solvent__Param % added /= 0) RETURN

       IF(.NOT. Solvent__Param % Build) RETURN


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
          IF(.NOT. Neighbors__Atoms(Total(:) % x&
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
            i=cellpi(Slv (nn) % Serial)
            j=cellpj(Slv (nn) % Serial)
            k=cellpk(Slv (nn) % Serial)
            DO o=1,SIZE(Ind_xyz)
               iv=Ind_xyz(o) % i
               jv=Ind_xyz(o) % j
               kv=Ind_xyz(o) % k
               nx=mod(mod(i+iv,nccx)+nccx,nccx)
               ny=mod(mod(j+jv,nccy)+nccy,nccy)
               nz=mod(mod(k+kv,nccz)+nccz,nccz)
               numcell=nz+nccz*(ny+nccy*nx)+1
               IF(numcell > nccx*nccy*nccz) STOP
               l=headp(numcell)
               DO WHILE(l > 0)
                  IF(l > nato_Slt) THEN
                     l=chainp(l)
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
                  l=chainp(l)
               END DO
               IF(.NOT. ok) EXIT
            END DO
            IF(.NOT. ok) EXIT
         END DO
         IF(.NOT. ok) THEN
            ok_mol(n) = .FALSE.
         END IF
      END DO

      nunits_SlvI=nunits_Slv
      nunits_Slv=COUNT(ok_mol)
      Size_Total=nato_Slt+nunits_Slv*nato_Slv
      WRITE(*,'(a)') ' Inserting solute into solvent ====>'
      WRITE(*,'(a,i5,a,i5,a,i5, a)') ' Eliminated '&
           &,nunits_SlvI-nunits_Slv,' solvent units over ', nunits_SlvI&
           &,' remain ',nunits_Slv,' units '
      ALLOCATE(Temp(SIZE(Total)))
      Temp=Total
      
      ALLOCATE(Atoms_InBox(Size_Total))
      Atoms_InBox(1:nato_Slt)=Temp(1:nato_Slt)

      count_a=0
      DO n=1,nunits_Slv
         IF(.NOT. ok_mol(n)) CYCLE
         DO m=1,nato_Slv
            count_a=count_a+1
            nn=(n-1)*nato_Slv + m
            Atoms_InBox(count_a) = Temp(nn+nato_Slt)
            Atoms_InBox(count_a) % Serial = count_a+nato_Slt
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
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE SimulationBox
