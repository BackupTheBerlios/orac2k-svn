MODULE PDBs_Mod

!!$***********************************************************************
!!$   Time-stamp: <2005-10-19 10:15:33 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Sat Oct 15 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

  USE INPUT_Mod, ONLY: Read_String, Parser, err_open,err_end,err_unr&
       &,err_fnf,err_args
  TYPE PDB_FIT
     LOGICAL :: fit=.FALSE.
     CHARACTER(3) :: type
  END TYPE PDB_FIT
  TYPE Atom_pdb
     CHARACTER(7) :: label
     CHARACTER(8) :: Res
     INTEGER :: ResNo
     REAL(8) :: charge
     REAL(8) :: xp,yp,zp
  END TYPE Atom_pdb

  LOGICAL, SAVE :: PDBs=.FALSE.
  INTEGER, DIMENSION(2), SAVE :: Residue_Range=(/-1,-1/)
  INTEGER, SAVE :: kpdb=0,n_write=1,n_compute=1
  INTEGER, SAVE :: natom, counter
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: index
  REAL(8), SAVE :: sunitc=2.68301916419150834D0
  CHARACTER(80) :: filename
  CHARACTER(3), DIMENSION(23), SAVE :: protein=(/ 'ala','arg','asn'&
       &,'asp','cys','gln','glu','gly','hsd','hse','hsp','ile'&
       &,'leu','lys','met','phe','pro','ser','thr','trp','tyr'&
       &,'val','cal'/)

  TYPE(PDB_FIT), SAVE  :: Template
  TYPE(Atom_pdb), DIMENSION(:), ALLOCATABLE, SAVE :: Coord_avg
CONTAINS
  SUBROUTINE Initialize(charge,prsymb,beta,resa,resb,ntap)

!!$======================== DECLARATIONS ================================*

    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*
    
    INTEGER :: ntap,resa(*),resb(*)
    CHARACTER(8) :: prsymb(*)
    CHARACTER(7) :: beta(*)
    REAL(8) :: charge(*)

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: i,j
    CHARACTER(3) :: Residue

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    natom=0
    IF(Residue_Range(1) /= -1 .AND. Residue_Range(2) /= -1) THEN
       DO i=1,ntap
          IF(resa(i) >= Residue_Range(1) .AND. resa(i) <=&
               & Residue_Range(2)) THEN
             natom=natom+1
          END IF
       END DO
       ALLOCATE(index(natom))
       natom=0
       DO i=1,ntap
          IF(resa(i) >= Residue_Range(1) .AND. resa(i) <=&
               & Residue_Range(2)) THEN
             natom=natom+1
             index(natom)=i
          END IF
       END DO
    ELSE
       DO i=1,ntap
          Residue=prsymb(resb(i))(1:3)
          DO j=1,LEN(protein)
             IF(Residue == protein(j) .OR. Residue ==&
                  & upper(protein(j))) THEN
                natom=natom+1
             END IF
          END DO
       END DO
       ALLOCATE(index(natom))
       natom=0
       DO i=1,ntap
          Residue=prsymb(resb(i))(1:3)
          DO j=1,LEN(protein)
             IF(Residue == protein(j) .OR. Residue ==&
                  & upper(protein(j))) THEN
                natom=natom+1
                index(natom)=i
             END IF
          END DO
       END DO
    END IF
    ALLOCATE(Coord_avg(natom))
    DO i=1,natom
       Coord_avg(i)%Charge=charge(index(i))*sunitc
       Coord_avg(i)%label=beta(index(i))
       Coord_avg(i)%Res=prsymb(resb(index(i)))
       Coord_avg(i)%Resno=resa(index(i))
       Coord_avg(i)%xp=0.0D0
       Coord_avg(i)%yp=0.0D0
       Coord_avg(i)%xp=0.0D0
    END DO

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*
  CONTAINS
    FUNCTION upper(string) RESULT(out)
      IMPLICIT NONE 
      CHARACTER(3) :: string,out
      INTEGER :: upper_to_lower,len_string,i
      
      upper_to_lower = IACHAR("a") - IACHAR("A")
      out=string
      len_string=LEN_TRIM(out)
      DO i=1,len_string
         SELECT CASE (out(i:i))
         CASE ('a':'z')
            out(i:i) = ACHAR(IACHAR(out(i:i)) - upper_to_lower)
         END SELECT
      END DO
    END FUNCTION upper
  END SUBROUTINE Initialize
  SUBROUTINE Compute(xp0,yp0,zp0)

!!$======================== DECLARATIONS ================================*

    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*
    
    INTEGER :: ntap
    REAL(8) :: xp0(*),yp0(*),zp0(*)

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: i

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    IF(.NOT. Template%fit) THEN
       DO i=1,natom
          Coord_Avg(i)%xp=Coord_Avg(i)%xp+xp0(index(i))
          Coord_Avg(i)%yp=Coord_Avg(i)%yp+yp0(index(i))
          Coord_Avg(i)%zp=Coord_Avg(i)%zp+zp0(index(i))
       END DO
    END IF
    counter=counter+1

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

  END SUBROUTINE Compute
  SUBROUTINE Write_it(fstep)

!!$======================== DECLARATIONS ================================*

    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*
    
    REAL(8) :: fstep

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: i,ResNo,k
    REAL(8) :: charge,xb,yb,zb
    CHARACTER(5) :: bet
    CHARACTER(3) :: Res
    

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    REWIND(kpdb)
    WRITE(kpdb,2) fstep
    k=0
    DO i=1,natom
       xb=Coord_avg(i)%xp/DBLE(counter)
       yb=Coord_avg(i)%yp/DBLE(counter)
       zb=Coord_avg(i)%zp/DBLE(counter)
       bet(1:5)=Coord_avg(i)%Label(1:5)
       CALL low_up(bet,5)
       resno=Coord_avg(i)%Resno
       res(1:3)=Coord_avg(i)%Res(1:3)
       charge=Coord_avg(i)%Charge
       WRITE(kpdb,1)'ATOM  ',i,bet(1:5),res,resno,xb,yb,zb&
            &,charge,DBLE(k)
    END DO

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

    WRITE(kpdb,'(a)')'TER  '
1   FORMAT(a5,i6,1x,a5,a3,1x,i5,4x,3f8.3,2f6.2)
2   FORMAT('COMMENT 1 Averaged Configuration at time step ',f11.2&
         &,'               ')
  END SUBROUTINE Write_it
  
  INCLUDE 'PDBs_Read.f90'
  
END MODULE PDBs_Mod
