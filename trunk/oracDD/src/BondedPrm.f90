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
MODULE ImpropersPrm

!!$***********************************************************************
!!$   Time-stamp: <2007-01-11 18:22:18 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Jan 11 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program oracS ----*

  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f
  USE Resid
  USE SystemTpg
  USE PrmUtilities
  USE Strings
  USE Node
  USE Parameters
  IMPLICIT none
  PRIVATE
  PUBLIC :: ImpropersPrm_,ImpropersPrm__Param, ImpropersPrm__Write&
       &, ImpropersPrm__Read
  TYPE(SystemPrm__Chain), ALLOCATABLE, TARGET, SAVE :: ImpropersPrm__Param(:)
  TYPE(SystemPrm__Chain2), POINTER, SAVE :: share(:)=>NULL()
CONTAINS
  FUNCTION ImpropersPrm_() RESULT(out)
    TYPE(SystemPrm__Chain), POINTER :: out(:)
    INTEGER :: n,i_iTors,nword,nn
    out=>NULL()
    i_iTors=-1
    DO n=1,SIZE(Paras)
       IF(My_Fxm('IMPRO', Paras (n) % Residue )) THEN
          i_iTors=n; EXIT
       END IF
    END DO
    IF(i_iTors == -1) THEN
       errmsg_f='Can''t find improper dihedral parameters in the input file'
       CALL Add_Errors(-1,errmsg_f)
       RETURN
    END IF
    
    IF(.NOT. PrmUtilities__Type(Paras (i_iTors) % line, Tpg % imph)) RETURN
    IF(.NOT. PrmUtilities__Gather('Improper Torsions ',Tpg % imph, share)) RETURN
    IF(ASSOCIATED(share)) THEN
       ALLOCATE(ImpropersPrm__Param(SIZE(share)))
       DO n=1,SIZE(Share)
          IF(ASSOCIATED(share(n) % g)) THEN
             ALLOCATE(ImpropersPrm__Param(n) % g (SIZE(share(n) % g)))
             ImpropersPrm__Param(n) % g = share(n) % g
          END IF
          ImpropersPrm__Param(n) % pt = share(n) % pt
          DEALLOCATE(share(n) % g)
       END DO
       DEALLOCATE(share)
    END IF
    out=>ImpropersPrm__Param
    IF(ALLOCATED(ImpropersPrm__Param)) THEN
       WRITE(*,*) 'Total improper dihedral Parameters No. =====>'&
            &,SIZE(ImpropersPrm__Param)
    END IF
  END FUNCTION ImpropersPrm_
  SUBROUTINE ImpropersPrm__Write
    INTEGER :: n,p,m    
    p=0
    IF(ALLOCATED(ImpropersPrm__Param)) p=SIZE(ImpropersPrm__Param)
    WRITE(kbinary) p
    IF(p == 0) RETURN
    DO n=1,p
       m=0
       IF(ALLOCATED(ImpropersPrm__Param(n) % g)) m=SIZE(ImpropersPrm__Param(n) % g)
       WRITE(kbinary) m
       IF(m == 0) CYCLE
       WRITE(kbinary) ImpropersPrm__Param(n) % pt, ImpropersPrm__Param(n) % g
    END DO
  END SUBROUTINE ImpropersPrm__Write
  FUNCTION ImpropersPrm__Read() RESULT(out)
    TYPE(SystemPrm__Chain), POINTER :: out(:)
    INTEGER :: n,p,m
    out=>NULL()
    READ(kbinary,ERR=100,END=200) p
    IF(p == 0) RETURN
    ALLOCATE(ImpropersPrm__Param(p))
    DO n=1,p
       READ(kbinary,ERR=100,END=200) m
       IF(m == 0) CYCLE
       ALLOCATE(ImpropersPrm__Param(n) % g (m))
       READ(kbinary,ERR=100,END=200) ImpropersPrm__Param(n) % pt&
            &, ImpropersPrm__Param(n) % g
    END DO
    out=>ImpropersPrm__Param
    RETURN
100 errmsg_f='Error while reading Lennard-Jones Parameters'
    CALL Add_Errors(-1,errmsg_f)
    RETURN
200 errmsg_f='End of file found while reading Lennard-Jones Parameters'
    CALL Add_Errors(-1,errmsg_f)
    RETURN    
  END FUNCTION ImpropersPrm__Read
  
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE ImpropersPrm
MODULE TorsionsPrm

!!$***********************************************************************
!!$   Time-stamp: <2007-01-11 17:30:57 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Jan 11 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program oracS ----*

  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f
  USE Resid
  USE SystemTpg
  USE PrmUtilities
  USE Strings
  USE Node
  USE Parameters
  IMPLICIT none
  PRIVATE
  PUBLIC :: TorsionsPrm_,TorsionsPrm__Param, TorsionsPrm__Read&
       &, Torsionsprm__Write
  TYPE(SystemPrm__Chain), ALLOCATABLE, TARGET, SAVE :: TorsionsPrm__Param(:)
  TYPE(SystemPrm__Chain2), POINTER, SAVE :: share(:)=>NULL()
CONTAINS
  FUNCTION TorsionsPrm_() RESULT(out)
    TYPE(SystemPrm__Chain), POINTER :: out(:)
    INTEGER :: n,i_Tors,nword,nn
    out=>NULL()
    i_Tors=-1
    DO n=1,SIZE(Paras)
       IF(My_Fxm('DIHED', Paras (n) % Residue )) THEN
          i_Tors=n; EXIT
       END IF
    END DO
    
    DO n=1,SIZE(Paras)
       IF(My_Fxm('TORSIO', Paras (n) % Residue )) THEN
          i_Tors=n; EXIT
       END IF
    END DO
    
    IF(i_Tors == -1) THEN
       errmsg_f='Can''t find dihedral parameters in the input file'
       CALL Add_Errors(-1,errmsg_f)
       RETURN
    END IF
    
    IF(.NOT. PrmUtilities__Type(Paras (i_Tors) % line, Tpg % dihed)) RETURN
    IF(.NOT. PrmUtilities__Gather('Torsions ',Tpg % dihed, share)) RETURN
    IF(ASSOCIATED(share)) THEN
       ALLOCATE(TorsionsPrm__Param(SIZE(share)))
       DO n=1,SIZE(Share)
          IF(ASSOCIATED(share(n) % g)) THEN
             ALLOCATE(TorsionsPrm__Param(n) % g (SIZE(share(n) % g)))
             TorsionsPrm__Param(n) % g = share(n) % g
          END IF
          TorsionsPrm__Param(n) % pt = share(n) % pt
          DEALLOCATE(share(n) % g)
       END DO
       DEALLOCATE(share)
    END IF
    out=>TorsionsPrm__Param
    IF(ALLOCATED(TorsionsPrm__Param)) THEN
       WRITE(*,*) 'Total Torsion Parameters No. =====>',SIZE(TorsionsPrm__Param)
    END IF
  END FUNCTION TorsionsPrm_
  SUBROUTINE TorsionsPrm__Write
    INTEGER :: n,p,m    
    p=0
    IF(ALLOCATED(TorsionsPrm__Param)) p=SIZE(TorsionsPrm__Param)
    WRITE(kbinary) p
    IF(p == 0) RETURN
    DO n=1,p
       m=0
       IF(ALLOCATED(TorsionsPrm__Param(n) % g)) m=SIZE(TorsionsPrm__Param(n) % g)
       WRITE(kbinary) m
       IF(m == 0) CYCLE
       WRITE(kbinary) TorsionsPrm__Param(n) % pt, TorsionsPrm__Param(n) % g
    END DO
  END SUBROUTINE TorsionsPrm__Write
  FUNCTION TorsionsPrm__Read() RESULT(out)
    TYPE(SystemPrm__Chain), POINTER :: out(:)
    INTEGER :: n,p,m
    out=>NULL()
    READ(kbinary,ERR=100,END=200) p
    IF(p == 0) RETURN
    ALLOCATE(TorsionsPrm__Param(p))
    DO n=1,p
       READ(kbinary,ERR=100,END=200) m
       IF(m == 0) CYCLE
       ALLOCATE(TorsionsPrm__Param(n) % g (m))
       READ(kbinary,ERR=100,END=200) TorsionsPrm__Param(n) % pt&
            &, TorsionsPrm__Param(n) % g
    END DO
    out=>TorsionsPrm__Param
    RETURN
100 errmsg_f='Error while reading Lennard-Jones Parameters'
    CALL Add_Errors(-1,errmsg_f)
    RETURN
200 errmsg_f='End of file found while reading Lennard-Jones Parameters'
    CALL Add_Errors(-1,errmsg_f)
    RETURN    
  END FUNCTION TorsionsPrm__Read
  
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE TorsionsPrm
MODULE AnglesPrm

!!$***********************************************************************
!!$   Time-stamp: <2007-01-11 17:24:19 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Jan 11 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program oracS ----*

  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f
  USE Resid
  USE SystemTpg
  USE PrmUtilities
  USE Strings
  USE Node
  USE Parameters
  IMPLICIT none
  PRIVATE
  PUBLIC :: AnglesPrm_,AnglesPrm__Param,AnglesPrm__Read, AnglesPrm__Write
  TYPE(SystemPrm__Chain), ALLOCATABLE, TARGET, SAVE :: AnglesPrm__Param(:)
  TYPE(SystemPrm__Chain2), POINTER, SAVE :: share(:)=>NULL()
CONTAINS
  FUNCTION AnglesPrm_() RESULT(out)
    TYPE(SystemPrm__Chain), POINTER :: out(:)
    INTEGER :: n,i_bends,nword,nn
    out=>NULL()
    i_Bends=-1
    DO n=1,SIZE(Paras)
       IF(My_Fxm('ANGLES', Paras (n) % Residue )) THEN
          i_bends=n; EXIT
       END IF
    END DO
    IF(i_Bends == -1) THEN
       DO n=1,SIZE(Paras)
          IF(My_Fxm('BENDIN', Paras (n) % Residue )) THEN
             i_bends=n; EXIT
          END IF
       END DO
    END IF
    IF(i_Bends == -1) THEN
       errmsg_f='Can''t find bending parameters in the input file'
       CALL Add_Errors(-1,errmsg_f)
       RETURN
    END IF

    IF(.NOT. PrmUtilities__Type(Paras (i_Bends) % line, Tpg % angles)) RETURN
    IF(.NOT. PrmUtilities__Gather('Angles ',Tpg % angles, share)) RETURN
    IF(ASSOCIATED(share)) THEN
       ALLOCATE(AnglesPrm__Param(SIZE(share)))
       DO n=1,SIZE(Share)
          IF(ASSOCIATED(share(n) % g)) THEN
             ALLOCATE(AnglesPrm__Param(n) % g (SIZE(share(n) % g)))
             AnglesPrm__Param(n) % g = share(n) % g
          END IF
          AnglesPrm__Param(n) % pt = share(n) % pt
          DEALLOCATE(share(n) % g)
       END DO
       DEALLOCATE(share)
    END IF
    out=>AnglesPrm__Param
    IF(ALLOCATED(AnglesPrm__Param)) THEN
       WRITE(*,*) 'Total Angle Parameters No. =====>',SIZE(AnglesPrm__Param)
    END IF
  END FUNCTION AnglesPrm_
  SUBROUTINE AnglesPrm__Write
    INTEGER :: n,p,m    
    p=0
    IF(ALLOCATED(AnglesPrm__Param)) p=SIZE(AnglesPrm__Param)
    WRITE(kbinary) p
    IF(p == 0) RETURN
    DO n=1,p
       m=0
       IF(ALLOCATED(AnglesPrm__Param(n) % g)) m=SIZE(AnglesPrm__Param(n) % g)
       WRITE(kbinary) m
       IF(m == 0) CYCLE
       WRITE(kbinary) AnglesPrm__Param(n) % pt, AnglesPrm__Param(n) % g
    END DO
  END SUBROUTINE AnglesPrm__Write
  FUNCTION AnglesPrm__Read() RESULT(out)
    TYPE(SystemPrm__Chain), POINTER :: out(:)
    INTEGER :: n,p,m
    out=>NULL()
    READ(kbinary,ERR=100,END=200) p
    IF(p == 0) RETURN
    ALLOCATE(AnglesPrm__Param(p))
    DO n=1,p
       READ(kbinary,ERR=100,END=200) m
       IF(m == 0) CYCLE
       ALLOCATE(AnglesPrm__Param(n) % g (m))
       READ(kbinary,ERR=100,END=200) AnglesPrm__Param(n) % pt&
            &, AnglesPrm__Param(n) % g
    END DO
    out=>AnglesPrm__Param
    RETURN
100 errmsg_f='Error while reading Lennard-Jones Parameters'
    CALL Add_Errors(-1,errmsg_f)
    RETURN
200 errmsg_f='End of file found while reading Lennard-Jones Parameters'
    CALL Add_Errors(-1,errmsg_f)
    RETURN    
  END FUNCTION AnglesPrm__Read
  
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE AnglesPrm
MODULE BondsPrm

!!$***********************************************************************
!!$   Time-stamp: <2007-01-11 17:14:09 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Jan 11 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program oracS ----*

  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f
  USE Resid
  USE SystemTpg
  USE PrmUtilities
  USE Strings
  USE Node
  USE Parameters
  IMPLICIT none
  PRIVATE
  PUBLIC :: BondsPrm_,BondsPrm__Param, BondsPrm__Read, BondsPrm__Write
  TYPE(SystemPrm__Chain), ALLOCATABLE, TARGET, SAVE :: BondsPrm__Param(:)
  TYPE(SystemPrm__Chain2), POINTER, SAVE :: share(:)=>NULL()
CONTAINS
  FUNCTION BondsPrm_() RESULT(out)
    TYPE(SystemPrm__Chain), POINTER :: out(:)
    INTEGER :: n,i_bonds,nword,nn
    out=>NULL()
    i_Bonds=-1
    DO n=1,SIZE(Paras)
       IF(My_Fxm('BOND', Paras (n) % Residue )) THEN
          i_bonds=n; EXIT
       END IF
    END DO
    IF(i_Bonds == -1) THEN
       errmsg_f='Can''t find stretching parameters in the input file'
       CALL Add_Errors(-1,errmsg_f)
       RETURN
    END IF
    IF(.NOT. PrmUtilities__Type(Paras (i_Bonds) % line, Tpg % bonds)) RETURN
    IF(.NOT. PrmUtilities__Gather('Bonds',Tpg % bonds, share)) RETURN
    IF(ASSOCIATED(share)) THEN
       ALLOCATE(BondsPrm__Param(SIZE(share)))
       DO n=1,SIZE(Share)
          IF(ASSOCIATED(share(n) % g)) THEN
             ALLOCATE(BondsPrm__Param(n) % g (SIZE(share(n) % g)))
             BondsPrm__Param(n) % g = share(n) % g
          END IF
          BondsPrm__Param(n) % pt = share(n) % pt
          DEALLOCATE(share(n) % g)
       END DO
       DEALLOCATE(share)
    END IF
    out=>BondsPrm__Param
    IF(ALLOCATED(BondsPrm__Param)) THEN
       WRITE(*,*) 'Total Stretching Parameters No. =====>',SIZE(BondsPrm__Param)
    END IF
  END FUNCTION BondsPrm_
  SUBROUTINE BondsPrm__Write
    INTEGER :: n,p,m    
    p=0
    IF(ALLOCATED(BondsPrm__Param)) p=SIZE(BondsPrm__Param)
    WRITE(kbinary) p
    IF(p == 0) RETURN
    DO n=1,p
       m=0
       IF(ALLOCATED(BondsPrm__Param(n) % g)) m=SIZE(BondsPrm__Param(n) % g)
       WRITE(kbinary) m
       IF(m == 0) CYCLE
       WRITE(kbinary) BondsPrm__Param(n) % pt, BondsPrm__Param(n) % g
    END DO
  END SUBROUTINE BondsPrm__Write
  FUNCTION BondsPrm__Read() RESULT(out)
    TYPE(SystemPrm__Chain), POINTER :: out(:)
    INTEGER :: n,p,m
    out=>NULL()
    READ(kbinary,ERR=100,END=200) p
    IF(p == 0) RETURN
    ALLOCATE(BondsPrm__Param(p))
    DO n=1,p
       READ(kbinary,ERR=100,END=200) m
       IF(m == 0) CYCLE
       ALLOCATE(BondsPrm__Param(n) % g (m))
       READ(kbinary,ERR=100,END=200) BondsPrm__Param(n) % pt&
            &, BondsPrm__Param(n) % g
    END DO
    out=>BondsPrm__Param
    RETURN
100 errmsg_f='Error while reading Lennard-Jones Parameters'
    CALL Add_Errors(-1,errmsg_f)
    RETURN
200 errmsg_f='End of file found while reading Lennard-Jones Parameters'
    CALL Add_Errors(-1,errmsg_f)
    RETURN    
  END FUNCTION BondsPrm__Read

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE BondsPrm
