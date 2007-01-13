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
  IMPLICIT none
  PRIVATE
  PUBLIC :: ImpropersPrm_,ImpropersPrm__Param
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
    WRITE(*,*) 'Total improper dihedral Parameters No. =====>',SIZE(ImpropersPrm__Param)
  END FUNCTION ImpropersPrm_
  
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
  IMPLICIT none
  PRIVATE
  PUBLIC :: TorsionsPrm_,TorsionsPrm__Param
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
    WRITE(*,*) 'Total Torsion Parameters No. =====>',SIZE(TorsionsPrm__Param)
  END FUNCTION TorsionsPrm_
  
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
  IMPLICIT none
  PRIVATE
  PUBLIC :: AnglesPrm_,AnglesPrm__Param
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
    WRITE(*,*) 'Total Angle Parameters No. =====>',SIZE(AnglesPrm__Param)
  END FUNCTION AnglesPrm_
  
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
  IMPLICIT none
  PRIVATE
  PUBLIC :: BondsPrm_,BondsPrm__Param
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
    WRITE(*,*) 'Total Stretching Parameters No. =====>',SIZE(BondsPrm__Param)
  END FUNCTION BondsPrm_
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE BondsPrm
