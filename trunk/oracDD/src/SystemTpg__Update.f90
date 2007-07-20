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

  SUBROUTINE SystemTpg__Update(New_Units)

!!$***********************************************************************
!!$   Time-stamp: <2007-01-22 14:52:06 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Mon Jan 22 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program  oracDD ----*
    INTEGER, INTENT(IN) :: New_Units
    INTEGER, POINTER :: SltSlv_Res(:,:)=>NULL()
    INTEGER, POINTER :: SltSlv_Grp(:,:)=>NULL()
    INTEGER :: n,m,l,p,q,offset,nato,Res_Begins,Res_Ends&
         &,Grp_Begins,Grp_Ends
    INTEGER, ALLOCATABLE :: Temp(:,:)
    TYPE(Chain), DIMENSION(:), ALLOCATABLE :: Temp1
    INTEGER :: o(2),o1,s,natop

    SltSlv_Res=>IndSequence__SltSlv_Res()
    Res_Begins=SltSlv_Res(1,2)
    Res_Ends=SltSlv_Res(2,2)
    natop=Res_Atm(2,Res_Ends)-Res_Atm(1,Res_Begins)+1

    SltSlv_Grp=>IndSequence__SltSlv_Grp()
    Grp_Begins=SltSlv_Grp(1,2)
    Grp_Ends=SltSlv_Grp(2,2)

    WRITE(kprint,*) '<---- Update topology Adding solvent molecules if any ---->' 
    WRITE(kprint,*) 'Atoms No. =====>',SIZE(atm_cnt)


    o=SHAPE(Tpg % Grp_Atm)
    ALLOCATE(Temp(o(1), o(2)))
    
    Temp=Tpg % Grp_Atm
    DEALLOCATE(Tpg % Grp_Atm)
    p=Grp_Ends-Grp_Begins+1
    ALLOCATE(Tpg % Grp_Atm (o(1), o(2)+(New_Units-1)*p))
    offset=0
    nato=0
    DO n=1,New_Units
       DO m=Grp_Begins,Grp_Ends
          Tpg % Grp_Atm(:,m+offset)=Temp(:,m)+nato
       END DO
       nato=nato+natop
       offset=offset+p
    END DO
    DEALLOCATE(Temp)

    o=SHAPE(Tpg % Res_Atm)
    ALLOCATE(Temp(o(1), o(2)))
    
    Temp=Tpg % Res_Atm
    DEALLOCATE(Tpg % Res_Atm)

    p=Res_Ends-Res_Begins+1
    ALLOCATE(Tpg % Res_Atm (o(1), o(2)+(New_Units-1)*p))
    offset=0
    nato=0
    DO n=1,New_Units
       DO m=Res_Begins,Res_Ends
          Tpg % Res_Atm(:,m+offset)=Temp(:,m)+nato
       END DO
       nato=nato+natop
       offset=offset+p
    END DO
    DEALLOCATE(Temp)

    WRITE(kprint,*) 'Residue No. =====>',SIZE(Tpg % Res_Atm,2)
    WRITE(kprint,*) 'Group No. =====>',SIZE(Tpg % Grp_Atm,2)

    DEALLOCATE(Tpg % atm)
    ALLOCATE(Tpg % atm(SIZE(atm_cnt)))
    DO n=1,SIZE(atm_cnt)
       Tpg % atm (n) % a => atm_cnt(n)
    END DO

    IF(ALLOCATED(Tpg % bonds)) THEN
       o=SHAPE(Tpg % bonds)
       ALLOCATE(Temp(o(1), o(2)))
       p=Tpg % S_bonds
       q=o(2)-Tpg % S_bonds+1
       Temp = Tpg % bonds
       DEALLOCATE(Tpg % bonds)
       ALLOCATE(Tpg % bonds(o(1), o(2)+(New_Units-1)*q))
       Tpg % bonds(:,1:o(2))=Temp
       offset=0
       nato=0
       DO n=1,New_Units
          DO m=p,p+q-1
             Tpg % bonds(:,m+offset)=Temp(:,m)+nato
          END DO
          offset=offset+q
          nato=nato+natop
       END DO
       WRITE(kprint,*) 'New Bonds No. =====>',SIZE(Tpg % bonds,2)
       DEALLOCATE(Temp)
    END IF
    IF(ALLOCATED(Tpg % Angles)) THEN
       o=SHAPE(Tpg % Angles)
       ALLOCATE(Temp(o(1), o(2)))
       p=Tpg % S_Angles
       q=o(2)-Tpg % S_Angles+1
       Temp = Tpg % Angles
       DEALLOCATE(Tpg % Angles)
       ALLOCATE(Tpg % Angles(o(1), o(2)+(New_Units-1)*q))
       Tpg % Angles(:,1:o(2))=Temp
       offset=0
       nato=0
       DO n=1,New_Units
          DO m=p,p+q-1
             Tpg % Angles(:,m+offset)=Temp(:,m)+nato
          END DO
          offset=offset+q
          nato=nato+natop
       END DO
       DEALLOCATE(Temp)
       WRITE(kprint,*) 'New Angles No. =====>',SIZE(Tpg % Angles,2)
    END IF
    IF(ALLOCATED(Tpg % Dihed)) THEN
       o=SHAPE(Tpg % Dihed)
       ALLOCATE(Temp(o(1), o(2)))
       p=Tpg % S_Dihed
       q=o(2)-Tpg % S_Dihed+1
       Temp = Tpg % Dihed
       DEALLOCATE(Tpg % Dihed)
       ALLOCATE(Tpg % Dihed(o(1), o(2)+(New_Units-1)*q))
       Tpg % Dihed(:,1:o(2))=Temp
       offset=0
       nato=0
       DO n=1,New_Units
          DO m=p,p+q-1
             Tpg % Dihed(:,m+offset)=Temp(:,m)+nato
          END DO
          offset=offset+q
          nato=nato+natop
       END DO
       DEALLOCATE(Temp)
       WRITE(kprint,*) 'New Torsions No. =====>',SIZE(Tpg % Dihed,2)
    END IF
    IF(ALLOCATED(Tpg % Imph)) THEN
       o=SHAPE(Tpg % Imph)
       ALLOCATE(Temp(o(1), o(2)))
       p=Tpg % S_Imph
       q=o(2)-Tpg % S_Imph+1
       Temp = Tpg % Imph
       DEALLOCATE(Tpg % Imph)
       ALLOCATE(Tpg % Imph(o(1), o(2)+(New_Units-1)*q))
       Tpg % Imph(:,1:o(2))=Temp
       offset=0
       nato=0
       DO n=1,New_Units
          DO m=p,p+q-1
             Tpg % Imph(:,m+offset)=Temp(:,m)+nato
          END DO
          offset=offset+q
          nato=nato+natop
       END DO
       DEALLOCATE(Temp)
       WRITE(kprint,*) 'New ITors No. =====>',SIZE(Tpg % imph ,2)
    END IF
    IF(ALLOCATED(Tpg % Int14)) THEN
       o=SHAPE(Tpg % Int14)
       ALLOCATE(Temp(o(1), o(2)))
       p=Tpg % S_Int14
       q=o(2)-Tpg % S_Int14+1
       Temp = Tpg % Int14
       DEALLOCATE(Tpg % Int14)
       ALLOCATE(Tpg % Int14(o(1), o(2)+(New_Units-1)*q))
       Tpg % Int14(:,1:o(2))=Temp
       offset=0
       nato=0
       DO n=1,New_Units
          DO m=p,p+q-1
             Tpg % Int14(:,m+offset)=Temp(:,m)+nato
          END DO
          offset=offset+q
          nato=nato+natop
       END DO
       DEALLOCATE(Temp)
       WRITE(kprint,*) 'New Int14 No. =====>',SIZE(Tpg % Int14,2)
    END IF

    IF(ALLOCATED(Tpg % Mol_Atm)) THEN
       o1=SIZE(Tpg % Mol_Atm)
       ALLOCATE(Temp1(o1))
       p=Tpg % S_Mol_Atm
       q=o1-p+1

       Temp1 = Tpg % Mol_Atm
       DEALLOCATE(Tpg % Mol_Atm)
       ALLOCATE(Tpg % Mol_Atm (o1+(New_Units-1)*q))

       DO n=1,p-1
          s=SIZE(Temp1 (n) % g)
          ALLOCATE(Tpg % Mol_Atm (n) % g (s))
          Tpg % Mol_Atm (n) % g = Temp1 (n) % g 
       END DO

       offset=0
       nato=0
       DO n=1,New_Units
          DO m=p,p+q-1
             s=SIZE(Temp1 (m) % g)
             ALLOCATE(Tpg % Mol_Atm (m + offset) % g (s))
             Tpg % Mol_Atm (m+offset) % g  = Temp1 (m) % g + nato
          END DO
          offset=offset+q
          nato=nato+natop
       END DO       
       DEALLOCATE(Temp1)
       WRITE(kprint,*) 'New Molecule No. =====>',SIZE(Tpg % Mol_Atm)
    END IF

  END SUBROUTINE SystemTpg__Update
