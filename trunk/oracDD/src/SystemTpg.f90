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
MODULE SystemTpg

!!$***********************************************************************
!!$   Time-stamp: <2007-01-12 18:49:51 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Tue Nov 28 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program ORAC ----*


  USE Types
  USE IndSequence
  USE SecondarySeq
  USE Parameters
  USE Errors, ONLY: Add_error=>Add,Print_Errors, errmsg_f
  USE AtomCnt, atm_cnt=> AtomCnts
  USE ResidueTpg
  USE Constants
  USE Node
  USE Tops

  IMPLICIT none
  PRIVATE
  PUBLIC :: SystemTpg_, SystemTpg__Type, Tpg, Atom2Cnt, SystemTpg__Update
  INTEGER, PARAMETER :: max_atms=7
  TYPE Atom2Cnt
     TYPE(AtomCnt__Type), POINTER  :: a=>NULL()
  END TYPE Atom2Cnt
  TYPE SystemTpg__Type
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: bonds,angles,&
          &dihed,imph,int14,Grp_Atm,Res_Atm
     INTEGER :: s_bonds,s_angles,s_dihed,s_imph,s_int14,s_Mol_Atm
     TYPE(Atom2Cnt), DIMENSION(:), ALLOCATABLE  :: atm
     TYPE(Chain), DIMENSION(:), ALLOCATABLE :: Mol_Atm
  END TYPE SystemTpg__Type
  TYPE(SystemTpg__Type), SAVE :: Tpg
  CHARACTER(len=max_pars), DIMENSION(:), ALLOCATABLE :: strngs
  INTEGER, DIMENSION(:,:), POINTER :: T_bonds,T_angles,T_tors,T_Int14
  TYPE(ResidueTpg__Type), DIMENSION(:), POINTER :: Res_Used=>NULL()
  INTEGER, DIMENSION(:,:), SAVE, POINTER :: Res_Atm=>NULL()
  INTEGER, DIMENSION(:,:), SAVE, POINTER :: Grp_Atm=>NULL()
CONTAINS
  SUBROUTINE SystemTpg_
    INTEGER :: n,m,level,c_Perc,i_F,count_a
    LOGICAL, SAVE :: Called=.FALSE.

    
    Grp_Atm=>IndSequence__Grp()
    Res_Atm=>IndSequence__Res()

    CALL Used_Residues
    CALL Atoms

    CALL Bonds
    CALL Angles
    CALL Torsions
    CALL ITorsions
    CALL Molecules
    DEALLOCATE(Res_Tpg)
    NULLIFY(Res_Used)
  CONTAINS
    SUBROUTINE Used_Residues
      LOGICAL, DIMENSION(:), ALLOCATABLE :: oks
      INTEGER :: n,m,i_F
      CHARACTER(len=max_char) :: res_i
      ALLOCATE(oks(SIZE(App_Char)))
      oks=.FALSE.
      DO n=1,SIZE(Secondary)
         IF(.NOT. ALLOCATED(Secondary(n) % line)) CYCLE
         DO m=1,SIZE(Secondary(n) % line)
            res_i=Secondary(n) % line(m)
            i_f=IndSequence__PickRes(res_i)
            oks(i_f)=.TRUE.
         END DO
      END DO
      DO n=1,SIZE(App_Char)
         IF(oks(n)) THEN
            res_i=App_Char (n) % Type
            CALL ResidueTpg_(res_i)
         END IF
      END DO
      Res_Used=>Res_Tpg
      DEALLOCATE(oks)
    END SUBROUTINE Used_Residues
    SUBROUTINE Atoms
      INTEGER :: n
      WRITE(*,*)

      WRITE(*,*) 'Atoms No. =====>',SIZE(atm_cnt)
      ALLOCATE(Tpg % atm(SIZE(atm_cnt)))
      ALLOCATE(Tpg % Grp_Atm(SIZE(Grp_Atm,1),SIZE(Grp_Atm,2)))
      ALLOCATE(Tpg % Res_Atm(SIZE(Res_Atm,1),SIZE(Res_Atm,2)))
      
      WRITE(*,*) 'Residue No. =====>',SIZE(Tpg % Res_Atm,2)
      WRITE(*,*) 'Group No. =====>',SIZE(Tpg % Grp_Atm,2)
      Tpg % Grp_Atm = Grp_Atm
      Tpg % Res_Atm = Res_Atm
      DO n=1,SIZE(atm_cnt)
         Tpg % atm (n) % a => atm_cnt(n)
      END DO
    END SUBROUTINE Atoms
    SUBROUTINE Bonds
      LOGICAL ::  end_of_list
      INTEGER :: count_a, count_out,i,m,jj,j,ia,iaa,ja,jaa
      INTEGER :: s_bonds(2)
      INTEGER, POINTER :: p_bonds(:)=>NULL()
      
      IF(.NOT. Node_()) STOP
      DO i=1,SIZE(atm_cnt)
         s_bonds(1)=i
         m=SIZE(atm_cnt(i) % cnt)
         DO jj=1,m
            j=atm_cnt(i) % cnt (jj)
            s_bonds(2)=j
            CALL Node__Push(s_bonds)
         END DO
      END DO
      count_out=Node__Size()
      ALLOCATE(t_bonds(2,count_out))
      count_A=0
      DO WHILE(Node__Pop(p_bonds))
         count_A=count_A+1
         t_bonds(1,count_A)=p_bonds(1)
         t_bonds(2,count_A)=p_bonds(2)
      END DO
      DO i=1,count_A
         ia=t_bonds(1,i)
         ja=t_bonds(2,i)
         DO j=1,i-1
            iaa=t_bonds(1,j)
            jaa=t_bonds(2,j)
            IF(iaa == -1) CYCLE
            IF(iaa == ja .AND. jaa == ia) THEN
               t_bonds(1,i)=-1
               EXIT
            END IF
         END DO
      END DO
      count_a=COUNT(t_bonds(1,:) /= -1)

      ALLOCATE(Tpg % bonds(2,count_A))
      count_a=0
      DO i=1,SIZE(t_bonds,2)
         IF(t_bonds(1,i) /= -1) THEN
            count_a=count_a+1
            Tpg % bonds(:,count_a)=t_bonds(:,i)
         END IF
      END DO
      DEALLOCATE(t_bonds)
      Tpg % s_bonds=get_Slv(Tpg % bonds)
      WRITE(*,*) 'Bonds No. =====>',count_A
      
      IF(ASSOCIATED(p_bonds)) DEALLOCATE(p_bonds)
    END SUBROUTINE Bonds
    SUBROUTINE Angles
      IMPLICIT NONE 
      INTEGER :: count_a, count_out,i,m1,m2,j,jj,kk,k,m3,oo,o,ia,iaa&
           &,ja,jaa,ka,kaa,ii,ll,l,i_f
      INTEGER :: bends(3),nstart,nend,offset
      LOGICAL :: end_of_list,ok
      INTEGER, POINTER :: p_bends(:)=>NULL()
      CHARACTER(len=max_char) :: res_i

      IF(.NOT. Node_()) STOP
      offset=0
      DO n=1,SIZE(Secondary)
         IF(.NOT. ALLOCATED(Secondary(n) % line)) CYCLE
         DO m=1,SIZE(Secondary(n) % line)
            res_i=Secondary(n) % line(m)
            i_f=IndSequence__PickRes(res_i)
            IF(ALLOCATED(Res_Used(i_F) % Angles)) THEN
               DO i=1,SIZE(Res_Used(i_F) % Angles,2)
                  bends(1)=Res_Used(i_F) % Angles (1, i)+offset
                  bends(2)=Res_Used(i_F) % Angles (2, i)+offset
                  bends(3)=Res_Used(i_F) % Angles (3, i)+offset
                  CALL Node__Push(bends)
               END DO
            END IF
            offset=offset+SIZE(Res_Used(i_F) % atm)
         END DO
      END DO
      count_out=Node__Size()
      nstart=count_out+1
      DO ii=1,SIZE(Tpg % bonds,2)
         i=Tpg % bonds (1,ii)
         j=Tpg % bonds (2,ii)
         bends(2)=i
         bends(3)=j
         m1=0
         IF(ALLOCATED(atm_cnt(i) % cnt)) THEN
            m1=SIZE(atm_cnt(i) % cnt)
         END IF
         DO ll=1,m1
            l=atm_cnt(i) % cnt (ll) 
            IF(l == j) CYCLE
            IF((atm_cnt(i) % Res_No == atm_cnt(j) % Res_No) &
                 &.AND. (atm_cnt(l) % Res_No == atm_cnt(j) % Res_No)) CYCLE
            bends(1)=l 
            CALL Node__Push(bends)
         END DO

         bends(1)=i
         bends(2)=j
         m2=0
         IF(ALLOCATED(atm_cnt(j) % cnt)) THEN
            m2=SIZE(atm_cnt(j) % cnt)
         END IF
         DO ll=1,m2
            l=atm_cnt(j) % cnt (ll) 
            IF(l == i) CYCLE
            IF((atm_cnt(i) % Res_No == atm_cnt(j) % Res_No) &
                 &.AND. (atm_cnt(l) % Res_No == atm_cnt(j) % Res_No)) CYCLE
            bends(3)=l
            CALL Node__Push(bends)
         END DO
      END DO
      count_out=Node__Size()
      nend=count_out
      ALLOCATE(t_angles(3,count_out))
      count_A=0
      DO WHILE(Node__Pop(p_bends))
         count_A=count_A+1
         t_angles(1,count_A)=p_bends(1)
         t_angles(2,count_A)=p_bends(2)
         t_angles(3,count_A)=p_bends(3)
      END DO

      DO i=nstart,nend
         ia=t_angles(1,i)
         ja=t_angles(2,i)
         ka=t_angles(3,i)
         DO j=nstart,i-1
            iaa=t_angles(1,j)
            jaa=t_angles(2,j)
            kaa=t_angles(3,j)
            IF(iaa == -1) CYCLE
            IF( (iaa == ia .AND. jaa == ja .AND. kaa == ka) .OR. &
                 &(iaa == ka .AND. jaa == ja .AND. kaa == ia)) THEN
               t_angles(1,i)=-1
               EXIT
            END IF
         END DO
      END DO
      count_a=COUNT(t_angles(1,:) /= -1)

      ALLOCATE(Tpg % angles(3,count_A))
      count_a=0
      DO i=1,SIZE(t_angles,2)
         IF(t_angles(1,i) /= -1) THEN
            count_a=count_a+1
            Tpg % angles (:,count_a)=t_angles(:,i)
         END IF
      END DO
      DEALLOCATE(t_angles)
      WRITE(*,*) 'Angles No. =====>',count_A
      IF(ASSOCIATED(p_bends)) DEALLOCATE(p_bends)
      Tpg % s_angles=get_Slv(Tpg % angles)

    END SUBROUTINE Angles
    SUBROUTINE Torsions
      IMPLICIT NONE 
      INTEGER :: count_a, count_out,i,m1,m2,j,k,ii,ll,l,n,ia,iaa,ja&
           &,jaa,ka,kaa,la,laa,nend,nstart,i_F
      INTEGER :: tors(4),new_Int14(2),k1,k2,k3,offset
      INTEGER, DIMENSION(:), POINTER :: p_tors=>NULL(),p_int14=>NULL()
      LOGICAL :: end_of_list,oka1,oka2,oka3,oka4,okb1,okb2,okb3,okb4,ok,dyn
      CHARACTER(len=max_char) :: res_i
      INTEGER :: mn,m,nstart_it

      IF(.NOT. Node_()) STOP
      offset=0
      DO n=1,SIZE(Secondary)
         IF(.NOT. ALLOCATED(Secondary(n) % line)) CYCLE
         DO m=1,SIZE(Secondary(n) % line)
            res_i=Secondary(n) % line(m)
            i_f=IndSequence__PickRes(res_i)
            IF(ALLOCATED(Res_Used(i_F) % dihed )) THEN
               DO i=1,SIZE(Res_Used(i_F) % dihed,2)
                  tors(1)=Res_Used(i_F) % dihed (1, i)+offset
                  tors(2)=Res_Used(i_F) % dihed (2, i)+offset
                  tors(3)=Res_Used(i_F) % dihed (3, i)+offset
                  tors(4)=Res_Used(i_F) % dihed (4, i)+offset
                  CALL Node__Push(tors)
               END DO
            END IF
            offset=offset+SIZE(Res_Used(i_F) % atm)
         END DO
      END DO
      count_out=Node__Size()
      nstart=count_out+1

      DO ii=1,SIZE(Tpg % angles,2)
         i=Tpg % angles (1,ii)
         j=Tpg % angles (2,ii)
         k=Tpg % angles (3,ii)
         m1=0
         IF(ALLOCATED(atm_cnt(i) % cnt)) THEN
            m1=SIZE(atm_cnt(i) % cnt)
         END IF
         tors(2)=i
         tors(3)=j
         tors(4)=k
         DO ll=1,m1
            l=atm_cnt(i) % cnt (ll) 
            IF(l == j .OR. l == k) CYCLE
            IF((atm_cnt(i) % Res_No == atm_cnt(j) % Res_No) &
              &.AND. (atm_cnt(i) % Res_No  == atm_cnt(k) % Res_No)&
              &.AND. (atm_cnt(i) % Res_No == atm_cnt(l) % Res_No) ) CYCLE
            tors(1)=l 
            CALL Node__Push(tors)
         END DO
         tors(1)=i
         tors(2)=j
         tors(3)=k
         m2=0
         IF(ALLOCATED(atm_cnt(k) % cnt)) THEN
            m2=SIZE(atm_cnt(k) % cnt)
         END IF
         DO ll=1,m2
            l=atm_cnt(k) % cnt (ll) 
            IF(l == j .OR. l == i) CYCLE
            IF((atm_cnt(i) % Res_No == atm_cnt(j) % Res_No) &
              &.AND. (atm_cnt(i) % Res_No  == atm_cnt(k) % Res_No)&
              &.AND. (atm_cnt(i) % Res_No == atm_cnt(l) % Res_No) ) CYCLE
            tors(4)=l
            CALL Node__Push(tors)
         END DO
      END DO

      count_out=Node__Size()
      nend=count_out
      ALLOCATE(t_tors(4,count_out))
      count_A=0
      DO WHILE(Node__Pop(p_tors)) 
         count_A=count_A+1
         t_tors(1,count_A)=p_tors(1)
         t_tors(2,count_A)=p_tors(2)
         t_tors(3,count_A)=p_tors(3)
         t_tors(4,count_A)=p_tors(4)
      END DO

      DO i=nstart,nend
         ia=t_tors(1,i)
         ja=t_tors(2,i)
         ka=t_tors(3,i)
         la=t_tors(4,i)
         DO j=nstart,i-1
            iaa=t_tors(1,j)
            jaa=t_tors(2,j)
            kaa=t_tors(3,j)
            laa=t_tors(4,j)
            IF(iaa == -1) CYCLE
            IF( (iaa == ia .AND. jaa == ja .AND. kaa == ka .AND. laa == la) .OR. &
                 &(iaa == la .AND. jaa == ka .AND. kaa == ja .AND. laa == ia)) THEN
               t_tors(1,i)=-1
               EXIT
            END IF
         END DO
      END DO
      count_a=COUNT(t_tors(1,:) /= -1)

      ALLOCATE(Tpg % dihed(4,count_A))
      count_a=0
      DO i=1,SIZE(t_tors,2)
         IF(t_tors(1,i) /= -1) THEN
            count_a=count_a+1
            Tpg % dihed (:,count_a)=t_tors(:,i)
         END IF
      END DO
      DEALLOCATE(t_tors)
      WRITE(*,*) 'Torsions No. =====>',SIZE(Tpg % dihed,2)
      Tpg % s_dihed=get_Slv(Tpg % dihed)

      IF(.NOT. Node_()) STOP
      offset=0
      DO n=1,SIZE(Secondary)
         IF(.NOT. ALLOCATED(Secondary(n) % line)) CYCLE
         DO m=1,SIZE(Secondary(n) % line)
            res_i=Secondary(n) % line(m)
            i_f=IndSequence__PickRes(res_i)
            IF(ALLOCATED(Res_Used(i_F) % Int14 )) THEN
               DO i=1,SIZE(Res_Used(i_F) % Int14,2)
                  tors(1)=Res_Used(i_F) % Int14 (1, i)+offset
                  tors(2)=Res_Used(i_F) % Int14 (2, i)+offset
                  CALL Node__Push(tors)
               END DO
            END IF
            offset=offset+SIZE(Res_Used(i_F) % atm)
         END DO
      END DO
      nstart_it=Node__Size()
      DO ii=nstart,SIZE(Tpg % dihed,2)
         i=Tpg % dihed (1,ii)
         j=Tpg % dihed (4,ii)
         ok=.TRUE.
         m1=0
         IF(ALLOCATED(Tpg % angles)) THEN
            m1=SIZE(Tpg % angles,2)
         END IF
         DO k=1,m1
            k1=Tpg % angles(1,k)
            k2=Tpg % angles(2,k)
            k3=Tpg % angles(3,k)
!!$
!!$--- eliminate 1-4 separated by an angle
!!$          
            IF((k1 == i .AND. k3 == j) .OR. (k1 == j .AND. k3 == i)) ok=.FALSE.
!!$
!!$--- eliminate 1-4 separated by an bond
!!$          
            IF((k1 == i .AND. k2 == j) .OR. (k1 == j .AND. k2 == i)) ok=.FALSE.
            IF((k2 == i .AND. k3 == j) .OR. (k2 == j .AND. k3 == i)) ok=.FALSE.
         END DO
         IF(ok) THEN
            new_Int14(1)=i
            new_Int14(2)=j
            CALL Node__Push(new_Int14)
         END IF
      END DO
      count_out=Node__Size()
      ALLOCATE(t_Int14(2,count_out))
      count_A=0
      DO WHILE(Node__Pop(p_int14))
         count_A=count_A+1
         t_Int14(1,count_A)=p_Int14(1)
         t_Int14(2,count_A)=p_Int14(2)
      END DO

      mn=0
      DO m=nstart_it+1,count_A
         i=t_Int14(1,m)
         j=t_Int14(2,m)
         ok=.TRUE.
         DO n=nstart_it,m-1
            k1=t_Int14(1,n)
            k2=t_Int14(2,n)
            IF((i == k1 .AND. j == k2) .OR. (i == k2 .AND. j == k1)) ok=.FALSE.
         END DO
         IF(ok) THEN
            mn=mn+1
            t_int14(:,nstart_it+mn)=t_int14(:,m)
         END IF
      END DO
      count_a=nstart_it+mn

      ALLOCATE(Tpg % int14(2,count_A))
      Tpg % Int14 =t_Int14(:,1:count_a)

      Tpg % s_Int14=get_Slv(Tpg % int14)
      DEALLOCATE(t_Int14)
      IF(ASSOCIATED(p_tors)) DEALLOCATE(p_tors)
      IF(ASSOCIATED(p_Int14)) DEALLOCATE(p_Int14)      

      WRITE(*,*) 'Int14 No. =====>',SIZE(Tpg % Int14,2)
    END SUBROUTINE Torsions
    SUBROUTINE ITorsions
      LOGICAL ::  end_of_list
      INTEGER :: count_a, count_out,i,n,m,i_p,i_n,Res_No,i_F,j,ic
      CHARACTER(len=max_char) :: res_i,lab0,res_p,res_n
      INTEGER :: itors(4)
      INTEGER, POINTER :: p_itors(:)=>NULL()
      
      IF(.NOT. Node_()) STOP
      Res_No=0
      DO n=1,SIZE(Secondary)
         IF(.NOT. ALLOCATED(Secondary(n) % line)) CYCLE
         DO m=1,SIZE(Secondary(n) % line)
            Res_No=Res_No+1
            res_i=Secondary(n) % line(m)
            i_p=-1
            i_n=-1
            i_f=IndSequence__PickRes(res_i)
            IF(m /= 1) THEN
               res_p=Secondary(n) % line(m-1)
               i_p=IndSequence__PickRes(res_p)
            END IF
            IF(m /= SIZE(Secondary(n) % line)) THEN
               res_n=Secondary(n) % line(m+1)
               i_N=IndSequence__PickRes(res_n)
            END IF
            IF(.NOT. ALLOCATED(App_Char(i_F) % imph)) CYCLE
            DO i=1,SIZE(App_Char(i_F) % imph, 2)
               DO j=1,4
                  lab0=TRIM(App_Char(i_F) % imph (j,i))
                  IF(lab0(1:1) == '+') THEN
                     IF(i_n == -1) THEN
                        errmsg_f='Improper torsion label contains a ''+'''&
                             &//' whereas we are at the end of the '&
                             &//'secondary sequence.'
                        CALL Add_error(-1,errmsg_f)
                        CALL Print_Errors()
                     END IF
                     DO ic=1,LEN_TRIM(lab0)
                        lab0(ic:ic)=lab0(ic+1:ic+1)
                     END DO
                     itors(j)=AtomCnt__Find(Res_No+1,lab0)
                  ELSE IF( lab0(1:1) == '-') THEN
                     IF(i_p == -1) THEN
                        errmsg_f='Improper torsion label contains a ''-'''&
                             &//' whereas we are at the beginning of the '&
                             &//'secondary sequence.'
                        CALL Add_error(-1,errmsg_f)
                        CALL Print_Errors()
                     END IF
                     DO ic=1,LEN_TRIM(lab0)
                        lab0(ic:ic)=lab0(ic+1:ic+1)
                     END DO
                     itors(j)=AtomCnt__Find(Res_No-1,lab0)
                  ELSE
                     itors(j)=AtomCnt__Find(Res_No,lab0)
                  END IF
               END DO
               CALL Node__Push(itors)
            END DO
         END DO
      END DO
      count_out=Node__Size()
      ALLOCATE(Tpg % imph (4,count_out))

      count_A=0
      DO WHILE(Node__Pop(p_itors))
         count_A=count_A+1
         Tpg % imph (:, count_a) = p_itors
      END DO
      WRITE(*,*) 'ITors No. =====>',count_A
      Tpg % s_imph=get_Slv(Tpg % imph)
      IF(ASSOCIATED(p_itors)) DEALLOCATE(p_itors)
    END SUBROUTINE ITorsions
    FUNCTION Get_Slv(tp) RESULT(out)
      INTEGER :: out
      INTEGER :: tp(:,:)
      INTEGER :: n
      out=SIZE(tp,2)+1
      DO n=1,SIZE(tp,2)
         IF(COUNT(Atm_Cnt(tp(:,n)) % Id_Slv == 2) /= 0) THEN
            out=n
            EXIT
         END IF
      END DO
    END FUNCTION Get_Slv

  END SUBROUTINE SystemTpg_
  SUBROUTINE Molecules
    LOGICAL, DIMENSION(:), ALLOCATABLE, SAVE :: oks,old__oks
    INTEGER, DIMENSION(:), ALLOCATABLE :: mols
    INTEGER, DIMENSION(:), POINTER :: p_mols=>NULL()
    INTEGER :: n,m,nmol,count_out,Mol_Atm, count_a
    LOGICAL :: Dyn=.FALSE.,end_of_list
    
    IF(.NOT. ALLOCATED(atm_cnt)) THEN
       errmsg_f=' List of atomic connections not yet defined'
       CALL Add_error(-1,errmsg_f)
       CALL Print_Errors()
    END IF
    ALLOCATE(oks(SIZE(atm_cnt)),old__oks(SIZE(atm_cnt)))
    oks=.TRUE.
    old__oks=.TRUE.
    nmol=0

    IF(.NOT. Node_()) STOP
    DO n=1,SIZE(atm_cnt)
       IF(oks(n)) THEN
          CALL Next_Connection(n)
          nmol=nmol+1
          Mol_Atm=COUNT(old__oks .NEQV. oks)
          ALLOCATE(mols(Mol_Atm))
          Mol_Atm=0
          DO m=1,SIZE(atm_cnt)
             IF(old__oks(m) .NEQV. oks(m)) THEN
                Mol_Atm=Mol_Atm+1
                mols(Mol_Atm)=m
             END IF
          END DO
          CALL Node__Push(mols)
          old__oks=oks
          DEALLOCATE(mols)
       END IF
    END DO
    count_out=Node__Size()
    ALLOCATE(Tpg % Mol_Atm (count_out))

    count_A=0
    DO WHILE(Node__Pop(p_mols)) 
       count_A=count_A+1
       ALLOCATE(Tpg % Mol_Atm (count_A) % g (SIZE(p_mols)))
       Tpg % Mol_Atm (count_A) % g = p_mols
    END DO
    
    DO n=1,SIZE(Tpg % Mol_Atm)
       IF(COUNT(Atm_Cnt(Tpg % Mol_Atm(n) % g(:)) % Id_Slv == 2) /= 0) THEN
          Tpg % s_Mol_Atm=n
          EXIT
       END IF
    END DO
    IF(ASSOCIATED(p_mols)) DEALLOCATE(p_mols)
    WRITE(*,*) 'Molecule No. =====>',SIZE(Tpg % Mol_Atm)
  CONTAINS
    RECURSIVE SUBROUTINE Next_Connection(ia)
      INTEGER, INTENT(IN) :: ia
      
      INTEGER :: i,ii

      IF(.NOT. oks(ia)) RETURN

      oks(ia)=.FALSE.
      DO ii=1,SIZE(atm_cnt(ia) % cnt)
         i=atm_cnt(ia) % cnt(ii)
         CALL Next_Connection(i)         
      END DO
    END SUBROUTINE Next_Connection
  END SUBROUTINE Molecules
  INCLUDE 'SystemTpg__Update.f90'
  SUBROUTINE SystemTpg__Write
    INTEGER :: n
    WRITE(kbinary) SHAPE(Tpg % bonds),SHAPE(Tpg % angles), SHAPE(Tpg % dihed)&
         &,SHAPE(Tpg % imph),SHAPE(Tpg % int14),SHAPE(Tpg % Grp_atm)&
         &,SHAPE(Tpg % res_atm),SIZE(Tpg %Mol_Atm),Tpg % s_bonds,Tpg % s_angles&
         &,Tpg % s_dihed,Tpg % s_imph,Tpg % s_int14,Tpg % s_Mol_Atm
    WRITE(kbinary) Tpg % bonds,Tpg % angles,Tpg % dihed,Tpg % imph&
         &,Tpg % int14,Tpg % Grp_Atm,Tpg % Res_Atm
    DO n=1,SIZE(Tpg % Mol_atm) 
       WRITE(kbinary) SIZE(Tpg % Mol_Atm(n) % g)
       WRITE(kbinary) Tpg % Mol_Atm(n) % g
    END DO
  END SUBROUTINE SystemTpg__Write
  SUBROUTINE SystemTpg__Read
    INTEGER :: n,s,o_Mol_Atm
    INTEGER, DIMENSION(2) :: o_bonds, o_angles,  o_dihed&
         &, o_imph, o_int14, o_Grp_atm&
         &, o_res_atm
    

    READ(kbinary)  o_bonds, o_angles,  o_dihed&
         &, o_imph, o_int14, o_Grp_atm, o_res_atm,o_Mol_Atm&
         &,Tpg % s_bonds,Tpg % s_angles&
         &,Tpg % s_dihed,Tpg % s_imph,Tpg % s_int14,Tpg % s_Mol_Atm
    ALLOCATE(Tpg % bonds(o_bonds(1),o_bonds(2)))
    ALLOCATE(Tpg % angles(o_angles(1),o_angles(2)))
    ALLOCATE(Tpg % dihed(o_dihed(1),o_dihed(2)))
    ALLOCATE(Tpg % imph(o_imph(1),o_imph(2)))
    ALLOCATE(Tpg % int14(o_int14(1),o_int14(2)))
    ALLOCATE(Tpg % Grp_atm(o_Grp_atm(1),o_Grp_atm(2)))
    ALLOCATE(Tpg % Res_Atm(o_Res_Atm(1),o_Res_Atm(2)))
    ALLOCATE(Tpg % Mol_Atm(o_Mol_Atm))
    READ(kbinary) Tpg % bonds,Tpg % angles,Tpg % dihed,Tpg % imph&
         &,Tpg % int14,Tpg % Grp_Atm,Tpg % Res_Atm
    DO n=1,SIZE(Tpg % Mol_atm)
       READ(kbinary) s; ALLOCATE(Tpg % Mol_Atm(n) % g(s))
       READ(kbinary) Tpg % Mol_Atm(n) % g
    END DO
  END SUBROUTINE SystemTpg__Read

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE SystemTpg
