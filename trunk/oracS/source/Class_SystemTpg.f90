MODULE Class_SystemTpg

!!$***********************************************************************
!!$   Time-stamp: <2006-12-14 18:32:03 marchi>                           *
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

  
  USE PARAMETERS_GLOBALS
  USE TOPOLOGY_GLOBALS
  USE ERROR_Mod,ONLY: Add_error=>Add,Print_Errors
  USE Class_AtomCnt, ONLY: AtomCnt, Atom_cnt=>Atoms, atm_cnt=>atm&
       &, Pick_Res, Find_Atom, AtomCnt__Grp_Atm=>Grp_Atm, &
       & AtomCnt__Res_Atm=>Res_Atm
  USE CONSTANTS
  USE Linked_Int4_D
  IMPLICIT none
  PRIVATE
  PUBLIC :: Init, Atom2Cnt, SystemTpg, Tpg
  INTEGER, PARAMETER :: max_atms=7
  TYPE Atom2Cnt
     TYPE(AtomCnt), POINTER  :: a=>NULL()
  END TYPE Atom2Cnt
  TYPE SystemTpg
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: bonds,angles,&
          &dihed,imph,int14,Grp_Atm,Res_Atm
     TYPE(Atom2Cnt), DIMENSION(:), ALLOCATABLE  :: atm
     TYPE(Chain), DIMENSION(:), ALLOCATABLE :: Mol_Atm
  END TYPE SystemTpg
  TYPE(SystemTpg), SAVE :: Tpg
  CHARACTER(len=max_pars), DIMENSION(:), ALLOCATABLE :: strngs
  INTEGER, DIMENSION(:,:), POINTER :: T_bonds,T_angles,T_tors
  CHARACTER(len=max_char_long) :: errmsg_f
CONTAINS
  SUBROUTINE Init
    INTEGER :: n,m,level,c_Perc
    LOGICAL, SAVE :: Called=.FALSE.

    IF(Called) RETURN
    Called=.TRUE.

    CALL Atom_Cnt
    WRITE(*,*) 'Atoms No. =====>',SIZE(atm_cnt)

    ALLOCATE(Tpg % atm(SIZE(atm_cnt)))
    ALLOCATE(Tpg % Grp_Atm(SIZE(AtomCnt__Grp_Atm,1),SIZE(AtomCnt__Grp_Atm,2)))
    ALLOCATE(Tpg % Res_Atm(SIZE(AtomCnt__Res_Atm,1),SIZE(AtomCnt__Res_Atm,2)))

    WRITE(*,*) 'Residue No. =====>',SIZE(Tpg % Res_Atm,2)
    WRITE(*,*) 'Group No. =====>',SIZE(Tpg % Grp_Atm,2)
    Tpg % Grp_Atm = AtomCnt__Grp_Atm
    Tpg % Res_Atm = AtomCnt__Res_Atm
    DO n=1,SIZE(atm_cnt)
       Tpg % atm (n) % a => atm_cnt(n)
    END DO

    CALL Bonds
    CALL Angles
    CALL Torsions
    CALL ITorsions
    CALL Molecules
  CONTAINS
    SUBROUTINE Bonds
      LOGICAL ::  end_of_list
      INTEGER :: count_a, count_out,i,m,jj,j,ia,iaa,ja,jaa
      INTEGER :: p_bonds(2)
      
      CALL Start()
      DO i=1,SIZE(atm_cnt)
         p_bonds(1)=i
         m=SIZE(atm_cnt(i) % cnt)
         DO jj=1,m
            j=atm_cnt(i) % cnt (jj)
            p_bonds(2)=j
            CALL Add(p_bonds, count_out)
         END DO
      END DO
      ALLOCATE(t_bonds(2,count_out))
      end_of_list=.FALSE.
      count_A=0
      CALL Extract(p_bonds, end_of_list, 0)
      DO WHILE(.NOT. end_of_list) 
         count_A=count_A+1
         t_bonds(1,count_A)=p_bonds(1)
         t_bonds(2,count_A)=p_bonds(2)
         CALL Extract(p_bonds, end_of_list)
      END DO
      CALL Cleanup()
      
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
      WRITE(*,*) 'Bonds No. =====>',count_A
    END SUBROUTINE Bonds
    SUBROUTINE Angles
      IMPLICIT NONE 
      INTEGER :: count_a, count_out,i,m1,m2,j,jj,kk,k,m3,oo,o,ia,iaa&
           &,ja,jaa,ka,kaa
      INTEGER :: bends(3)
      LOGICAL :: end_of_list,ok

      CALL Start()

      DO i=1,SIZE(atm_cnt)
         m1=SIZE(atm_cnt(i) % cnt)
         DO jj=1,m1
            j=atm_cnt(i) % cnt (jj) 
            m2=SIZE(atm_cnt(j) % cnt)
            DO kk=1,m2
               k=atm_cnt(j) % cnt (kk) 
               IF(k == i) CYCLE
               bends(1)=i
               bends(2)=j
               bends(3)=k
               ok=.TRUE.
               m3=SIZE(atm_cnt(k) % cnt)
               DO oo=1,m3
                  o=atm_cnt(k) % cnt(oo)
                  IF(o == i) ok=.FALSE.
               END DO
               IF(i /= k .AND. ok) CALL Add(bends,count_out)
            END DO
         END DO
      END DO

      ALLOCATE(t_angles(3,count_out))
      end_of_list=.FALSE.
      count_A=0
      CALL Extract(bends, end_of_list, 0)
      DO WHILE(.NOT. end_of_list) 
         count_A=count_A+1
         t_angles(1,count_A)=bends(1)
         t_angles(2,count_A)=bends(2)
         t_angles(3,count_A)=bends(3)
         CALL Extract(bends, end_of_list)
      END DO
      CALL Cleanup()

      DO i=1,count_A
         ia=t_angles(1,i)
         ja=t_angles(2,i)
         ka=t_angles(3,i)
         DO j=1,i-1
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

    END SUBROUTINE Angles
    SUBROUTINE Torsions
      IMPLICIT NONE 
      INTEGER :: count_a, count_out,i,m1,m2,j,k,ii,ll,l,n,ia,iaa,ja&
           &,jaa,ka,kaa,la,laa
      INTEGER :: tors(4)
      LOGICAL :: end_of_list,oka1,oka2,oka3,oka4,okb1,okb2,okb3,okb4

      CALL Start()

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
            tors(1)=l 
            IF(l == j .OR. l == k) CYCLE
            CALL Add(tors,count_out)
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
            tors(4)=l
            IF(l == j .OR. l == i) CYCLE
            CALL Add(tors,count_out)
         END DO
      END DO

      ALLOCATE(t_tors(4,count_out))

      end_of_list=.FALSE.
      count_A=0
      CALL Extract(tors, end_of_list, 0)
      DO WHILE(.NOT. end_of_list) 
         count_A=count_A+1
         t_tors(1,count_A)=tors(1)
         t_tors(2,count_A)=tors(2)
         t_tors(3,count_A)=tors(3)
         t_tors(4,count_A)=tors(4)
         CALL Extract(tors, end_of_list)
      END DO
      CALL Cleanup()


      DO i=1,count_A
         ia=t_tors(1,i)
         ja=t_tors(2,i)
         ka=t_tors(3,i)
         la=t_tors(4,i)
         DO j=1,i-1
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
      WRITE(*,*) 'Torsions No. =====>',count_A

    END SUBROUTINE Torsions
    SUBROUTINE ITorsions
      LOGICAL ::  end_of_list
      INTEGER :: count_a, count_out,i,n,m,i_p,i_n,Res_No,i_F,j,ic
      CHARACTER(len=max_char) :: res_i,lab0,res_p,res_n
      INTEGER :: p_itors(4)
      
      CALL Start()

      
      Res_No=0
      DO n=1,SIZE(Secondary_Seq)
         DO m=1,SIZE(Secondary_Seq(n) % line)
            Res_No=Res_No+1
            res_i=Secondary_Seq(n) % line(m)
            i_p=-1
            i_n=-1
            i_f=Pick_Res(res_i,App_Char)
            IF(m /= 1) THEN
               res_p=Secondary_Seq(n) % line(m-1)
               i_p=Pick_Res(res_p,App_Char)
            END IF
            IF(m /= SIZE(Secondary_Seq(n) % line)) THEN
               res_n=Secondary_Seq(n) % line(m+1)
               i_N=Pick_Res(res_n,App_Char)
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
                     p_itors(j)=Find_Atom(i_n,Res_No+1,lab0)
                  ELSE IF( lab0(1:1) == '-') THEN
                     IF(i_n == -1) THEN
                        errmsg_f='Improper torsion label contains a ''-'''&
                             &//' whereas we are at the begiining of the '&
                             &//'secondary sequence.'
                        CALL Add_error(-1,errmsg_f)
                        CALL Print_Errors()
                     END IF
                     DO ic=1,LEN_TRIM(lab0)
                        lab0(ic:ic)=lab0(ic+1:ic+1)
                     END DO
                     p_itors(j)=Find_Atom(i_p,Res_No-1,lab0)
                  ELSE
                     p_itors(j)=Find_Atom(i_F,Res_No,lab0)
                  END IF
               END DO
               CALL Add(p_itors,count_out)
            END DO
         END DO
      END DO
      ALLOCATE(Tpg % imph (4,count_out))

      end_of_list=.FALSE.
      count_A=0
      CALL Extract(p_itors, end_of_list, 0)
      DO WHILE(.NOT. end_of_list) 
         count_A=count_A+1
         Tpg % imph (:, count_a) = p_itors
         CALL Extract(p_itors, end_of_list)
      END DO
      CALL Cleanup()
      WRITE(*,*) 'ITors No. =====>',count_A
    END SUBROUTINE ITorsions
  END SUBROUTINE Init
  SUBROUTINE Molecules
    LOGICAL, DIMENSION(:), ALLOCATABLE, SAVE :: oks,old__oks
    INTEGER, DIMENSION(:), ALLOCATABLE :: mols
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
    CALL Start()
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
          CALL Add(mols,count_out)
          old__oks=oks
          DEALLOCATE(mols)
       END IF
    END DO
    ALLOCATE(Tpg % Mol_Atm (count_out))

    end_of_list=.FALSE.
    count_A=0
    CALL Extract(mols, end_of_list, Dyn, 0)
    DO WHILE(.NOT. end_of_list) 
       count_A=count_A+1
       ALLOCATE(Tpg % Mol_Atm (count_A) % g (SIZE(mols)))
       Tpg % Mol_Atm (count_A) % g = mols
       CALL Extract(mols, end_of_list, Dyn)
    END DO
    CALL Cleanup()
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

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE Class_SystemTpg
