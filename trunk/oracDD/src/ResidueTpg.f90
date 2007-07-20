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
MODULE ResidueTpg

!!$***********************************************************************
!!$   Time-stamp: <2007-01-10 16:04:11 marchi>                           *
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

  USE AtomCnt
  USE Strings, ONLY: My_Fxm
  USE STRPAK, ONLY: SP_Getnum
  USE Parameters_globals
  USE Errors,ONLY: Add_errors=>Add, Print_Errors, errmsg_f
  USE Myparse
  USE Node
  USE Tops
  USE Print_Defs
  IMPLICIT none
  PRIVATE
  PUBLIC :: ResidueTpg_,Res_Tpg, ResidueTpg__Type

  INTEGER, PARAMETER :: max_atms=7
  TYPE ResidueTpg__Type
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: grp_Atm,bonds,angles,&
          &dihed,int14
     TYPE(AtomCnt__Type), DIMENSION(:), ALLOCATABLE  :: atm
  END TYPE ResidueTpg__Type
  TYPE(ResidueTpg__Type), DIMENSION(:), ALLOCATABLE, TARGET, SAVE :: Res_Tpg
  INTEGER, SAVE :: i_F
CONTAINS
  SUBROUTINE  ResidueTpg_(res_i)
    CHARACTER(len=*) :: res_i
    
    INTEGER :: nato, n, m, i, j,iflag,i_mass,Grp_no,jm,nword
    LOGICAL, SAVE :: Called=.FALSE.
    CHARACTER(len=max_char) :: line
    LOGICAL :: ok
    LOGICAL, SAVE :: First_Time=.TRUE.
    
    IF(first_time) THEN
       ALLOCATE(Res_Tpg(SIZE(App_Char)))
    END IF
    
    i_f=Pick_Res(res_i,App_Char)
    nato=0
    DO i=1,SIZE(App_Char(i_F) % group)
       nato=nato+SIZE(App_Char(i_F) % group(i) % g)
    END DO
    ALLOCATE(Res_Tpg (i_F) % atm(nato))

    Grp_No=SIZE(App_Char(i_F) % group)
    ALLOCATE(Res_Tpg (i_f) % Grp_Atm(2,Grp_No))

    Grp_No=0
    nato=0
    DO i=1,SIZE(App_Char(i_F) % group)
       Grp_No=Grp_No+1
       jm=SIZE(App_Char(i_F)  % group (i) % g)
       DO j=1,jm
          nato=nato+1
          Res_Tpg (i_f) % atm(nato) % Res = App_Char(i_F) % Type
          Res_Tpg (i_f) % atm(nato) % Id_Res = i_F
          Res_Tpg (i_F) % atm(nato) % Grp_No = Grp_No
          IF(j == 1) Res_Tpg (i_F) % Grp_Atm(1,Grp_No)=nato

          nword=MyParse_(App_Char(i_F)  % group (i) % g(j))

          Res_Tpg (i_F) % atm(nato) % beta = TRIM(strngs(1))
          Res_Tpg (i_F) % atm(nato) % betab = TRIM(strngs(2))
          CALL SP_Getnum(strngs(3),Res_Tpg (i_F) % atm (nato) % chg, iflag)
       END DO
       Res_Tpg (i_F) % Grp_Atm(2,Grp_No)=nato
    END DO
    i=Find_Atom(0,' ', Res_Tpg (i_F) % atm)
    IF(First_Time) WRITE(kprint,'(''Computing residue topology ====>'')')
    CALL bonds
    CALL Angles
    CALL Torsions
    CALL Int14

    i_mass=-1
    DO n=1,SIZE(App_Char)
       IF(My_Fxm('mass',App_Char(n) % Type)) THEN
          i_mass=n
          EXIT
       END IF
    END DO
    

    IF(i_mass == -1) THEN
       errmsg_f='Atomic masses not found on Topology file: List of &
            &masses in CHARMM format must be added at the&
            & beginning of any Topology file'
       CALL Add_Errors(-1,errmsg_f)
       CALL Print_Errors()
    END IF
    
    DO n=1,SIZE(Res_Tpg (i_F) % atm)
       line=Res_Tpg (i_F) % atm(n) % betab
       ok=.FALSE.
       DO m=1,SIZE(App_Char(i_mass) % mass,2)
          IF(My_Fxm(TRIM(line), App_Char(i_mass) % mass (1,m))) THEN
             Res_Tpg (i_F) % atm(n) % Id_Type = m
             CALL SP_Getnum(App_Char (i_mass) % mass (2,m)&
                  & ,Res_Tpg (i_F) % atm (n) % mass, iflag)
             ok=.TRUE.
             EXIT
          END IF
       END DO
       
       IF(.NOT. ok) THEN
          WRITE(kprint,*) n,Res_Tpg (i_F) % atm(n) % Res
          WRITE(kprint,*) n,Res_Tpg (i_F) % atm(n) % beta
          WRITE(kprint,*) n,Res_Tpg (i_F) % atm(n) % betab
          STOP
       END IF
    END DO
    First_Time=.FALSE.
  CONTAINS
    SUBROUTINE bonds
      INTEGER :: n,m,i,j,count_A,i1,i2      
      LOGICAL :: end_of_list
      CHARACTER(len=max_char) :: label0,label1,lab0
      INTEGER :: new_bond(2),nato
      INTEGER, DIMENSION(:), ALLOCATABLE :: ind_x
      INTEGER, DIMENSION(:), POINTER :: p_bonds=>NULL()
      LOGICAL :: dyn
      
!!$
!!$--- Add bonds already in the residue i_F list
!!$
      IF(.NOT. Node_()) STOP
      DO i=1,SIZE(App_Char(i_F) % bonds,2)
         DO j=1,SIZE(App_Char(i_F) % bonds,1)
            label0=TRIM(App_Char(i_F) % bonds(j,i))
            IF(ICHAR(label0(1:1)) == 49 .OR. ICHAR(label0(1:1)) == 50) CYCLE
            new_bond(j)=Find_Atom(i_F,label0)
         END DO
         CALL Node__Push(new_bond)
      END DO
!!$
!!$--- Gather all bonds for the residue i_F
!!$              
      
      count_a=Node__Size()
      ALLOCATE(Res_Tpg (i_F) % bonds(2,count_a))
      
      count_A=0
      DO WHILE(Node__Pop(p_bonds))
         count_A=count_A+1
         Res_Tpg (i_F) % bonds(1,count_A)=p_bonds(1)
         Res_Tpg (i_F) % bonds(2,count_A)=p_bonds(2)
      END DO
      
      ALLOCATE(ind_x(SIZE(Res_Tpg (i_F) % atm)))      
      ind_x=0
      
      DO i=1,SIZE(Res_Tpg (i_F) % bonds,2)
         i1=Res_Tpg (i_F) % bonds(1,i)
         i2=Res_Tpg (i_F) % bonds(2,i)
         ind_x(i1)=ind_x(i1)+1
         ind_x(i2)=ind_x(i2)+1
      END DO
      nato=SIZE(Res_Tpg (i_F) % atm)
      DO i=1,nato
         ALLOCATE(Res_Tpg (i_F) % atm(i) % cnt (ind_x(i)))
      END DO
      
      ind_x=0
      DO i=1,SIZE(Res_Tpg (i_F) % bonds,2)
         i1=Res_Tpg (i_F) % bonds(1,i)
         i2=Res_Tpg (i_F) % bonds(2,i)
         ind_x(i1)=ind_x(i1)+1
         ind_x(i2)=ind_x(i2)+1
         Res_Tpg (i_F) % atm(i1) % cnt(ind_x(i1)) = i2
         Res_Tpg (i_F) % atm(i2) % cnt(ind_x(i2)) = i1
      END DO
      WRITE(kprint,'(''.'')',ADVANCE='NO')
      
      DEALLOCATE(ind_x)
      IF(ASSOCIATED(p_bonds)) DEALLOCATE(p_bonds)
    END SUBROUTINE Bonds
    SUBROUTINE Angles
      IMPLICIT NONE 
      INTEGER :: count_a, count_out,i,m1,m2,j,k,ii,ll,l,n,ia,iaa,ja&
           &,jaa,ka,kaa,la,laa,q,qq
      LOGICAL :: end_of_list,oka1,oka2,oka3,oka4,okb1,okb2,okb3,okb4
      INTEGER, DIMENSION(3) :: new_angles
      INTEGER, DIMENSION(:), POINTER :: p_angles=>NULL()
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: t_angles
      LOGICAL :: dyn
      
      IF(.NOT. ALLOCATED(Res_Tpg (i_F) % bonds)) RETURN
      
      IF(.NOT. Node_()) STOP    
      DO ii=1,SIZE(Res_Tpg (i_F) % bonds,2)
         i=Res_Tpg (i_F) % bonds (1,ii)
         j=Res_Tpg (i_F) % bonds (2,ii)
         m1=0
         IF(ALLOCATED(Res_Tpg (i_F) % atm(i) % cnt)) THEN
            m1=SIZE(Res_Tpg (i_F) % atm(i) % cnt)
         END IF
         new_angles(2)=i
         new_angles(3)=j
         DO ll=1,m1
            l=Res_Tpg (i_F) % atm(i) % cnt (ll) 
            IF(l == j) CYCLE
            new_angles(1)=l 
            CALL Node__Push(new_angles)
         END DO
         new_angles(1)=i
         new_angles(2)=j
         m2=0
         IF(ALLOCATED(Res_Tpg (i_F) % atm(j) % cnt)) THEN
            m2=SIZE(Res_Tpg (i_F) % atm(j) % cnt)
         END IF
         DO ll=1,m2
            l=Res_Tpg (i_F) % atm(j) % cnt (ll) 
            IF(l == i) CYCLE
            new_angles(3)=l
            CALL Node__Push(new_angles)
         END DO
      END DO
      count_out=Node__Size()
      ALLOCATE(t_angles(3,count_out))
      
      count_A=0
      DO WHILE(Node__Pop(p_angles)) 
         count_A=count_A+1
         t_angles(1,count_A)=p_angles(1)
         t_angles(2,count_A)=p_angles(2)
         t_angles(3,count_A)=p_angles(3)
      END DO
      DO i=1,count_A
         ia=t_angles(1,i)
         ja=t_angles(2,i)
         ka=t_angles(3,i)
!!$
!!$--- Eliminate bending in rigid three charge water models
!!$
         m1=0
         IF(ALLOCATED(Res_Tpg (i_F) % atm(ia) % cnt)) THEN
            m1=SIZE(Res_Tpg (i_F) % atm(ia) % cnt)
         END IF
         DO qq=1,m1
            q=Res_Tpg (i_F) % atm(ia) % cnt (qq) 
            IF(q == ka) THEN
               t_angles(1,i)=-1
            END IF
         END DO
         DO j=1,i-1
            iaa=t_angles(1,j)
            jaa=t_angles(2,j)
            kaa=t_angles(3,j)
            IF(iaa == -1) CYCLE
            IF( (iaa == ia .AND. jaa == ja .AND. kaa == ka)&
                 & .OR. (iaa == ka .AND. jaa == ja .AND. kaa == ia)) THEN
               t_angles(1,i)=-1
               EXIT
            END IF
         END DO
      END DO
      count_a=COUNT(t_angles(1,:) /= -1)
      
      IF(count_a /= 0) THEN
         ALLOCATE(Res_Tpg (i_F) % angles(3,count_A))    
         count_a=0
         DO i=1,SIZE(t_angles,2)
            IF(t_angles(1,i) /= -1) THEN
               count_a=count_a+1
               Res_Tpg (i_F) % angles (:,count_a)=t_angles(:,i)
            END IF
         END DO
      END IF
      
      DEALLOCATE(t_angles)
      IF(ASSOCIATED(p_angles)) DEALLOCATE(p_angles)
      
      WRITE(kprint,'(''.'')',ADVANCE='NO') 
      
    END SUBROUTINE Angles
    SUBROUTINE Torsions
      IMPLICIT NONE 
      INTEGER :: count_a, count_out,i,m1,m2,j,k,ii,ll,l,n,ia,iaa,ja&
           &,jaa,ka,kaa,la,laa,q,qq
      LOGICAL :: end_of_list,oka1,oka2,oka3,oka4,okb1,okb2,okb3,okb4
      INTEGER, DIMENSION(4) :: new_Tors
      INTEGER, DIMENSION(:), POINTER :: p_Tors=>NULL()
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: t_Tors
      LOGICAL :: dyn
      
      IF(.NOT. ALLOCATED(Res_Tpg (i_F) % Angles)) RETURN
      
      IF(.NOT. Node_()) STOP    
      DO ii=1,SIZE(Res_Tpg (i_F) % angles,2)
         i=Res_Tpg (i_F) % angles (1,ii)
         j=Res_Tpg (i_F) % angles (2,ii)
         k=Res_Tpg (i_F) % angles (3,ii)
         m1=0
         IF(ALLOCATED(Res_Tpg (i_F) % atm(i) % cnt)) THEN
            m1=SIZE(Res_Tpg (i_F) % atm(i) % cnt)
         END IF
         new_Tors(2)=i
         new_Tors(3)=j
         new_Tors(4)=k
         DO ll=1,m1
            l=Res_Tpg (i_F) % atm(i) % cnt (ll) 
            IF(l == j) CYCLE
            new_Tors(1)=l 
            CALL Node__Push(new_Tors)
         END DO
         new_Tors(1)=i
         new_Tors(2)=j
         new_Tors(3)=k
         m2=0
         IF(ALLOCATED(Res_Tpg (i_F) % atm(k) % cnt)) THEN
            m2=SIZE(Res_Tpg (i_F) % atm(k) % cnt)
         END IF
         DO ll=1,m2
            l=Res_Tpg (i_F) % atm(k) % cnt (ll) 
            IF(l == j) CYCLE
            new_Tors(4)=l
            CALL Node__Push(new_Tors)
         END DO
      END DO
      
      count_out=Node__Size()
      ALLOCATE(t_Tors(4,count_out))
      
      count_A=0
      DO WHILE(Node__Pop(p_tors)) 
         count_A=count_A+1
         t_Tors(1,count_A)=p_Tors(1)
         t_Tors(2,count_A)=p_Tors(2)
         t_Tors(3,count_A)=p_Tors(3)
         t_Tors(4,count_A)=p_Tors(4)
      END DO
      
      DO i=1,count_A
         ia=t_Tors(1,i)
         ja=t_Tors(2,i)
         ka=t_Tors(3,i)
         la=t_Tors(4,i)
         DO j=1,i-1
            iaa=t_Tors(1,j)
            jaa=t_Tors(2,j)
            kaa=t_Tors(3,j)
            laa=t_Tors(4,j)
            IF(iaa == -1) CYCLE
            IF( (iaa == ia .AND. jaa == ja .AND. kaa == ka .AND. laa == la) .OR. &
                 &(iaa == la .AND. jaa == ka .AND. kaa == ja .AND. laa == ia)) THEN
               t_tors(1,i)=-1
               EXIT
            END IF
         END DO
      END DO
      count_a=COUNT(t_Tors(1,:) /= -1)
      
      ALLOCATE(Res_Tpg (i_F) % dihed(4,count_A))
      
      count_a=0
      DO i=1,SIZE(t_Tors,2)
         IF(t_Tors(1,i) /= -1) THEN
            count_a=count_a+1
            Res_Tpg (i_F) % dihed (:,count_a)=t_Tors(:,i)
         END IF
      END DO
      DEALLOCATE(t_Tors)
      IF(ASSOCIATED(p_tors)) DEALLOCATE(p_tors)
      WRITE(kprint,'(''.'')',ADVANCE='NO') 
      
    END SUBROUTINE Torsions
    SUBROUTINE Int14
      IMPLICIT NONE 
      INTEGER :: count_a, count_out,i,m1,m2,j,k,ii,ll,l,n,m,k1,k2,k3,mn
      LOGICAL :: end_of_list,oka1,oka2,oka3,oka4,okb1,okb2,okb3,okb4
      INTEGER, DIMENSION(2) :: new_int14
      INTEGER, DIMENSION(:), POINTER :: p_int14=>NULL()
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: t_int14
      LOGICAL :: dyn,ok
      
      IF(.NOT. ALLOCATED(Res_Tpg (i_F) % Dihed)) RETURN
      
      IF(.NOT. Node_()) STOP    
      DO ii=1,SIZE(Res_Tpg (i_F) % Dihed,2)
         i=Res_Tpg (i_F) % dihed (1,ii)
         j=Res_Tpg (i_F) % dihed (4,ii)
         ok=.TRUE.
         m1=0
         IF(ALLOCATED(Res_Tpg (i_F) % angles)) THEN
            m1=SIZE(Res_Tpg (i_F) % angles,2)
         END IF
         DO k=1,m1
            k1=Res_Tpg (i_F) % angles(1,k)
            k2=Res_Tpg (i_F) % angles(2,k)
            k3=Res_Tpg (i_F) % angles(3,k)
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
      ALLOCATE(t_int14(2,count_out))
      
      count_A=0
      DO WHILE(Node__Pop(p_int14))
         count_A=count_A+1
         t_Int14(1,count_A)=p_Int14(1)
         t_Int14(2,count_A)=p_Int14(2)
      END DO

      mn=0
      DO m=1,count_A
         i=t_Int14(1,m)
         j=t_Int14(2,m)
         ok=.TRUE.
         DO n=1,m-1
            k1=t_Int14(1,n)
            k2=t_Int14(2,n)
            IF((i == k1 .AND. j == k2) .OR. (i == k2 .AND. j == k1)) ok=.FALSE.
         END DO
         IF(ok) THEN
            mn=mn+1
            t_int14(:,mn)=t_int14(:,m)
         END IF
      END DO
      ALLOCATE(Res_Tpg (i_F) % int14(2,mn))
      Res_Tpg (i_F) % Int14 =t_Int14(:,1:mn)
      DEALLOCATE(t_Int14)
      IF(ASSOCIATED(p_Int14)) DEALLOCATE(p_Int14)
      
      WRITE(kprint,'(''.'')',ADVANCE='NO') 
    END SUBROUTINE Int14
    FUNCTION Pick_Res(res,Res_C) RESULT (out)
      TYPE(Tops__Type), DIMENSION(:) :: Res_C
      CHARACTER(len=*) :: res
      INTEGER :: i,i_found,out
      
      DO i=1,SIZE(Res_C)
         IF(res  == Res_C(i)%Type) THEN
            i_found=i
            EXIT
         END IF
      END DO
      out=i_found
    END FUNCTION Pick_Res
    FUNCTION Find_Atom(i_Z,label,atm) RESULT(out)
      INTEGER ::  out, i_Z
      CHARACTER(len=*) :: label
      TYPE(AtomCnt__Type), DIMENSION(:), OPTIONAL :: atm
      INTEGER :: n
      CHARACTER(len=max_char) :: label0
      CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE, SAVE :: beta
      
      IF(PRESENT(atm))THEN
         IF(ALLOCATED(beta)) DEALLOCATE(beta)
         ALLOCATE(beta(SIZE(atm)))
         beta=atm(:) % beta
         out=0
         RETURN
      END IF
      DO n=1,SIZE(beta)
         IF(TRIM(beta(n)) == TRIM(label)) THEN
            out=n
            RETURN
         END IF
      END DO
      errmsg_f='Atom label '//TRIM(label)//' not found on residue of type '&
           &//TRIM(App_Char(i_Z) % Type)
      CALL Add_Errors(-1,errmsg_f)
      CALL Print_Errors()
    END FUNCTION Find_Atom
  END SUBROUTINE ResidueTpg_
  
  
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*
  
END MODULE ResidueTpg
