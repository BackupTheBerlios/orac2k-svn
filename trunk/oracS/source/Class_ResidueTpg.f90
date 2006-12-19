MODULE Class_ResidueTpg

!!$***********************************************************************
!!$   Time-stamp: <2006-12-19 14:16:54 marchi>                           *
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


  USE STRINGS_Mod, ONLY: My_Fxm
  USE STRPAK, ONLY: SP_Getnum
  USE PARAMETERS_GLOBALS
  USE TOPOLOGY_GLOBALS
  USE ERROR_Mod,ONLY: Add_error=>Add, Print_Errors
  USE MYPARSE_Mod, MY_Parse=>parse
  USE Linked_Int4_D
  IMPLICIT none
  PRIVATE
  PUBLIC :: Init,ResidueTpg,AtomTpg
  INTEGER, PARAMETER :: max_atms=7
  TYPE AtomTpg
     CHARACTER(len=max_atms) :: Res,beta, betab
     INTEGER :: Grp_No,Id_res,Id_Type
     REAL(8) :: chg,mass
     INTEGER, DIMENSION(:), ALLOCATABLE :: cnt
  END TYPE AtomTpg
  TYPE ResidueTpg
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: grp_Atm,bonds,angles,&
          &dihed,int14
     TYPE(AtomTpg), DIMENSION(:), ALLOCATABLE  :: atm
  END TYPE ResidueTpg
  CHARACTER(len=max_pars), DIMENSION(:), ALLOCATABLE :: strngs
  CHARACTER(len=max_char_long) :: errmsg_f
CONTAINS
  FUNCTION Init(Res_i) RESULT(Res_Tpg)
    TYPE(ResidueTpg) :: Res_Tpg
    CHARACTER(len=max_char) :: res_i
    
    INTEGER :: nato, n, m, i, j,i_F,iflag,i_mass,Grp_no,jm
    LOGICAL, SAVE :: Called=.FALSE.
    CHARACTER(len=max_char) :: line
    LOGICAL :: ok
    LOGICAL, SAVE :: First_Time=.TRUE.
    
    
    i_f=Pick_Res(res_i,App_Char)
    nato=0
    DO i=1,SIZE(App_Char(i_F) % group)
       nato=nato+SIZE(App_Char(i_F) % group(i) % g)
    END DO
    ALLOCATE(Res_Tpg % atm(nato))

    Grp_No=SIZE(App_Char(i_F) % group)
    ALLOCATE(Res_Tpg % Grp_Atm(2,Grp_No))

    Grp_No=0
    nato=0
    DO i=1,SIZE(App_Char(i_F) % group)
       Grp_No=Grp_No+1
       jm=SIZE(App_Char(i_F)  % group (i) % g)
       DO j=1,jm
          nato=nato+1
          Res_Tpg % atm(nato) % Res = App_Char(i_F) % Type
          Res_Tpg % atm(nato) % Id_Res = i_F
          Res_Tpg % atm(nato) % Grp_No = Grp_No
          IF(j == 1) Res_Tpg % Grp_Atm(1,Grp_No)=nato
          CALL My_Parse(App_Char(i_F)  % group (i) % g(j), strngs)
          Res_Tpg % atm(nato) % beta = TRIM(strngs(1))
          Res_Tpg % atm(nato) % betab = TRIM(strngs(2))
          CALL SP_Getnum(strngs(3),Res_Tpg % atm (nato) % chg, iflag)
       END DO
       Res_Tpg % Grp_Atm(2,Grp_No)=nato
    END DO
    i=Find_Atom(0,' ', Res_Tpg % atm)
    IF(First_Time) WRITE(*,'(''Computing residue topology ====>'')')
    CALL bonds(Res_Tpg,i_F)
    CALL Angles(Res_Tpg)
    CALL Torsions(Res_Tpg)
    CALL Int14(Res_Tpg)

    i_mass=-1
    DO n=1,SIZE(App_Char)
       IF(My_Fxm('mass',App_Char(n) % Type)) THEN
          i_mass=n
          EXIT
       END IF
    END DO
    
    IF(i_mass == -1) STOP
    
    DO n=1,SIZE(Res_Tpg % atm)
       line=Res_Tpg % atm(n) % betab
       ok=.FALSE.
       DO m=1,SIZE(App_Char(i_mass) % mass,2)
          IF(My_Fxm(TRIM(line), App_Char(i_mass) % mass (1,m))) THEN
             Res_Tpg % atm(n) % Id_Type = m
             CALL SP_Getnum(App_Char (i_mass) % mass (2,m) ,Res_Tpg % atm (n) % mass, iflag)
             ok=.TRUE.
             EXIT
          END IF
       END DO
       
       IF(.NOT. ok) THEN
          WRITE(*,*) n,Res_Tpg % atm(n) % Res
          WRITE(*,*) n,Res_Tpg % atm(n) % beta
          WRITE(*,*) n,Res_Tpg % atm(n) % betab
          STOP
       END IF
    END DO
    First_Time=.FALSE.
  END FUNCTION Init
  SUBROUTINE Change(Res_Tpg,Pot)
    TYPE(ResidueTpg) :: Res_Tpg
    CHARACTER(len=max_char), DIMENSION(:) :: Pot
    INTEGER :: n,m
    LOGICAL, DIMENSION(:), ALLOCATABLE :: oks

    ALLOCATE(oks(SIZE(Res_Tpg % atm)))
    oks=.TRUE.
    DO n=1,SIZE(Pot)
       DO m=1,SIZE(Res_Tpg % atm)
          IF(.NOT. oks(m)) CYCLE
          IF(My_Fxm(TRIM(Res_Tpg % atm(m) % betab), TRIM(Pot(n))) ) THEN
             oks(m)=.FALSE.
             Res_Tpg % atm(n) % Id_Type = n
          END IF
       END DO
    END DO
    IF(COUNT(oks) /= 0) THEN
       errmsg_f='Atom type unknown. Catastrophic error! Should not be here.'
       CALL Add_Error(-1, errmsg_f)
       CALL Print_Errors()
    END IF
    DEALLOCATE(oks)
  END SUBROUTINE Change
  SUBROUTINE bonds(Res_Tpg,i_F)
    INTEGER :: i_F
    TYPE(ResidueTpg) :: Res_Tpg
    
    INTEGER :: n,m,i,j,count_A,i1,i2      
    LOGICAL :: end_of_list
    CHARACTER(len=max_char) :: label0,label1,lab0
    INTEGER :: new_bond(2),nato
    INTEGER, DIMENSION(:), ALLOCATABLE :: ind_x,p_bonds
    LOGICAL :: dyn
    
!!$
!!$--- Add bonds already in the residue i_F list
!!$
    CALL Start()
    DO i=1,SIZE(App_Char(i_F) % bonds,2)
       DO j=1,SIZE(App_Char(i_F) % bonds,1)
          label0=TRIM(App_Char(i_F) % bonds(j,i))
          IF(ICHAR(label0(1:1)) == 49 .OR. ICHAR(label0(1:1)) == 50) CYCLE
          new_bond(j)=Find_Atom(i_F,label0)
       END DO
       CALL Add(new_bond, count_A)
    END DO
!!$
!!$--- Gather all bonds for the residue i_F
!!$              
    end_of_list=.FALSE.
    ALLOCATE(Res_Tpg % bonds(2,count_a))
    count_A=0
    CALL Extract(p_bonds, end_of_list, dyn, 0)
    DO WHILE(.NOT. end_of_list) 
       count_A=count_A+1
       Res_Tpg % bonds(1,count_A)=p_bonds(1)
       Res_Tpg % bonds(2,count_A)=p_bonds(2)
       CALL Extract(p_bonds, end_of_list, dyn)
    END DO
    CALL CleanUp()
    
    ALLOCATE(ind_x(SIZE(Res_Tpg % atm)))      
    ind_x=0
    
    DO i=1,SIZE(Res_Tpg % bonds,2)
       i1=Res_Tpg % bonds(1,i)
       i2=Res_Tpg % bonds(2,i)
       ind_x(i1)=ind_x(i1)+1
       ind_x(i2)=ind_x(i2)+1
    END DO
    nato=SIZE(Res_Tpg % atm)
    DO i=1,nato
       ALLOCATE(Res_Tpg % atm(i) % cnt (ind_x(i)))
    END DO
    
    ind_x=0
    DO i=1,SIZE(Res_Tpg % bonds,2)
       i1=Res_Tpg % bonds(1,i)
       i2=Res_Tpg % bonds(2,i)
       ind_x(i1)=ind_x(i1)+1
       ind_x(i2)=ind_x(i2)+1
       Res_Tpg % atm(i1) % cnt(ind_x(i1)) = i2
       Res_Tpg % atm(i2) % cnt(ind_x(i2)) = i1
    END DO
    WRITE(*,'(''.'')',ADVANCE='NO')

    DEALLOCATE(ind_x)
    IF(ALLOCATED(p_bonds)) DEALLOCATE(p_bonds)
  END SUBROUTINE Bonds
  SUBROUTINE Angles(Res_Tpg)
    IMPLICIT NONE 
    TYPE(ResidueTpg) :: Res_Tpg
    INTEGER :: count_a, count_out,i,m1,m2,j,k,ii,ll,l,n,ia,iaa,ja&
         &,jaa,ka,kaa,la,laa,q,qq
    LOGICAL :: end_of_list,oka1,oka2,oka3,oka4,okb1,okb2,okb3,okb4
    INTEGER, DIMENSION(3) :: new_angles
    INTEGER, DIMENSION(:), ALLOCATABLE :: p_angles
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: t_angles
    LOGICAL :: dyn
    
    IF(.NOT. ALLOCATED(Res_Tpg % bonds)) RETURN
    CALL Start()
    
    DO ii=1,SIZE(Res_Tpg % bonds,2)
       i=Res_Tpg % bonds (1,ii)
       j=Res_Tpg % bonds (2,ii)
       m1=0
       IF(ALLOCATED(Res_Tpg % atm(i) % cnt)) THEN
          m1=SIZE(Res_Tpg % atm(i) % cnt)
       END IF
       new_angles(2)=i
       new_angles(3)=j
       DO ll=1,m1
          l=Res_Tpg % atm(i) % cnt (ll) 
          IF(l == j) CYCLE
          new_angles(1)=l 
          CALL Add(new_angles,count_out)
       END DO
       new_angles(1)=i
       new_angles(2)=j
       m2=0
       IF(ALLOCATED(Res_Tpg % atm(j) % cnt)) THEN
          m2=SIZE(Res_Tpg % atm(j) % cnt)
       END IF
       DO ll=1,m2
          l=Res_Tpg % atm(j) % cnt (ll) 
          IF(l == i) CYCLE
          new_angles(3)=l
          CALL Add(new_angles,count_out)
       END DO
    END DO

    ALLOCATE(t_angles(3,count_out))
    
    end_of_list=.FALSE.
    count_A=0
    CALL Extract(p_angles, end_of_list, dyn, 0)
    DO WHILE(.NOT. end_of_list) 
       count_A=count_A+1
       t_angles(1,count_A)=p_angles(1)
       t_angles(2,count_A)=p_angles(2)
       t_angles(3,count_A)=p_angles(3)
       CALL Extract(p_angles, end_of_list, dyn)
    END DO
    CALL Cleanup()

    DO i=1,count_A
       ia=t_angles(1,i)
       ja=t_angles(2,i)
       ka=t_angles(3,i)
!!$
!!$--- Eliminate bending in rigid three charge water models
!!$
       m1=0
       IF(ALLOCATED(Res_Tpg % atm(ia) % cnt)) THEN
          m1=SIZE(Res_Tpg % atm(ia) % cnt)
       END IF
       DO qq=1,m1
          q=Res_Tpg % atm(ia) % cnt (qq) 
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
       ALLOCATE(Res_Tpg % angles(3,count_A))    
       count_a=0
       DO i=1,SIZE(t_angles,2)
          IF(t_angles(1,i) /= -1) THEN
             count_a=count_a+1
             Res_Tpg % angles (:,count_a)=t_angles(:,i)
          END IF
       END DO
    END IF
       
    DEALLOCATE(t_angles)
    IF(ALLOCATED(p_angles)) DEALLOCATE(p_angles)
    
    WRITE(*,'(''.'')',ADVANCE='NO') 
    
  END SUBROUTINE Angles
  SUBROUTINE Torsions(Res_Tpg)
    IMPLICIT NONE 
    TYPE(ResidueTpg) :: Res_Tpg
    INTEGER :: count_a, count_out,i,m1,m2,j,k,ii,ll,l,n,ia,iaa,ja&
         &,jaa,ka,kaa,la,laa,q,qq
    LOGICAL :: end_of_list,oka1,oka2,oka3,oka4,okb1,okb2,okb3,okb4
    INTEGER, DIMENSION(4) :: new_Tors
    INTEGER, DIMENSION(:), ALLOCATABLE :: p_Tors
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: t_Tors
    LOGICAL :: dyn
    
    IF(.NOT. ALLOCATED(Res_Tpg % Angles)) RETURN
    CALL Start()
    
    DO ii=1,SIZE(Res_Tpg % angles,2)
       i=Res_Tpg % angles (1,ii)
       j=Res_Tpg % angles (2,ii)
       k=Res_Tpg % angles (3,ii)
       m1=0
       IF(ALLOCATED(Res_Tpg % atm(i) % cnt)) THEN
          m1=SIZE(Res_Tpg % atm(i) % cnt)
       END IF
       new_Tors(2)=i
       new_Tors(3)=j
       new_Tors(4)=k
       DO ll=1,m1
          l=Res_Tpg % atm(i) % cnt (ll) 
          IF(l == j) CYCLE
          new_Tors(1)=l 
          CALL Add(new_Tors,count_out)
       END DO
       new_Tors(1)=i
       new_Tors(2)=j
       new_Tors(3)=k
       m2=0
       IF(ALLOCATED(Res_Tpg % atm(k) % cnt)) THEN
          m2=SIZE(Res_Tpg % atm(k) % cnt)
       END IF
       DO ll=1,m2
          l=Res_Tpg % atm(k) % cnt (ll) 
          IF(l == j) CYCLE
          new_Tors(4)=l
          CALL Add(new_Tors,count_out)
       END DO
    END DO

    ALLOCATE(t_Tors(4,count_out))
    
    end_of_list=.FALSE.
    count_A=0
    CALL Extract(p_Tors, end_of_list, dyn, 0)
    DO WHILE(.NOT. end_of_list) 
       count_A=count_A+1
       t_Tors(1,count_A)=p_Tors(1)
       t_Tors(2,count_A)=p_Tors(2)
       t_Tors(3,count_A)=p_Tors(3)
       t_Tors(4,count_A)=p_Tors(4)
       CALL Extract(p_Tors, end_of_list, Dyn)
    END DO
    CALL Cleanup()

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
    
    ALLOCATE(Res_Tpg % dihed(4,count_A))

    count_a=0
    DO i=1,SIZE(t_Tors,2)
       IF(t_Tors(1,i) /= -1) THEN
          count_a=count_a+1
          Res_Tpg % dihed (:,count_a)=t_Tors(:,i)
       END IF
    END DO
    DEALLOCATE(t_Tors)
    IF(ALLOCATED(p_tors)) DEALLOCATE(p_tors)
    WRITE(*,'(''.'')',ADVANCE='NO') 
   
  END SUBROUTINE Torsions
  SUBROUTINE Int14(Res_Tpg)
    IMPLICIT NONE 
    TYPE(ResidueTpg) :: Res_Tpg
    INTEGER :: count_a, count_out,i,m1,m2,j,k,ii,ll,l,n,k1,k2,k3
    LOGICAL :: end_of_list,oka1,oka2,oka3,oka4,okb1,okb2,okb3,okb4
    INTEGER, DIMENSION(2) :: new_int14
    INTEGER, DIMENSION(:), ALLOCATABLE :: p_int14
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: t_int14
    LOGICAL :: dyn,ok
    
    IF(.NOT. ALLOCATED(Res_Tpg % Dihed)) RETURN

    CALL Start()
    DO ii=1,SIZE(Res_Tpg % Dihed,2)
       i=Res_Tpg % dihed (1,ii)
       j=Res_Tpg % dihed (4,ii)
       ok=.TRUE.
       m1=0
       IF(ALLOCATED(Res_Tpg % angles)) THEN
          m1=SIZE(Res_Tpg % angles,2)
       END IF
       DO k=1,m1
          k1=Res_Tpg % angles(1,k)
          k2=Res_Tpg % angles(2,k)
          k3=Res_Tpg % angles(3,k)
!!$
!!$--- eliminate 1-4 separated by an angle
!!$          
          IF((k1 == i .AND. k3 == j) .OR. (k1 == j .AND. k3 == i)) ok=.FALSE.
!!$
!!$--- eliminate 1-4 separated by an bond
!!$          
          IF((k1 == i .AND. k2 == j) .OR. (k2 == j .AND. k1 == i)) ok=.FALSE.
          IF((k2 == i .AND. k3 == j) .OR. (k3 == j .AND. k2 == i)) ok=.FALSE.
       END DO
       IF(ok) THEN
          new_Int14(1)=i
          new_Int14(2)=j
          CALL Add(new_Int14,count_out)
       END IF
    END DO

    ALLOCATE(t_int14(2,count_out))
    
    end_of_list=.FALSE.
    count_A=0
    CALL Extract(p_Int14, end_of_list, dyn, 0)
    DO WHILE(.NOT. end_of_list) 
       count_A=count_A+1
       t_Int14(1,count_A)=p_Int14(1)
       t_Int14(2,count_A)=p_Int14(2)
       CALL Extract(p_Int14, end_of_list, dyn)
    END DO
    CALL Cleanup()

    ALLOCATE(Res_Tpg % int14(2,count_A))
    Res_Tpg % Int14 =t_Int14
    DEALLOCATE(t_Int14)
    IF(ALLOCATED(p_Int14)) DEALLOCATE(p_Int14)

    WRITE(*,'(''.'')',ADVANCE='NO') 
  END SUBROUTINE Int14
  FUNCTION Pick_Res(res,Res_C) RESULT (out)
    TYPE(Unit_Char), DIMENSION(:) :: Res_C
    CHARACTER(len=max_char) :: res
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
    TYPE(AtomTpg), DIMENSION(:), OPTIONAL :: atm
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
    CALL Add_Error(-1,errmsg_f)
    CALL Print_Errors()
  END FUNCTION Find_Atom


!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE Class_ResidueTpg
