MODULE Class_ResidueTpg

!!$***********************************************************************
!!$   Time-stamp: <2006-12-18 12:23:34 marchi>                           *
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
  USE ERROR_Mod,ONLY: Add_error=>Add, Print_Errors
  USE STRPAK
  USE STRINGS_Mod, ONLY: My_Fxm
  USE MYPARSE_Mod, MY_Parse=>parse
  USE Linked_Int4_D
  IMPLICIT none
  PRIVATE
  PUBLIC :: Init
  INTEGER, PARAMETER :: max_atms=7
  TYPE AtomTpg
     CHARACTER(len=max_atms) :: Res,beta, betab
     INTEGER :: Grp_No,Id_res
     REAL(8) :: chg,mass
     INTEGER, DIMENSION(:), ALLOCATABLE :: cnt
  END TYPE AtomTpg
  TYPE ResidueTpg
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: bonds,angles,&
          &dihed,imph,int14
     TYPE(AtomTpg), DIMENSION(:), ALLOCATABLE  :: atm
  END TYPE ResidueTpg
  CHARACTER(len=max_pars), DIMENSION(:), ALLOCATABLE :: strngs
  CHARACTER(len=max_char_long) :: errmsg_f
CONTAINS
  SUBROUTINE Init(Res_Tpg,Res_i)
    TYPE(ResidueTpg) :: Resa
    CHARACTER(len=max_char) :: res_i
    
    INTEGER :: nato, n, m, i, j,i_F
    LOGICAL, SAVE :: Called=.FALSE.
    
    IF(Called) RETURN
    Called=.TRUE.
    
    i_f=Pick_Res(res_i,App_Char)
    nato=SUM(SIZE(App_Char(i_F) % group(:) % g))
    ALLOCATE(Res_Tpg % atm(nato))
    
    Grp_No
    nato=0
    DO i=1,SIZE(App_Char(i_F) % group)
       Grp_No=Grp_No+1
       jm=SIZE(App_Char(i_F)  % group (i) % g)
       DO j=1,jm
          nato=nato+1
          Res_Tpg % atm(nato) % Res = App_Char(i_F) % Type
          Res_Tpg % atm(nato) % Id_Res = i_F
          Res_Tpg % atm(nato) % Grp_No = Grp_No
          IF(j == 1) Grp_Atm(1,Grp_No)=nato
          CALL My_Parse(App_Char(i_F)  % group (i) % g(j), strngs)
          Res_Tpg % atm(nato) % beta = TRIM(strngs(1))
          Res_Tpg % atm(nato) % betab = TRIM(strngs(2))
          CALL SP_Getnum(strngs(3),Res_Tpg % atm (nato) % chg, iflag)
       END DO
       Grp_Atm(2,Grp_No)=nato
    END DO
    CALL bonds(i_F)
    
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
          WRITE(*,*) n,atm(n) % Res
          WRITE(*,*) n,atm(n) % beta
          WRITE(*,*) n,atm(n) % betab
          STOP
       END IF
    END DO
  END SUBROUTINE Init
  SUBROUTINE bonds(i_F)
    INTEGER :: i_F
    
    INTEGER :: n,m,i,j,count_A,i1,i2      
    LOGICAL :: end_of_list
    CHARACTER(len=max_char) :: label0,label1,lab0
    INTEGER :: new_bond(2)
    INTEGER, DIMENSION(:), ALLOCATABLE :: ind_x
    
!!$
!!$--- Add bonds already in the residue i_F list
!!$
    CALL Start()
    DO i=1,SIZE(App_Char(i_F) % bonds,2)
       DO j=1,SIZE(App_Char(i_F) % bonds,1)
          label0=TRIM(App_Char(i_F) % bonds(j,i))
          IF(ICHAR(label0(1:1)) == 49 .OR. ICHAR(label0(1:1)) == 50) CYCLE
          new_bond(j)=Find_Atom(i_F,Res_No,label0)
       END DO
       CALL Add(new_bond, count_A)
    END DO
!!$
!!$--- Gather all bonds for the residue i_F
!!$              
    
    end_of_list=.FALSE.
    ALLOCATE(Res_Tpg % bonds(2,count_a))
    count_A=0
    CALL Extract(bonds, end_of_list, 0)
    DO WHILE(.NOT. end_of_list) 
       count_A=count_A+1
       Res_Tpg % bonds(1,count_A)=bonds(1)
       Res_Tpg % bonds(2,count_A)=bonds(2)
       CALL Extract(bonds, end_of_list)
    END DO
    CALL CleanUp()
    
    ALLOCATE(ind_x(nato))      
    ind_x=0
    
    DO i=1,SIZE(Tpg % bonds,2)
       i1=Res_Tpg % bonds(1,i)
       i2=Res_Bonds(2,i)
       ind_x(i1)=ind_x(i1)+1
       ind_x(i2)=ind_x(i2)+1
    END DO
    
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
    DEALLOCATE(ind_x)
  END SUBROUTINE Bonds
  SUBROUTINE Bendings(i_F)
    IMPLICIT NONE 
    INTEGER :: count_a, count_out,i,m1,m2,j,k,ii,ll,l,n,ia,iaa,ja&
         &,jaa,ka,kaa,la,laa
    INTEGER :: angles(3)
    LOGICAL :: end_of_list,oka1,oka2,oka3,oka4,okb1,okb2,okb3,okb4
    
    CALL Start()
    
    DO ii=1,SIZE(Res_Tpg % bondss,2)
       i=Tpg % bonds (1,ii)
       j=Tpg % bonds (2,ii)
       m1=0
       IF(ALLOCATED(Res_Tpg % atm(i) % cnt)) THEN
          m1=SIZE(Res_Tpg % atm(i) % cnt)
       END IF
       angles(2)=i
       angles(3)=j
       DO ll=1,m1
          l=Res_Tpg % atm(i) % cnt (ll) 
          angles(1)=l 
          IF(l == j) CYCLE
!!$
!!$--- Eliminate bending in rigid three charge water models
!!$
          CALL Add(angles,count_out)
       END DO
       angles(1)=i
       angles(2)=j
       m2=0
       IF(ALLOCATED(Res_Tpg % atm(j) % cnt)) THEN
          m2=SIZE(Res_Tpg % atm(j) % cnt)
       END IF
       DO ll=1,m2
          l=Res_Tpg % atm(j) % cnt (ll) 
          angles(3)=l
          IF(l == j) CYCLE
!!$
!!$--- Eliminate bending in rigid three charge water models
!!$
          CALL Add(tors,count_out)
       END DO
    END DO
    
    ALLOCATE(t_angles(4,count_out))
    
    end_of_list=.FALSE.
    count_A=0
    CALL Extract(angles, end_of_list, 0)
    DO WHILE(.NOT. end_of_list) 
       count_A=count_A+1
       t_angles(1,count_A)=angles(1)
       t_angles(2,count_A)=angles(2)
       t_angles(3,count_A)=angles(3)
       CALL Extract(tors, end_of_list)
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

!!$------------------------------------------------------------------------

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
    
    ALLOCATE(Tpg % angles(3,count_A))
    count_a=0
    DO i=1,SIZE(t_angles,2)
       IF(t_angles(1,i) /= -1) THEN
          count_a=count_a+1
          Tpg % angles (:,count_a)=t_angles(:,i)
       END IF
    END DO
    DEALLOCATE(t_angles)
    WRITE(*,*) 'Bendings No. =====>',count_A
    
  END SUBROUTINE Bendings
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
  FUNCTION Find_Atom(i_Z,No,label) RESULT(out)
    INTEGER :: No, out, i_Z
    CHARACTER(len=*) :: label
    INTEGER :: n
    CHARACTER(len=max_char) :: label0
    IF((.NOT. ALLOCATED(Res_Atm)) .OR. (.NOT. ALLOCATED(App_Char))) THEN
       errmsg_f='Cannot find atom without initialization'
       CALL Add_Error(-1,errmsg_f)
       CALL Print_Errors()
    END IF
    DO n=Res_Atm(1,No),Res_Atm(2,No)
       IF(TRIM(atm(n) % beta) == TRIM(label)) THEN
          out=n
          RETURN
       END IF
    END DO
    WRITE(label0,'(1x,i4,1x)') No
    errmsg_f='Inter-residue connection of secondary sequence not found. Atom '&
         &//TRIM(label)//' not on residue No. '//TRIM(label0)&
         &//'[ of type '//TRIM(App_Char(i_Z) % Type)//'] &
         &of the secondary sequence.'
    CALL Add_Error(-1,errmsg_f)
    CALL Print_Errors()
  END FUNCTION Find_Atom


!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE Class_ResidueTpg
