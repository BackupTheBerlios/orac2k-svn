MODULE Class_AtomCnt

!!$***********************************************************************
!!$   Time-stamp: <2006-12-19 13:38:44 marchi>                           *
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
  PUBLIC :: Atoms,AtomCnt,atm, Res_Atm, Grp_atm, Find_Atom, Pick_Res
  INTEGER, PARAMETER :: max_atms=7
  TYPE AtomCnt
     CHARACTER(len=max_atms) :: Res,beta, betab
     INTEGER :: Res_No, Grp_No, Id_Res, Id_Type
     REAL(8) :: chg,mass
     INTEGER, DIMENSION(:), ALLOCATABLE :: cnt
  END TYPE AtomCnt
  TYPE(AtomCnt), DIMENSION(:), ALLOCATABLE, TARGET, SAVE :: atm
  INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE :: Res_Atm
  INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE :: Grp_Atm
  CHARACTER(len=max_pars), DIMENSION(:), ALLOCATABLE :: strngs
  CHARACTER(len=max_char_long) :: errmsg_f
CONTAINS
  SUBROUTINE Atoms
    TYPE(Unit_Char), POINTER :: Res_Inst
    CHARACTER(len=max_char) :: res_i,line
    INTEGER :: nato, n, m, i, j, Res_No, ii, iflag, Grp_No,i_F,jm,i_mass
    TYPE index
       INTEGER, DIMENSION(:), POINTER :: i
    END TYPE index
    TYPE(index), DIMENSION(2) :: inds
    LOGICAL :: ok
    LOGICAL, SAVE :: Called=.FALSE.
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: nind_ptch
    CHARACTER(len=max_char), DIMENSION(:,:), ALLOCATABLE :: ind_ptch

    IF(Called) RETURN
    Called=.TRUE.

    nato=0
    Res_No=0
    Grp_No=0
    ALLOCATE(inds(1) % i (SIZE(Secondary_Seq(1) % line)))
    ALLOCATE(inds(2) % i (SIZE(Secondary_Seq(2) % line)))
    m=0
    DO n=1,SIZE(patches)
       IF(patches (n) % Type == 'link') THEN
          m=m+1
       END IF
    END DO
    ALLOCATE(ind_ptch(2,m))
    ALLOCATE(nind_ptch(2,m))
    m=0
    DO n=1,SIZE(patches)
       IF(patches (n) % Type == 'link') THEN
          m=m+1
          nind_ptch(1,m)=patches(n)%one
          nind_ptch(2,m)=patches(n)%two
          ind_Ptch(1,m)='Link '//TRIM(patches(n)%Res_l(1))
          ind_Ptch(2,m)='Link '//TRIM(patches(n)%Res_l(2))
       END IF
    END DO
    DO n=1,SIZE(Secondary_Seq)
       DO m=1,SIZE(Secondary_Seq(n) % line)
          Res_No=Res_No+1
          IF(n == 1 .AND. ALLOCATED(ind_ptch)) THEN
             res_i=Secondary_Seq(n) % line(m)
             ok=.FALSE.
             DO i=1,SIZE(ind_ptch,2)
                DO j=1,2
                   IF(nind_ptch(j,i) == m) THEN
                      res_i=ind_ptch(j,i)
                      Secondary_Seq(n) % line(m)=TRIM(res_i)
                      ok=.TRUE.
                   END IF
                END DO
             END DO
          END IF
          i_f=Pick_Res(res_i,App_Char)
          inds(n) % i (m) = i_f
!!$
!!$--- Count atoms
!!$
          IF(inds(n) % i (m) == -1) THEN
             errmsg_f='Resdidue '//TRIM(res_i)//&
                  &' not found in force field database'
             CALL Add_Error(-1,errmsg_f)
             RETURN
          END IF
          DO i=1,SIZE(App_Char(i_F) % group)
             Grp_No=Grp_No+1
             nato=nato+SIZE(App_Char(i_F)  % group (i) % g)
          END DO
       END DO
    END DO
    ALLOCATE(atm(nato))
    ALLOCATE(Res_Atm(2,Res_No))
    ALLOCATE(Grp_Atm(2,Grp_No))

    nato=0
    Res_No=0
    Grp_No=0
    DO n=1,SIZE(Secondary_Seq)
       DO m=1,SIZE(Secondary_Seq(n) % line)

          Res_No=Res_No+1
          
          i_F=inds(n) % i (m) 
!!$
!!$--- Count atoms
!!$
          DO i=1,SIZE(App_Char(i_F) % group)
             Grp_No=Grp_No+1
             jm=SIZE(App_Char(i_F)  % group (i) % g)
             DO j=1,jm
                nato=nato+1
                atm(nato) % Res = App_Char(i_F) % Type
                atm(nato) % Res_No = Res_No
                atm(nato) % Id_Res = i_F
                atm(nato) % Grp_No = Grp_No
                IF(i == 1 .AND. j == 1) Res_Atm(1,Res_No)=nato
                IF(j == 1) Grp_Atm(1,Grp_No)=nato
                CALL My_Parse(App_Char(i_F)  % group (i) % g(j), strngs)
                atm(nato) % beta = TRIM(strngs(1))
                atm(nato) % betab = TRIM(strngs(2))
                CALL SP_Getnum(strngs(3),atm (nato) % chg, iflag)
             END DO
             Grp_Atm(2,Grp_No)=nato
          END DO
          Res_Atm(2,Res_No)=nato
       END DO
    END DO

    CALL Get_Connections

    i_mass=-1
    DO n=1,SIZE(App_Char)
       IF(My_Fxm('mass',App_Char(n) % Type)) THEN
          i_mass=n
          EXIT
       END IF
    END DO

    IF(i_mass == -1) STOP

    DO n=1,SIZE(atm)
       line=atm(n) % betab
       ok=.FALSE.
       DO m=1,SIZE(App_Char(i_mass) % mass,2)
          IF(My_Fxm(TRIM(line), App_Char(i_mass) % mass (1,m))) THEN
             atm(n) % Id_Type = m
             CALL SP_Getnum(App_Char (i_mass) % mass (2,m) ,atm (n) % mass, iflag)
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

    DEALLOCATE(ind_ptch)
    DEALLOCATE(nind_ptch)
    CONTAINS
      SUBROUTINE Get_Connections
        INTEGER :: Res_No,n,m,i,j,i_N,i_P,i_F,i1,i2,ii1,ii2,count_a,count_extra,m1,o,o1
        INTEGER :: extra_Bond(2),new_bond(2),bonds(2),nind_x,n1,m2,i_2
        
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: Res_Bonds
        INTEGER, DIMENSION(:), ALLOCATABLE :: ind_x
        CHARACTER(len=max_char) :: label0,label1,lab0
        LOGICAL :: end_of_list,ok
        
        Res_No=0
        DO n=1,SIZE(Secondary_Seq)
           DO m=1,SIZE(Secondary_Seq(n) % line)
              Res_No=Res_No+1
              i_P=-1
              i_N=-1
              IF(m /= 1) i_P=inds(n) % i (m-1)
              IF(m /= SIZE(Secondary_Seq(n) % line)) i_N=inds(n) % i (m+1)
              i_F=inds(n) % i (m) 
!!$
!!$--- Add bonds already in the residue i_F list
!!$
              CALL Start()
              DO i=1,SIZE(App_Char(i_F) % bonds,2)
                 DO j=1,SIZE(App_Char(i_F) % bonds,1)
                    label0=TRIM(App_Char(i_F) % bonds(j,i))
                    IF(ICHAR(label0(1:1)) == 49 .OR. ICHAR(label0(1:1)) == 50) EXIT
                    new_bond(j)=Find_Atom(i_F,Res_No,label0)
                 END DO
                 label0=TRIM(App_Char(i_F) % bonds(1,i))
                 label1=TRIM(App_Char(i_F) % bonds(2,i))
                 IF((ICHAR(label0(1:1)) == 49 .AND. ICHAR(label1(1:1)) == 50)&
                         & .OR. (ICHAR(label0(1:1)) == 50 &
                         &.AND. ICHAR(label1(1:1)) == 49)) THEN
                    m2=-1
                    DO n1=1,SIZE(nind_ptch,2)
                       IF(nind_ptch(1,n1) == m) THEN
                          m2=nind_ptch(2,n1)
                       END IF
                    END DO
                    m1=-1
                    DO n1=1,SIZE(nind_ptch,2)
                       IF(nind_ptch(2,n1) == m) THEN
                          m1=nind_ptch(1,n1)
                       END IF
                    END DO
                    IF(m2 /= -1) THEN
                       DO j=1,SIZE(App_Char(i_F) % bonds,1)
                          label0=TRIM(App_Char(i_F) % bonds(j,i))
                          label1=label0(2:)
                          IF(ICHAR(label0(1:1)) == 49) THEN
                             new_bond(j)=Find_Atom(i_F,Res_No,label1)
                          ELSE IF(ICHAR(label0(1:1)) == 50) THEN
                             i_2=inds(n) % i (m2) 
                             new_bond(j)=Find_Atom(i_2,m2,label1)
                          END IF
                       END DO
                       CALL Add(new_bond, count_A)
                    ELSE IF(m1 == -1) THEN
                       WRITE(lab0,'(i4)') m
                       label0=TRIM(App_Char(i_F) % bonds(1,i))
                       label1=TRIM(App_Char(i_F) % bonds(2,i))
                       errmsg_f='Inter-residue connection '//TRIM(label0)&
                            &//' '//TRIM(label1)//' not found. Residue No. '&
                            &//TRIM(lab0)//' not on the list of linked residues.'
                       CALL Add_Error(-1,errmsg_f)
                       CALL Print_Errors()
                    END IF
                 ELSE
                    CALL Add(new_bond, count_A)
                 END IF
              END DO
!!$
!!$--- Add extra bonds if connections to other residues exist
!!$
              count_extra=0
              IF(ALLOCATED(App_Char(i_F)%ends)) THEN
                 IF((.NOT. My_Fxm('*',TRIM(App_Char(i_F)%ends(1,1))))&
                      & .AND. i_P == -1) THEN 
                    errmsg_f='While creating atom topology, residue '//&
                         &TRIM(App_Char(i_F) % Type)//' appers to be a&
                         &t the beginning of & 
                         &the sequence, but a dangling bond was found.'
                    CALL Add_error(-1,errmsg_f)
                    CALL Print_Errors()
                 ELSE IF((.NOT. My_Fxm('*',TRIM(App_Char(i_F)%ends(2&
                      &,1)))) .AND. i_N == -1) THEN 
                    errmsg_f='While creating atom topology, residue '//&
                         &TRIM(App_Char(i_F) % Type)//' appers to be at the end of &
                         &the sequence, but a dangling bond was found.'
                    CALL Add_error(-1,errmsg_f)
                    CALL Print_Errors()
                 END IF
                 IF(.NOT. My_Fxm('*',TRIM(App_Char(i_F)%ends(2,1)))) THEN
                    extra_bond(1)=Find_Atom(i_F,Res_No,App_Char(i_F) % ends(2,1))
                    extra_bond(2)=Find_Atom(i_N,Res_No+1,App_Char(i_N) % ends(1,1))
                    count_extra=1
                 END IF
              END IF
!!$
!!$--- Gather all bonds for the residue i_F
!!$              

              end_of_list=.FALSE.
              ALLOCATE(Res_bonds(2,count_a+Count_Extra))
              count_A=0
              CALL Extract(bonds, end_of_list, 0)
              DO WHILE(.NOT. end_of_list) 
                 count_A=count_A+1
                 Res_Bonds(1,count_A)=bonds(1)
                 Res_Bonds(2,count_A)=bonds(2)
                 CALL Extract(bonds, end_of_list)
              END DO
              IF(count_Extra /= 0) THEN
                 Res_bonds(1,count_A+1)=extra_bond(1)
                 Res_bonds(2,count_A+1)=extra_bond(2)
              END IF
              CALL CleanUp()
              nind_x=Res_Atm(2,Res_No)-Res_Atm(1,Res_No)+1
              ALLOCATE(ind_x(nind_x))      
              ind_x=0
              
              DO i=1,SIZE(Res_Bonds,2)
                 i1=Res_Bonds(1,i)
                 ii1=i1-Res_Atm(1,Res_No)+1
                 i2=Res_Bonds(2,i)
                 ii2=i2-Res_Atm(1,Res_No)+1
                 IF(ii1 <= nind_x) ind_x(ii1)=ind_x(ii1)+1
                 IF(ii2 <= nind_x) ind_x(ii2)=ind_x(ii2)+1
              END DO
              DO i=Res_Atm(1,Res_No),Res_Atm(2,Res_No)
                 ii1=i-Res_Atm(1,Res_No)+1
                 IF(ii1 <= nind_x) ALLOCATE(atm(i) % cnt (ind_x(ii1)))
              END DO
              ind_x=0
              DO i=1,SIZE(Res_Bonds,2)
                 i1=Res_Bonds(1,i)
                 ii1=i1-Res_Atm(1,Res_No)+1
                 i2=Res_Bonds(2,i)
                 ii2=i2-Res_Atm(1,Res_No)+1
                 IF(ii1 <= nind_x) THEN
                    ind_x(ii1)=ind_x(ii1)+1
                    atm(i1) % cnt(ind_x(ii1)) = i2
                 END IF
                 IF(ii2 <= nind_x) THEN
                    ind_x(ii2)=ind_x(ii2)+1
                    atm(i2) % cnt(ind_x(ii2)) = i1
                 END IF
              END DO
              DEALLOCATE(ind_x)
              DEALLOCATE(Res_bonds)
           END DO
        END DO
!!$
!!$--- Add additional links between residues
!!$



!!$
!!$--- Then add missing connections: So far, the bond between 
!!$--- the n and n+1 residues is counted only for the n-th residue
!!$
        DO n=1,SIZE(atm)
           DO m=1,SIZE(atm (n) % cnt)
              m1=atm (n) % cnt (m) 
              ok=.FALSE.
              DO o=1,SIZE(atm (m1) % cnt)
                 o1=atm (m1) % cnt (o) 
                 IF(o1 == n) THEN
                    ok=.TRUE.
                    EXIT
                 END IF
              END DO
              IF(.NOT. ok) THEN
                 o1=SIZE(atm (m1) % cnt)
                 ALLOCATE(ind_x(o1+1))
                 ind_x(1:o1)=atm (m1) % cnt
                 ind_x(o1+1)=n
                 DEALLOCATE(atm (m1) % cnt)
                 ALLOCATE(atm (m1) % cnt (o1+1))
                 atm (m1) % cnt=ind_x
                 DEALLOCATE(ind_x)
              END IF
           END DO
        END DO
      END SUBROUTINE Get_Connections
  END SUBROUTINE Atoms
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

END MODULE Class_AtomCnt
