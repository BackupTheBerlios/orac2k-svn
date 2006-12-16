SUBROUTINE New_Residues

!!$***********************************************************************
!!$   Time-stamp: <2006-12-15 17:19:04 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Mon Nov 27 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*

  IMPLICIT none
  
!!$----------------------- VARIABLES IN COMMON --------------------------*

!!$------------------------- LOCAL VARIABLES ----------------------------*

  INTEGER :: n,new_r,m,i_r,i_p,p1,p2,i
  CHARACTER(len=max_char) :: pres_l,res_l,lab0
  TYPE(Unit_Char), POINTER :: Res_p,Res_n,Res_r
  TYPE(Unit_Char), DIMENSION(:), ALLOCATABLE, SAVE, TARGET :: Res_Char1
  CHARACTER(len=max_pars), DIMENSION(:), ALLOCATABLE  :: strngs
  TYPE(list), DIMENSION(:), ALLOCATABLE, SAVE :: atoms,angles,acc
  CHARACTER(len=max_data) :: errmsg_w,errmsg_f
  LOGICAL :: ok0(2)
  INTEGER :: ip, ip1, ip0
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*


  IF(.NOT. ALLOCATED(Patches)) RETURN

  m=COUNT(patches % Type == 'resi')+2*COUNT(patches % Type == 'resi')
  ALLOCATE(Add_Char(m))
  CALL copy(Res_Char,Res_Char1)

  new_r=0
  DO n=1,SIZE(patches)
     IF(patches (n) % Type == 'resi') THEN
        new_r=new_r+1
        res_l=patches(n) % res
        pres_l=patches(n) % pres
        CALL Pick(i_r,i_p)
        Res_r=>Res_Char1(i_r)
        Res_p=>Res_Char1(i_p)
        Res_n=>Add_Char(new_r)
        CALL Get_Deleted(Res_p%dele)        

        CALL Dele_group(Res_r%group)
        CALL Dele_Tpg(Res_r%bonds)
        CALL Dele_Tpg(Res_r%imph)
        CALL Dele_Tpg(Res_r%acc)
        CALL Dele_Tpg(Res_r%acc_)
        CALL Dele_Tpg(Res_r%don)
        CALL Dele_Tpg(Res_r%don_)
        CALL Dele_Tpg(Res_r%ends)

        CALL Assemble_group(Res_r%group,Res_p % group ,Res_n % group)
        CALL Assemble_Tpg(Res_r%bonds, Res_p%bonds, Res_n%bonds)
        CALL Assemble_Tpg(Res_r%imph, Res_p%imph, Res_n%imph)
        CALL Assemble_Tpg(Res_r%acc, Res_p%acc, Res_n%acc)
        CALL Assemble_Tpg(Res_r%acc_, Res_p%acc_, Res_n%acc_)
        CALL Assemble_Tpg(Res_r%don, Res_p%don, Res_n%don)
        CALL Assemble_Tpg(Res_r%don_, Res_p%don_, Res_n%don_)
        CALL Assemble_Tpg(Res_r%ends, Res_p%ends, Res_n%ends)
        Res_n%FField=Res_r%FField; Res_n%Residue='RESIDUE'; Res_n%type=patches(n)%New_Res
        CALL Verify_Ends(Res_n%bonds,Res_n%ends)
        CALL Dealloc_Deleted
     END IF
  END DO

  DO n=1,SIZE(patches)
     IF(patches (n) % Type == 'link') THEN
        ok0=.TRUE.
        DO ip0=1,2
           lab0='Link '//TRIM(patches(n)%Res_l(ip0))
           DO ip1=1,SIZE(Add_Char)
              IF(TRIM(Add_Char(ip1) % Type) == TRIM(lab0)) ok0(ip0)=.FALSE.
           END DO
        END DO
        IF(COUNT(.NOT. ok0) == 2) CYCLE
        pres_l=patches(n) % pres
        DO ip=1,2
           res_l=patches(n) % Res_l(ip)
           CALL Pick(i_r,i_p)
           WRITE(*,*) 'Ippa 2',i_r,i_p
           Res_r=>Res_Char1(i_r)
           Res_p=>Res_Char1(i_p)
           CALL Get_Deleted_link(Res_p%dele,ip)

           new_r=new_r+1
           Res_n=>Add_Char(new_r)
           CALL Dele_group(Res_r%group)
           CALL Dele_Tpg(Res_r%bonds)
           CALL Dele_Tpg(Res_r%imph)
           CALL Dele_Tpg(Res_r%acc)
           CALL Dele_Tpg(Res_r%acc_)
           CALL Dele_Tpg(Res_r%don)
           CALL Dele_Tpg(Res_r%don_)
           CALL Dele_Tpg(Res_r%ends)

           CALL Assemble_group(Res_r%group,Res_p % group ,Res_n % group)
           CALL Assemble_Tpg(Res_r%bonds, Res_p%bonds, Res_n%bonds)
           CALL Assemble_Tpg(Res_r%imph, Res_p%imph, Res_n%imph)
           CALL Assemble_Tpg(Res_r%acc, Res_p%acc, Res_n%acc)
           CALL Assemble_Tpg(Res_r%acc_, Res_p%acc_, Res_n%acc_)
           CALL Assemble_Tpg(Res_r%don, Res_p%don, Res_n%don)
           CALL Assemble_Tpg(Res_r%don_, Res_p%don_, Res_n%don_)
           CALL Assemble_Tpg(Res_r%ends, Res_p%ends, Res_n%ends)
           Res_n%FField=Res_r%FField; Res_n%Residue='RESIDUE'
           Res_n%type='Link '//TRIM(patches(n)%Res_l(ip))
           CALL Verify_Ends(Res_n%bonds,Res_n%ends)
           CALL Dealloc_Deleted
        END DO
     END IF
  END DO
  CALL Append_Tpg
  DEALLOCATE(Res_Char)
  DEALLOCATE(Res_Char1)
  DEALLOCATE(Add_Char)
  
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*
CONTAINS
  SUBROUTINE Pick(i_r,i_p)
    IMPLICIT NONE 
    INTEGER, INTENT(OUT) :: i_r,i_p
    INTEGER :: i

    i_p=-1
    i_r=-1
    
    DO i=1,SIZE(Res_Char)
       IF(MY_Fxm('RESI',Res_Char(i)%Residue)) THEN
          IF(MY_Fxm(TRIM(res_l),Res_Char(i)%Type)) THEN
             i_r=i
          END IF
       ELSE IF(MY_Fxm('PRES',Res_Char(i)%Residue)) THEN
          IF(MY_Fxm(TRIM(pres_l),Res_Char(i)%Type)) THEN
             i_p=i
          END IF
       END IF
       IF(i_p /= -1 .AND. i_r /= -1) EXIT
    END DO
    IF(i_p == -1 .OR. i_r == -1) THEN
       errmsg_f='Generating new residues '//TRIM(ADJUSTL(res_l))&
            &//' failed. Required PRES or RESIDUE not found.'
       CALL Add_error(-1,errmsg_f)
    END IF
  END SUBROUTINE Pick
  SUBROUTINE Dele_Group(group)
    TYPE(list), DIMENSION(:) :: group
    INTEGER :: n,m,o,c,p

    DO n=1,SIZE(group)
       c=0
       DO m=1,SIZE(group(n) % g)
          CALL MY_Parse(group(n) % g(m),strngs)
          DO o=1,SIZE(atoms)
             DO p=1,SIZE(atoms(o) % g) 
                IF(MY_Fxm(atoms (o) % g (p),strngs(1))) THEN
                   group(n) % g(m) = ' '
                   c=c+1
                   EXIT
                END IF
             END DO
          END DO
       END DO
       IF(c /= 0 .AND. c == SIZE(group(n) % g)) group(n) % g(1) = 'deleted'
    END DO
  END SUBROUTINE Dele_Group
  SUBROUTINE Dele_Tpg(tpg)
    CHARACTER(len=max_char), DIMENSION (:,:) :: tpg
    INTEGER :: n,m,o,p,iflag
    LOGICAL :: found
    CHARACTER(len=max_char) :: line,linea
    DO n=1,SIZE(tpg,2)
       found=.FALSE.
       DO o=1,SIZE(atoms)
          DO p=1,SIZE(atoms(o) % g) 
             linea=atoms(o) % g (p)
             DO m=1,SIZE(tpg,1)
                line=Tpg(m,n)
                IF(MY_Fxm(TRIM(linea),line)) THEN
                   found=.TRUE.
                END IF
             END DO
          END DO
       END DO
       IF(found) THEN
          Tpg(:,n)=' '
          EXIT
       END IF
    END DO

  END SUBROUTINE Dele_Tpg
  SUBROUTINE Assemble_group(group_r,group_p,group_n)
    TYPE(list), DIMENSION(:), ALLOCATABLE :: group_r,group_p,group_n
    INTEGER :: n_r,n_p,m_r,m_p,new_g,nxt_p,nxt_r,iflag,c,p1,p2
    CHARACTER(len=*), PARAMETER :: lst='   '
    CHARACTER(len=max_char) :: line,tok_r,tok_p
    LOGICAL :: ok

    DO n_r=1,SIZE(group_r)
       DO m_r=1,SIZE(group_r (n_r) % g)
          ok=.TRUE.
          line=group_r (n_r) % g(m_r)
          IF(LEN_TRIM(line) == 0) CYCLE
          nxt_r=1
          CALL Token(0,lst,line,nxt_r,tok_r,iflag)
          DO n_p=1,SIZE(group_p)
             DO m_p=1,SIZE(group_p (n_p) % g)
                line=group_p (n_p) % g(m_p)
                nxt_p=1
                CALL Token(0,lst,line,nxt_p,tok_p,iflag)
                IF(TRIM(tok_p) == TRIM(tok_r)) THEN
                   ok=.FALSE.
                END IF
             END DO
          END DO
          IF(.NOT. ok) THEN
             group_r (n_r) % g(m_r)= ' '
          END IF
       END DO
    END DO

    new_g=0
    DO n_r=1,SIZE(group_r)
       c=0
       DO m_r=1,SIZE(group_r (n_r) % g)
          IF(LEN_TRIM(group_r (n_r) % g(m_r)) /= 0) THEN
             c=c+1
          END IF
       END DO
       IF(c /= 0) THEN
          new_g=new_g+1
       END IF
    END DO
    new_g=new_g+SIZE(group_p)
    ALLOCATE(group_n(new_g))
    
    new_g=0
    DO n_r=1,SIZE(group_r)
       c=0
       DO m_r=1,SIZE(group_r (n_r) % g)
          IF(LEN_TRIM(group_r (n_r) % g(m_r)) /= 0) THEN
             c=c+1
          END IF
       END DO
       IF(c /= 0) THEN
          new_g=new_g+1
          ALLOCATE(group_n(new_g) % g(c))
          c=0
          DO m_r=1,SIZE(group_r (n_r) % g)
             IF(LEN_TRIM(group_r (n_r) % g(m_r)) /= 0) THEN
                c=c+1
                group_n(new_g) % g (c) = group_r(n_r) % g (m_r)
             END IF
          END DO
       END IF
    END DO
    DO n_p=1,SIZE(group_p)
       new_g=new_g+1
       ALLOCATE(group_n(new_g) % g(SIZE(group_p (n_p) % g)))
       DO m_p=1,SIZE(group_p (n_p) % g)
          group_n(new_g) % g = group_p(n_p) % g 
       END DO
    END DO
  END SUBROUTINE Assemble_group
  SUBROUTINE Assemble_Tpg(Tpg_r,Tpg_p,Tpg_n)
    CHARACTER(len=max_char), DIMENSION (:,:), ALLOCATABLE  :: tpg_r,Tpg_p,Tpg_n
    INTEGER :: n_r,n_p,m_r,m_p,new_g,c
    
    new_g=0
    DO n_r=1,SIZE(Tpg_r,2)
       c=0
       DO m_r=1,SIZE(Tpg_r,1)
          IF(LEN_TRIM(Tpg_r (m_r,n_r)) /= 0) THEN
             c=c+1
          END IF
       END DO
       IF(c /= 0) THEN
          new_g=new_g+1
       END IF
    END DO
    IF(ALLOCATED(Tpg_p)) THEN
       new_g=new_g+SIZE(Tpg_p,2)
    END IF
    ALLOCATE(Tpg_n(SIZE(Tpg_r,1),new_g))

    new_g=0
    DO n_r=1,SIZE(Tpg_r,2)
       c=0
       DO m_r=1,SIZE(Tpg_r,1)
          IF(LEN_TRIM(Tpg_r (m_r,n_r)) /= 0) THEN
             c=c+1
          END IF
       END DO
       IF(c /= 0) THEN
          new_g=new_g+1
          Tpg_n(:,new_g)=Tpg_r(:,n_r)
       END IF
    END DO
    IF(ALLOCATED(Tpg_p)) THEN
       DO n_p=1,SIZE(Tpg_p,2)
          new_g=new_g+1
          Tpg_n(:,new_g)=Tpg_p(:,n_p)
       END DO
    END IF
  END SUBROUTINE Assemble_Tpg
  SUBROUTINE Get_Deleted(dele)
    CHARACTER(len=max_char), DIMENSION (:,:) :: dele
    INTEGER :: o,nword,n_at,n_ang,n_acc

    n_at=0
    n_ang=0
    n_acc=0
    DO o=1,SIZE(dele,2)
       CALL My_Parse(dele(1,o),strngs)
       IF(MY_Fxm('atom',strngs(1))) THEN
          n_at=n_at+1
       ELSE IF(MY_Fxm('angle',strngs(1))) THEN
          n_ang=n_ang+1
       ELSE IF(MY_Fxm('acc',strngs(1))) THEN
          n_acc=n_acc+1
       END IF
    END DO
    ALLOCATE(atoms(n_at),angles(n_ang),acc(n_acc))
    n_at=0
    n_ang=0
    n_acc=0
    DO o=1,SIZE(dele,2)
       CALL My_Parse(dele(1,o),strngs)
       IF(MY_Fxm('atom',strngs(1))) THEN
          n_at=n_at+1
          nword=SIZE(strngs)
          ALLOCATE(atoms(n_at)%g(nword-1))
          atoms(n_at) % g = strngs(2:)
       ELSE IF(MY_Fxm('angle',strngs(1))) THEN
          n_ang=n_ang+1
          nword=SIZE(strngs)
          ALLOCATE(angles(n_at)%g(nword-1))
          angles(n_ang) % g = strngs(2:)          
       ELSE IF(MY_Fxm('acc',strngs(1))) THEN
          n_acc=n_acc+1
          nword=SIZE(strngs)
          ALLOCATE(acc(n_at)%g(nword-1))
          acc(n_acc) % g = strngs(2:)
       END IF
    END DO
  END SUBROUTINE Get_Deleted
  SUBROUTINE Get_Deleted_link(dele,which)
    CHARACTER(len=max_char), DIMENSION (:,:) :: dele
    INTEGER :: which
    INTEGER :: o,nword,n_at,n_ang,n_acc,n_bend,p
    CHARACTER(len=max_char) :: C_Which,lab0
    n_at=0
    n_ang=0
    n_acc=0
    WRITE(c_which,'(i1)') which
    DO o=1,SIZE(dele,2)
       CALL My_Parse(dele(1,o),strngs)
       lab0=TRIM(strngs(2))
       IF(MY_Fxm('atom',strngs(1)) .AND. lab0(1:1) == C_Which(1:1)) THEN
          n_at=n_at+1
       ELSE IF(MY_Fxm('angle',strngs(1))) THEN
          nword=SIZE(strngs)
          n_bend=0
          DO p=2,nword
             lab0=TRIM(strngs(p))
             IF(lab0(1:1) == C_Which(1:1)) THEN
                n_bend=n_bend+1
             END IF
          END DO
          IF(n_Bend /= 0) n_ang=n_ang+1
       ELSE IF(MY_Fxm('acc',strngs(1))) THEN
          n_acc=n_acc+1
       END IF
    END DO
    ALLOCATE(atoms(n_at),angles(n_ang),acc(n_acc))
    n_at=0
    n_ang=0
    n_acc=0
    DO o=1,SIZE(dele,2)
       CALL My_Parse(dele(1,o),strngs)
       lab0=TRIM(strngs(2))
       IF(MY_Fxm('atom',strngs(1)) .AND. lab0(1:1) == C_Which(1:1)) THEN
          n_at=n_at+1
          nword=SIZE(strngs)
          ALLOCATE(atoms(n_at)%g(nword-1))
          DO p=2,nword
             atoms(n_at) % g (p-1)  = strngs(p)(2:)
          END DO
       ELSE IF(MY_Fxm('angle',strngs(1))) THEN
          nword=SIZE(strngs)
          n_bend=0
          DO p=2,nword
             lab0=TRIM(strngs(p))
             IF(lab0(1:1) == C_Which(1:1)) THEN
                n_bend=n_bend+1
             END IF
          END DO
          IF(n_Bend /= 0) n_ang=n_ang+1
          ALLOCATE(angles(n_ang)%g(n_bend))
          n_bend=0
          DO p=2,nword
             lab0=TRIM(strngs(p))
             IF(lab0(1:1) == C_Which(1:1)) THEN
                n_bend=n_bend+1
                angles(n_ang) % g (n_bend)=lab0(2:) 
             END IF
          END DO
          
       ELSE IF(MY_Fxm('acc',strngs(1))) THEN
          n_acc=n_acc+1
          nword=SIZE(strngs)
          ALLOCATE(acc(n_at)%g(nword-1))
          acc(n_acc) % g = strngs(2:)
       END IF
    END DO
  END SUBROUTINE Get_Deleted_link
  SUBROUTINE Verify_Ends(bonds,ends)
    CHARACTER(len=max_char), DIMENSION (:,:)  :: bonds,ends
    INTEGER, PARAMETER :: ntypes=2
    CHARACTER(len=7), DIMENSION(ntypes), SAVE :: atom_type=(/'n','c'/)
    INTEGER, DIMENSION(ntypes), SAVE :: atom_max=(/2,2/)
    INTEGER :: ind(2),n,m,p,n_bnd
    CHARACTER(len=max_char) :: line
    INTEGER, SAVE :: first_time=0
    
    ind=-1
    DO n=1,2
       line=TRIM(ADJUSTL(ends(n,1)))
       n_bnd=COUNT(ends(n,1) == bonds(:,:))
       DO p=1,ntypes
          IF(atom_type(p) == line .AND. n_bnd > atom_max(p)) THEN
             ind(n)=n
             EXIT
          END IF
       END DO
    END DO

    DO n=1,SIZE(ends,1)
       IF(ind(n) /= -1) THEN
          ends(n,1)=' * '
       END IF
    END DO

  END SUBROUTINE Verify_Ends
  SUBROUTINE Dealloc_deleted
    INTEGER :: o
    IF(ALLOCATED(atoms)) THEN
       DO o=1,SIZE(atoms)
          DEALLOCATE(atoms(o) % g)
       END DO
       DEALLOCATE(atoms)
    END IF
    IF(ALLOCATED(angles)) THEN
       DO o=1,SIZE(angles)
          DEALLOCATE(angles(o) % g)
       END DO
       DEALLOCATE(angles)
    END IF
    IF(ALLOCATED(acc)) THEN
       DO o=1,SIZE(acc)
          DEALLOCATE(acc(o) % g)
       END DO
       DEALLOCATE(acc)
    END IF
  END SUBROUTINE Dealloc_deleted
  SUBROUTINE Append_Tpg
    INTEGER :: o,n,m

    ALLOCATE(App_Char(SIZE(Res_Char)+SIZE(Add_Char)))
    DO n=1,SIZE(Res_Char)
       CALL Copy(Res_Char, App_Char, n, n)
    END DO
    o=SIZE(Res_Char)
    DO n=1,SIZE(Add_Char)
       o=o+1
       CALL Copy(Add_Char, App_Char, n, o)
    END DO
  END SUBROUTINE Append_Tpg

END SUBROUTINE New_Residues
