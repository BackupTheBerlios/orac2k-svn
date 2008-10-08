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
SUBROUTINE NewResidues_
  IMPLICIT none
  INTEGER :: n,new_r,m,i_r,i_p,p1,p2,i
  CHARACTER(len=max_char) :: pres_l,res_l,lab0
  TYPE(Tops__Type), POINTER :: Res_r,Res_p,Res_n
  TYPE(list), DIMENSION(:), ALLOCATABLE :: atoms,angles,acc
  LOGICAL :: ok0(2),ok_link,ok_resi
  INTEGER :: ip, ip1, ip0,ipp,ipp_end,o,p
  CHARACTER(len=max_char), DIMENSION(:,:), ALLOCATABLE :: share
  TYPE(list), DIMENSION(:), ALLOCATABLE :: shareg
  TYPE(Tops__Type), DIMENSION(:), POINTER :: Stores=>NULL()

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*


  IF(.NOT. ALLOCATED(Patches)) THEN
     CALL Append_Tpg
     DEALLOCATE(Res_Char)
     RETURN
  END IF

  stores=>Tops__Store(Res_Char)
  m=COUNT(patches % Type == 'resi')+2*COUNT(patches % Type == 'link')
  ALLOCATE(Add_Char(m))
  new_r=0
  DO n=1,SIZE(patches)
     ok_resi=.FALSE.; ok_Link=.FALSE.
     ipp_end=1
     IF(patches (n) % Type == 'resi') ok_resi=.TRUE.
     IF(patches (n) % Type == 'link') THEN
        ok0=.TRUE.
        DO ip0=1,2
           lab0='Link '//TRIM(patches(n)%Res_l(ip0))
           DO ip1=1,SIZE(Add_Char)
              IF(TRIM(Add_Char(ip1) % Type) == TRIM(lab0)) ok0(ip0)=.FALSE.
           END DO
        END DO
        IF(COUNT(.NOT. ok0) == 2) CYCLE
        ok_link=.TRUE.
        ipp_end=2
     END IF
     pres_l=patches(n) % pres
     DO ipp=1,ipp_end
        new_r=new_r+1
        IF(ok_resi) THEN
           res_l=patches(n) % res
        ELSE
           res_l=patches(n) % Res_l(ipp)
        END IF

        CALL Pick(i_r,i_p)
        

        Res_r=>Res_Char(i_r)
        Res_p=>Res_Char(i_p)
        Res_n=>Add_Char(new_r)

        IF(ok_resi) CALL Get_Deleted(Res_p%dele)

        IF(ok_link) CALL Get_Deleted_link(Res_p%dele,ipp)

        CALL Dele_group(Res_r%group)
        CALL Dele_Tpg(Res_r%bonds)
        CALL Dele_Tpg(Res_r%imph)
        CALL Dele_Tpg(Res_r%acc)
        
        CALL Dele_Tpg(Res_r%acc_)
        CALL Dele_Tpg(Res_r%don)
        CALL Dele_Tpg(Res_r%don_)
        CALL Dele_Tpg(Res_r%ends)

        IF(ok_resi) CALL Assemble_group(Res_r%group,Res_p % group)
        IF(ok_link) CALL Assemble_group(Res_r%group,Res_p % group, ipp)
        IF(ALLOCATED(shareg)) THEN
           IF(ALLOCATED(Add_Char(new_r) % group)) DEALLOCATE(Add_Char(new_r) % group)
           ALLOCATE(ADD_Char(new_r) % group (SIZE(shareg)))
           DO m=1,SIZE(shareg)
              IF(ALLOCATED(shareg(m) % g)) THEN
                 ALLOCATE(Add_Char(new_r) % group(m) % g(SIZE(shareg(m) % g)))
                 Add_Char(new_r) % group(m) % g=shareg(m) % g
              END IF
           END DO
        END IF

        CALL Assemble_Tpg(Res_r%bonds, Res_p%bonds)
        IF(ALLOCATED(share)) THEN
           IF(ALLOCATED(Add_Char(new_r) %bonds)) DEALLOCATE(Add_Char(new_r) %bonds)
           ALLOCATE(Add_Char(new_r) % bonds(SIZE(share,1),SIZE(share,2)))
           Add_Char(new_r) % bonds = share
        END IF
        CALL Assemble_Tpg(Res_r%imph, Res_p%imph)
        IF(ALLOCATED(share)) THEN
           IF(ALLOCATED(Add_Char(new_r) %imph)) DEALLOCATE(Add_Char(new_r) %imph)
           ALLOCATE(Add_Char(new_r) % imph(SIZE(share,1),SIZE(share,2)))
           Add_Char(new_r) % imph = share
        END IF
        CALL Assemble_Tpg(Res_r%acc, Res_p%acc)
        IF(ALLOCATED(share)) THEN
           IF(ALLOCATED(Add_Char(new_r) %acc)) DEALLOCATE(Add_Char(new_r) %acc)
           ALLOCATE(Add_Char(new_r) % acc(SIZE(share,1),SIZE(share,2)))
           Add_Char(new_r) % acc = share
        END IF
        CALL Assemble_Tpg(Res_r%acc_, Res_p%acc_)
        IF(ALLOCATED(share)) THEN
           IF(ALLOCATED(Add_Char(new_r) %acc_)) DEALLOCATE(Add_Char(new_r) %acc_)
           ALLOCATE(Add_Char(new_r) % acc_(SIZE(share,1),SIZE(share,2)))
           Add_Char(new_r) % acc_ = share
        END IF
        CALL Assemble_Tpg(Res_r%don, Res_p%don)
        IF(ALLOCATED(share)) THEN
           IF(ALLOCATED(Add_Char(new_r) %don)) DEALLOCATE(Add_Char(new_r) %don)
           ALLOCATE(Add_Char(new_r) % don(SIZE(share,1),SIZE(share,2)))
           Add_Char(new_r) % don = share
        END IF
        CALL Assemble_Tpg(Res_r%don_, Res_p%don_)
        IF(ALLOCATED(share)) THEN
           IF(ALLOCATED(Add_Char(new_r) % don_)) DEALLOCATE(Add_Char(new_r) % don_)
           ALLOCATE(Add_Char(new_r) %  don_(SIZE(share,1),SIZE(share,2)))
           Add_Char(new_r) %  don_ = share
        END IF
        CALL Assemble_Tpg(Res_r%ends, Res_p%ends)
        IF(ALLOCATED(share)) THEN
           IF(ALLOCATED(Add_Char(new_r) % ends)) DEALLOCATE(Add_Char(new_r) % ends)
           ALLOCATE(Add_Char(new_r) %  ends(SIZE(share,1),SIZE(share,2)))
           Add_Char(new_r) %  ends = share
        END IF
        Res_n%FField=Res_r%FField; Res_n%Residue='RESIDUE'
        IF(ok_resi) Res_n%type=patches(n)%New_Res
        IF(ok_link) Res_n%type='Link '//TRIM(patches(n)%Res_l(ipp))
        CALL Verify_Ends(Res_n%bonds,Res_n%ends)
        CALL Dealloc_Deleted
        Res_Char(i_r)=Stores(i_r)
     END DO
  END DO

  CALL Append_Tpg
  DEALLOCATE(Res_Char)
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
          IF(MY_Fxm(TRIM(res_l),TRIM(Res_Char(i)%Type))) THEN
             i_r=i
          END IF
       ELSE IF(MY_Fxm('PRES',Res_Char(i)%Residue)) THEN
          IF(MY_Fxm(TRIM(pres_l),TRIM(Res_Char(i)%Type))) THEN
             i_p=i
          END IF
       END IF
       IF(i_p /= -1 .AND. i_r /= -1) EXIT
    END DO

    IF(i_p == -1 .OR. i_r == -1) THEN
       errmsg_f='Generating new residues '//TRIM(ADJUSTL(res_l))&
            &//' failed. Required PRES or RESIDUE not found.'
       CALL Add_errors(-1,errmsg_f)
       CALL Print_Errors()
    END IF
  END SUBROUTINE Pick
  SUBROUTINE Dele_Group(group)
    TYPE(list), DIMENSION(:), ALLOCATABLE :: group
    INTEGER :: n,m,o,c,p,nword
    IF(.NOT. ALLOCATED(group)) RETURN

    DO n=1,SIZE(group)
       c=0
       DO m=1,SIZE(group(n) % g)
          nword=MYParse_(group(n) % g(m))
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
    CHARACTER(len=max_char), DIMENSION (:,:), ALLOCATABLE :: tpg
    INTEGER :: n,m,o,p,iflag,nword
    LOGICAL :: found
    CHARACTER(len=max_char) :: line,linea

    IF(.NOT. ALLOCATED(tpg)) RETURN
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
  SUBROUTINE Assemble_group(group_r,group_pa,which)
    TYPE(list), DIMENSION(:), ALLOCATABLE :: group_r,group_pa
    INTEGER, OPTIONAL :: which
    INTEGER :: n_r,n_p,m_r,m_p,new_g,nxt_p,nxt_r,iflag,c,p1,p2,p_dim,nword,nn_p,mm_p
    CHARACTER(len=*), PARAMETER :: lst='   '
    CHARACTER(len=max_char) :: line,tok_r,tok_p
    LOGICAL :: ok
    CHARACTER(len=max_char) :: C_Which,lab0
    TYPE(list), DIMENSION(:), ALLOCATABLE :: group_t
    TYPE(list), DIMENSION(:), ALLOCATABLE :: group_p

    IF(ALLOCATED(shareg)) DEALLOCATE(shareg)
    IF(ALLOCATED(group_p)) DEALLOCATE(group_p)

    IF(ALLOCATED(group_pa)) THEN
       ALLOCATE(group_p(SIZE(group_pa)))
       DO n_r=1,SIZE(group_pa)
          IF(ALLOCATED(group_pa(n_r) % g)) THEN
             ALLOCATE(group_p(n_r) % g(SIZE(group_pa(n_r) % g)))
             group_p(n_r) % g=group_pa(n_r) % g
          END IF
       END DO
    END IF

    IF(present(which)) THEN
       WRITE(c_which,'(i1)') which
       ALLOCATE(group_t(SIZE(group_p)))
       DO n_p=1,SIZE(group_p)
          ALLOCATE(group_t(n_p) % g (SIZE(group_p (n_p) % g)))
       END DO
       DO n_p=1,SIZE(group_p)
          p_dim=0
          DO m_p=1,SIZE(group_p (n_p) % g)
             lab0=ADJUSTL(TRIM(group_p (n_p) % g(m_p)))
             IF(lab0(1:1) == C_which(1:1)) THEN
                group_t (n_p) % g(m_p)=lab0(2:)
                p_dim=p_dim+1
             ELSE
                group_t (n_p) % g(m_p)=' '
             END IF
          END DO
          DEALLOCATE(group_p (n_p) % g)
          ALLOCATE(group_p(n_p) % g (p_dim))
          p_dim=0
          DO m_p=1,SIZE(group_t (n_p) % g)
             lab0=TRIM(group_t (n_p) % g(m_p))
             IF(LEN_TRIM(lab0) /= 0 ) THEN
                p_dim=p_dim+1
                group_p (n_p) % g(p_dim)=lab0
             END IF
          END DO
       END DO
       DEALLOCATE(group_t)
    END IF
    
    DO n_r=1,SIZE(group_r)
       DO m_r=1,SIZE(group_r (n_r) % g)
          ok=.TRUE.
          line=group_r (n_r) % g(m_r)
          IF(LEN_TRIM(line) == 0) CYCLE
          nword=MyParse_(line)
          tok_r=strngs(1)
          DO n_p=1,SIZE(group_p)
             DO m_p=1,SIZE(group_p (n_p) % g)
                line=group_p (n_p) % g(m_p)
                nword=MyParse_(line)
                tok_p=strngs(1)
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
    ALLOCATE(shareg(new_g))
    
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
          ALLOCATE(shareg(new_g) % g(c))
          c=0
          DO m_r=1,SIZE(group_r (n_r) % g)
             IF(LEN_TRIM(group_r (n_r) % g(m_r)) /= 0) THEN
                c=c+1
                shareg(new_g) % g (c) = group_r(n_r) % g (m_r)
             END IF
          END DO
       END IF
    END DO
    DO n_p=1,SIZE(group_p)
       new_g=new_g+1
       ALLOCATE(shareg(new_g) % g(SIZE(group_p (n_p) % g)))
       DO m_p=1,SIZE(group_p (n_p) % g)
          shareg(new_g) % g = group_p(n_p) % g 
       END DO
    END DO
  END SUBROUTINE Assemble_group
  SUBROUTINE Assemble_Tpg(Tpg_r,Tpg_p)
    CHARACTER(len=max_char), DIMENSION (:,:), ALLOCATABLE  :: tpg_r,Tpg_p
    INTEGER :: n_r,n_p,m_r,m_p,new_g,c,nn,p_dim

    CHARACTER(len=max_char), DIMENSION (:,:), ALLOCATABLE  :: tpg_t
    CHARACTER(len=max_char), DIMENSION (:), ALLOCATABLE  :: lab0
    CHARACTER(len=max_char) :: C_Which

    IF(ALLOCATED(share)) DEALLOCATE(share)
    new_g=0
    IF(ALLOCATED(Tpg_r)) THEN
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
    END IF
    IF(ALLOCATED(Tpg_p)) THEN
       new_g=new_g+SIZE(Tpg_p,2)
    END IF
    IF(new_g == 0) RETURN

    ALLOCATE(share(SIZE(Tpg_r,1),new_g))

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
          share(:,new_g)=Tpg_r(:,n_r)
       END IF
    END DO
    IF(ALLOCATED(Tpg_p)) THEN
       DO n_p=1,SIZE(Tpg_p,2)
          new_g=new_g+1
          share(:,new_g)=Tpg_p(:,n_p)
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
       nword=MyParse_(dele(1,o))
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
       nword=MyParse_(dele(1,o))
       IF(MY_Fxm('atom',strngs(1))) THEN
          n_at=n_at+1
          ALLOCATE(atoms(n_at)%g(nword-1))
          atoms(n_at) % g = strngs(2:)
       ELSE IF(MY_Fxm('angle',strngs(1))) THEN
          n_ang=n_ang+1
          ALLOCATE(angles(n_at)%g(nword-1))
          angles(n_ang) % g = strngs(2:)          
       ELSE IF(MY_Fxm('acc',strngs(1))) THEN
          n_acc=n_acc+1
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
       nword=MyParse_(dele(1,o))
       lab0=TRIM(strngs(2))
       IF(MY_Fxm('atom',strngs(1)) .AND. lab0(1:1) == C_Which(1:1)) THEN
          n_at=n_at+1
       ELSE IF(MY_Fxm('angle',strngs(1))) THEN
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
       nword=MyParse_(dele(1,o))
       lab0=TRIM(strngs(2))
       IF(MY_Fxm('atom',strngs(1)) .AND. lab0(1:1) == C_Which(1:1)) THEN
          n_at=n_at+1
          ALLOCATE(atoms(n_at)%g(nword-1))
          DO p=2,nword
             atoms(n_at) % g (p-1)  = strngs(p)(2:)
          END DO
       ELSE IF(MY_Fxm('angle',strngs(1))) THEN
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
    INTEGER :: n,m,o(2),p,q
    TYPE(Tops__Type), DIMENSION(:), POINTER :: Tops__

    
    q=0
    IF(ALLOCATED(Add_Char)) THEN
       DO n=1,SIZE(Add_Char)
          IF(ALLOCATED(Add_Char (n) % group)) q=q+1
       END DO
    END IF
    ALLOCATE(App_Char(SIZE(Res_Char)+q))

    Tops__=>Res_Char
    q=0
    DO n=1,SIZE(App_Char)
       p=n
       IF(n > SIZE(Res_Char)) THEN
          p=n-SIZE(Res_Char)
          Tops__=>Add_Char
          IF(.NOT. ALLOCATED(Add_Char (p) % group)) CYCLE
          q=q+1
       ELSE
          q=q+1
       END IF
       App_Char(q) % FField=Tops__(p) % FField
       App_Char(q) % Type=Tops__(p) % Type
       App_Char(q) % Residue=Tops__(p) % Residue
       IF(ALLOCATED(Tops__(p) % bonds)) THEN
          o=SHAPE(Tops__(p) % bonds)
          ALLOCATE(App_Char(q) % bonds(o(1), o(2)))
          App_Char(q) % Bonds=Tops__(p) % bonds
       END IF
       IF(ALLOCATED(Tops__(p) % imph)) THEN
          o=SHAPE(Tops__(p) % imph)
          ALLOCATE(App_Char(q) % imph(o(1), o(2)))
          App_Char(q) % imph=Tops__(p) % imph
       END IF
       IF(ALLOCATED(Tops__(p) % acc)) THEN
          o=SHAPE(Tops__(p) % acc)
          ALLOCATE(App_Char(q) % acc(o(1), o(2)))
          App_Char(q) % acc=Tops__(p) % acc
       END IF
       IF(ALLOCATED(Tops__(p) % don)) THEN
          o=SHAPE(Tops__(p) % don)
          ALLOCATE(App_Char(q) % don(o(1), o(2)))
          App_Char(q) % don=Tops__(p) % don
       END IF
       IF(ALLOCATED(Tops__(p) % acc_)) THEN
          o=SHAPE(Tops__(p) % acc_)
          ALLOCATE(App_Char(q) % acc_(o(1), o(2)))
          App_Char(q) % acc_=Tops__(p) % acc_
       END IF
       IF(ALLOCATED(Tops__(p) % don_)) THEN
          o=SHAPE(Tops__(p) % don_)
          ALLOCATE(App_Char(q) % don_(o(1), o(2)))
          App_Char(q) % don_=Tops__(p) % don_
       END IF
       IF(ALLOCATED(Tops__(p) % dele)) THEN
          o=SHAPE(Tops__(p) % dele)
          ALLOCATE(App_Char(q) % dele(o(1), o(2)))
          App_Char(q) % dele=Tops__(p) % dele
       END IF
       IF(ALLOCATED(Tops__(p) % ends)) THEN
          o=SHAPE(Tops__(p) % ends)
          ALLOCATE(App_Char(q) % ends(o(1), o(2)))
          App_Char(q) % ends=Tops__(p) % ends
       END IF
       IF(ALLOCATED(Tops__(p) % mass)) THEN
          o=SHAPE(Tops__(p) % mass)
          ALLOCATE(App_Char(q) % mass(o(1), o(2)))
          App_Char(q) % mass=Tops__(p) % mass
       END IF
       IF(ALLOCATED(Tops__(p) % group)) THEN
          ALLOCATE(App_Char(q) % group (SIZE(Tops__(p) % group)))
          DO m=1,SIZE(Tops__(p) % group)
             IF(ALLOCATED(Tops__(p) % group (m) % g)) THEN
                ALLOCATE(App_Char(q) % group (m) % g (SIZE(Tops__(p) % group (m) % g)))
                App_Char(q) % group(m) % g =Tops__(p) % group(m) % g
             END IF
          END DO
       END IF
    END DO
  END SUBROUTINE Append_Tpg
  FUNCTION Tops__Store(Tops__, delete) RESULT(out)
    INTEGER, OPTIONAL :: delete
    TYPE(Tops__type), DIMENSION(:), ALLOCATABLE :: Tops__
    TYPE(Tops__type), DIMENSION(:), ALLOCATABLE, SAVE, TARGET :: store
    TYPE(Tops__type), DIMENSION(:), POINTER  :: out
    INTEGER :: o(2),n,m

    out=>NULL()
    IF(PRESENT(delete)) THEN
       IF(ALLOCATED(Store)) DEALLOCATE(Store)
       RETURN
    END IF
    IF(ALLOCATED(Store)) DEALLOCATE(Store)
    IF(ALLOCATED(Tops__)) THEN
       ALLOCATE(Store(SIZE(Tops__)))
       DO n=1,SIZE(Tops__)
          o=0
          IF(ALLOCATED(Tops__(n) % bonds)) THEN
             o=SHAPE(Tops__(n) % bonds)
             ALLOCATE(Store(n) % bonds(o(1), o(2)))
             Store(n) % bonds=Tops__(n) % bonds
          END IF
          IF(ALLOCATED(Tops__(n) % imph)) THEN
             o=SHAPE(Tops__(n) % imph)
             ALLOCATE(Store(n) % imph(o(1), o(2)))
             Store(n) % imph=Tops__(n) %  imph
          END IF
          IF(ALLOCATED(Tops__(n) % acc)) THEN
             o=SHAPE(Tops__(n) % acc)
             ALLOCATE(Store(n) % acc(o(1), o(2)))
             Store(n) % acc=Tops__(n) %  acc
          END IF
          IF(ALLOCATED(Tops__(n) % don)) THEN
             o=SHAPE(Tops__(n) % don)
             ALLOCATE(Store(n) % don(o(1), o(2)))
             Store(n) % don=Tops__(n) %  don
          END IF
          IF(ALLOCATED(Tops__(n) % acc_)) THEN
             o=SHAPE(Tops__(n) % acc_)
             ALLOCATE(Store(n) % acc_(o(1), o(2)))
             Store(n) % acc_=Tops__(n) %  acc_
          END IF
          IF(ALLOCATED(Tops__(n) % don_)) THEN
             o=SHAPE(Tops__(n) % don_)
             ALLOCATE(Store(n) % don_(o(1), o(2)))
             Store(n) % don_=Tops__(n) %  don_
          END IF
          IF(ALLOCATED(Tops__(n) % dele)) THEN
             o=SHAPE(Tops__(n) % dele)
             ALLOCATE(Store(n) % dele(o(1), o(2)))
             Store(n) % dele=Tops__(n) %  dele
          END IF
          IF(ALLOCATED(Tops__(n) % ends)) THEN
             o=SHAPE(Tops__(n) % ends)
             ALLOCATE(Store(n) % ends(o(1), o(2)))
             Store(n) % ends=Tops__(n) %  ends
          END IF
          IF(ALLOCATED(Tops__(n) % mass)) THEN
             o=SHAPE(Tops__(n) % mass)
             ALLOCATE(Store(n) % mass(o(1), o(2)))
             Store(n) % mass=Tops__(n) %  mass
          END IF
          IF(ALLOCATED(Tops__(n) % group)) THEN
             ALLOCATE(Store(n) % group (SIZE(Tops__(n) % group)))
             DO m=1,SIZE(Tops__(n) % group)
                IF(ALLOCATED(Tops__(n) % group (m) % g)) THEN
                   ALLOCATE(Store(n) % group (m) % g (SIZE(Tops__(n) % group (m) % g)))
                   Store(n) % group(m) % g =Tops__(n) % group(m) % g
                END IF
             END DO
          END IF
          Store(n) % FField=Tops__(n) % FField
          Store(n) % Type=Tops__(n) % Type
          Store(n) % Residue=Tops__(n) % Residue
       END DO
       out=>Store
    END IF
  END FUNCTION Tops__Store

END SUBROUTINE NewResidues_
