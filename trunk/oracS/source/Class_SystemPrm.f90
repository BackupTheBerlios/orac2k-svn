MODULE Class_SystemPrm

!!$***********************************************************************
!!$   Time-stamp: <2006-12-18 22:20:35 marchi>                           *
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
  USE ERROR_Mod, ONLY: Add_error=>Add,Print_Errors
  USE MYPARSE_Mod, ONLY: MY_Parse=>parse
  USE STRINGS_Mod, ONLY: MY_FxM
  USE CONSTANTS
  USE UNITS
  USE STRPAK, ONLY: SP_getnum
  USE Linked_Real8_D
  USE Class_SystemTpg, ONLY: SystemTpg__Init=>Init, SystemTpg, Tpg, Atom2Cnt
  IMPLICIT none
  PRIVATE
  PUBLIC :: Init, SystemPrm, Prm
  INTEGER, PARAMETER :: max_atms=7

  TYPE SystemPrm
     TYPE(Chain_M), DIMENSION(:), ALLOCATABLE :: bonds,angles,&
          &dihed,imph,int14,LJs
     REAL(8), DIMENSION (:), ALLOCATABLE :: c12,c6,c12_,c6_
     CHARACTER(len=max_atms), DIMENSION(:), ALLOCATABLE :: Types
  END TYPE SystemPrm
  TYPE Par_Pot
     INTEGER, DIMENSION(:), ALLOCATABLE :: b
     REAL(8), DIMENSION(:), ALLOCATABLE :: p
  END TYPE Par_Pot

  TYPE(SystemPrm), SAVE :: Prm
  CHARACTER(len=max_pars), DIMENSION(:), ALLOCATABLE :: strngs
  CHARACTER(len=max_char_long) :: errmsg_f
  TYPE(Par_Pot), DIMENSION(:), ALLOCATABLE, SAVE ::  p_Tpg
  REAL(8), DIMENSION(:), ALLOCATABLE  :: Params
  LOGICAL, DIMENSION(:), ALLOCATABLE  :: oks
CONTAINS
  SUBROUTINE Init

    CALL SystemTpg__Init    
    
    CALL Types
    CALL Lennard_Jones
    CALL Bonds
    CALL Angles
    CALL Dihedrals
    CALL Improper
    CALL Print_Errors()
    
  CONTAINS
    SUBROUTINE Types
      INTEGER :: n,ntypes,i_bonds,count_out,count_A,iflags,nlines,nword
      CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE :: rewrite
      CHARACTER(len=max_char) :: lab0

      i_Bonds=-1
      DO n=1,SIZE(Topology)
         IF(My_Fxm('HEADER', Topology (n) % Residue )) THEN
            i_bonds=n; EXIT
         END IF
      END DO
      IF(i_Bonds == -1) THEN
         errmsg_f='Can''t find atom type labels in the input file'
         CALL Add_Error(-1,errmsg_f)
         RETURN
      END IF
      ntypes=SIZE(Topology(i_bonds) % line)
      count_a=0
      DO n=1,ntypes
         IF(My_FxM('mass',Topology (i_Bonds) % line(n))) THEN
            count_a=count_a+1
         END IF
      END DO
      ALLOCATE(Prm % Types(count_a))
      count_a=0
      DO n=1,ntypes
         IF(My_FxM('mass',Topology (i_Bonds) % line(n))) THEN
            count_a=count_a+1
            CALL My_Parse(Topology (i_Bonds) % line(n), strngs)
            nword=SIZE(strngs)
            IF(nword < 3) THEN
               errmsg_f='Wrong format on the atom type labels file.&
                    & At least four fields are expected. It was found: '//&
                    &TRIM(Topology (i_Bonds) % line(n))
               CALL Add_Error(-1,errmsg_f)
               CALL Print_Errors()
            END IF
            Prm % Types(count_a) = TRIM(strngs(3))
         END IF
      END DO

      ALLOCATE(rewrite(count_a),oks(count_a))
      oks=.FALSE.
      DO n=1,SIZE(Tpg % atm)
         lab0=TRIM(Tpg % atm(n) % a % betab)
         oks(Type_Number(lab0))=.TRUE.
      END DO
      count_a=0
      DO n=1,SIZE(Prm % Types)
         IF(oks(n)) THEN
            count_a=count_a+1
            rewrite(count_a)=Prm % Types(n)
         END IF
      END DO
      DEALLOCATE(Prm % Types)
      ALLOCATE(Prm % Types (count_a))
      Prm % Types=rewrite(1:count_a)
      DEALLOCATE(oks,rewrite)
      
      WRITE(*,*) 'Atomic types in use No. =====>',SIZE(Prm % Types)
    END SUBROUTINE Types
    SUBROUTINE Lennard_Jones
      INTEGER :: n,m,ij,ntypes,i,i_bonds,count_out,count_A,iflags&
           &,nlines,nword,m1,ip
      LOGICAL :: ok,noerror
      INTEGER, DIMENSION(:), ALLOCATABLE :: ind_x
      REAL(8) :: eps,sigma,eps14,sigma14,tmpi
      TYPE(ChainR8C), DIMENSION(:), ALLOCATABLE :: LJ
      CHARACTER(len=max_char) :: lab0

      i_Bonds=-1
      DO n=1,SIZE(Paras)
         IF(My_Fxm('NONBOND', Paras (n) % Residue )) THEN
            i_bonds=n; EXIT
         END IF
      END DO
      IF(i_Bonds == -1) THEN
         errmsg_f='Can''t find atom type labels in the input file'
         CALL Add_Error(-1,errmsg_f)
         RETURN
      END IF

      ALLOCATE(ind_x(SIZE(Paras(i_Bonds) % line)))

      ALLOCATE(Prm % LJs (SIZE(Prm % Types)))
      ALLOCATE(LJ(SIZE(Paras(i_Bonds) % line)))

      count_a=0
      DO n=1,SIZE(Paras(i_Bonds) % line)
         CALL MY_Parse(Paras(i_Bonds) % line (n), strngs)
         nword=SIZE(strngs)
         ok=.TRUE.
         IF(nword == 1) ok=.FALSE.
         IF(ALLOCATED(Params)) DEALLOCATE(Params)
         ALLOCATE(Params(nword-1))
         DO i=2,nword
            CALL SP_Getnum(strngs(i),Params(i-1),iflags)
            IF(iflags /= 0) THEN
               ok=.FALSE.
               EXIT
            END IF
         END DO
         IF(ok) THEN
            count_a=count_a+1
            LJ(count_a) % lab =TRIM(strngs(1))
            ALLOCATE(LJ(count_a) % g(nword-1))
            LJ(count_a) % g = Params
            DEALLOCATE(Params)
         END IF
      END DO
      ALLOCATE(oks(SIZE(Prm % LJs)))
      
      oks=.FALSE.
      DO n=1,count_a
         ip=Type_Number(LJ(n) % lab, noerror)
         
         IF(ip > 0) THEN
            oks(ip)=.TRUE.
            Prm % LJs (ip) % pt = ip

            m1=SIZE(LJ(n) % g)
            SELECT CASE(m1)
            CASE(3)
               sigma=LJ(n) % g (3)*LJ_fact
               IF(LJ(n) % g (1) == 0.0D0) THEN
                  eps=-LJ(n) % g (2)
               ELSE IF(LJ(n) % g (2) == 0.0D0) THEN
                  eps=LJ(n) % g (1)
               ELSE
                  errmsg_f='Unrecognized Lennard-Jones format: '&
                       &//TRIM(Paras(i_Bonds) % line(ind_x(count_a)))
                  CALL Add_Error(-1,errmsg_f)
               END IF
               ALLOCATE(Prm % LJs (ip) % g (2))
               Prm % LJs (ip) % g =(/sigma,eps/)
            CASE(6)
               sigma=LJ(n) % g (3)*LJ_fact
               IF(LJ(n) % g (1) == 0.0D0) THEN
                  eps=-LJ(n) % g (2)
               ELSE IF(LJ(n) % g (2) == 0.0D0) THEN
                  eps=LJ(n) % g (1)
               ELSE
                  errmsg_f='Unrecognized Lennard-Jones format: '&
                       &//TRIM(Paras(i_Bonds) % line(ind_x(count_a)))
                  CALL Add_Error(-1,errmsg_f)
               END IF
               sigma14=LJ(n) % g (6)*LJ_fact
               IF(LJ(n) % g (4) == 0.0D0) THEN
                  eps14=-LJ(n) % g (5)
               ELSE IF(LJ(n) % g (5) == 0.0D0) THEN
                  eps14=LJ(n) % g (4)
               ELSE
                  errmsg_f='Unrecognized Lennard-Jones format: '&
                       &//TRIM(Paras(i_Bonds) % line(ind_x(count_a)))
               END IF
               ALLOCATE(Prm % LJs (ip) % g (4))
               Prm % LJs (ip) % g =(/sigma,eps,sigma14,eps14/)
            CASE DEFAULT
               errmsg_f='Unrecognized Lennard-Jones format: '&
                    &//TRIM(Paras(i_Bonds) % line(ind_x(count_a)))
               CALL Add_Error(-1,errmsg_f)
            END SELECT
         END IF
      END DO
      IF(COUNT(oks) /= SIZE(oks)) THEN
         DO n=1,SIZE(oks)
            IF(.NOT. oks(n)) lab0=TRIM(lab0)//' '//TRIM(Prm % Types (n))
         END DO
         errmsg_f='Parameter types '//TRIM(lab0)//' not found in&
              & the LJ parameter list'
         CALL Add_Error(-1,errmsg_f)
         CALL Print_Errors()
      END IF

      m=SIZE(Prm % LJs)*(SIZE(Prm % LJs)+1)/2
      ALLOCATE(Prm % c12 (m),Prm % c6 (m))
      ALLOCATE(Prm % c12_ (m),Prm % c6_ (m))

      DO n=1,SIZE(Prm % LJs)
         DO m=n,SIZE(Prm % LJs)
            sigma=(Prm % LJs (n) % g(1) + Prm % LJs (m) % g(1))*0.5D0
            eps=SQRT(Prm % LJs (n) % g(2)*Prm % LJs (m) % g(2))
            ij=m*(m-1)/2+n
            tmpi=1000.0D0*4.184D0/(unite*avogad)
            Prm % c6(ij) = tmpi*4.0D0*eps*sigma**6
            Prm % c12(ij) = tmpi*4.0D0*eps*sigma**12
            IF(SIZE(Prm % LJs (n) % g) == 4) THEN
               sigma=(Prm % LJs (n) % g(3) + Prm % LJs (m) % g(3))*0.5D0
               eps=SQRT(Prm % LJs (n) % g(4)*Prm % LJs (m) % g(4))
               ij=m*(m-1)/2+n
               tmpi=1000.0D0*4.184D0/(unite*avogad)
               Prm % c6_(ij) = tmpi*4.0D0*eps*sigma**6
               Prm % c12_(ij) = tmpi*4.0D0*eps*sigma**12
            ELSE
               Prm % c6_(ij) = Prm % c6(ij) 
               Prm % c12_(ij) = Prm % c12(ij) 
            END IF
         END DO
      END DO
      CALL Print_Errors()

      DEALLOCATE(LJ,oks)
    END SUBROUTINE Lennard_Jones
    SUBROUTINE Bonds
      LOGICAL :: end_of_list, Dyn
      INTEGER :: n,m_bonds,ia,ib,i_bonds,count_out,count_A,iflags,nlines,nword

      i_Bonds=-1
      DO n=1,SIZE(Paras)
         IF(My_Fxm('BOND', Paras (n) % Residue )) THEN
            i_bonds=n; EXIT
         END IF
      END DO
      IF(i_Bonds == -1) THEN
         errmsg_f='Can''t find bending parameters in the input file'
         CALL Add_Error(-1,errmsg_f)
         RETURN
      END IF

      CALL Start()

      m_bonds=SIZE(Tpg % bonds,2)

      CALL Associate_Type(Paras(i_Bonds) % line, Tpg % bonds, Tpg % atm, count_out)
      CALL Gather_Prm('Bonds ',Tpg % bonds, Tpg % atm, m_bonds, Prm % bonds,  count_out)

    END SUBROUTINE Bonds
    SUBROUTINE Angles
      CHARACTER(len=max_char) :: labs(3)
      LOGICAL :: end_of_list, Dyn
      INTEGER :: n,m_angles,ia,ib,ic,id,i_bonds,count_out&
           &,iflags,nlines,nword
      
      i_Bonds=-1
      DO n=1,SIZE(Paras)
         IF(My_Fxm('ANGLES', Paras (n) % Residue )) THEN
            i_bonds=n; EXIT
         END IF
      END DO
      IF(i_Bonds == -1) THEN
         DO n=1,SIZE(Paras)
            IF(My_Fxm('BENDIN', Paras (n) % Residue )) THEN
               i_bonds=n; EXIT
            END IF
         END DO
      END IF

      IF(i_Bonds == -1) THEN
         errmsg_f='Cant find bending parameters in the input file'
         CALL Add_Error(-1,errmsg_f)
         RETURN
      END IF
      CALL Start()

      m_angles=SIZE(Tpg % angles,2)

      CALL Associate_Type(Paras(i_Bonds) % line, Tpg % angles, Tpg % atm, count_out)

      CALL Gather_Prm('Angles ',Tpg % angles, Tpg % atm, m_angles, Prm % angles,  count_out)

      WRITE(*,*) 'Angle Parameters No. =====>',SIZE(Prm % angles)
    END SUBROUTINE Angles
    SUBROUTINE Dihedrals
      CHARACTER(len=max_char) :: labs(4)
      LOGICAL :: end_of_list, Dyn
      INTEGER :: n,m_dihed,i,ia,ib,ic,id,i_bonds,count_out,o1,o2&
           &,iflags,nlines,nword
      
      i_Bonds=-1
      DO n=1,SIZE(Paras)
         IF(My_Fxm('DIHED', Paras (n) % Residue )) THEN
            i_bonds=n; EXIT
         END IF
      END DO

      DO n=1,SIZE(Paras)
         IF(My_Fxm('TORSIO', Paras (n) % Residue )) THEN
            i_bonds=n; EXIT
         END IF
      END DO

      IF(i_Bonds == -1) THEN
         errmsg_f='Cant find dihedral parameters in the input file'
         CALL Add_Error(-1,errmsg_f)
         CALL Print_Errors()
      END IF
      CALL Start()

      m_dihed=SIZE(Tpg % dihed,2)

      CALL Associate_Type(Paras(i_Bonds) % line, Tpg % dihed, Tpg % atm, count_out)
      CALL Gather_Prm('Torsions ',Tpg % dihed, Tpg % atm, m_dihed, Prm % dihed,  count_out)
      WRITE(*,*) 'Dihedral Parameters No. =====>',SIZE(Prm % dihed)

    END SUBROUTINE Dihedrals

    SUBROUTINE Improper
      CHARACTER(len=max_char) :: labs(4)
      LOGICAL :: end_of_list, Dyn
      INTEGER :: n,m_dihed,i,ia,ib,ic,id,i_bonds,count_out,o1,o2&
           &,iflags,nlines,nword
      
      i_Bonds=-1
      DO n=1,SIZE(Paras)
         IF(My_Fxm('IMPRO', Paras (n) % Residue )) THEN
            i_bonds=n; EXIT
         END IF
      END DO

      IF(i_Bonds == -1) THEN
         errmsg_f='Cant find improper dihedral parameters in the input file'
         CALL Add_Error(-1,errmsg_f)
         CALL Print_Errors()
      END IF
      CALL Start()

      m_dihed=SIZE(Tpg % imph,2)

      CALL Associate_Type(Paras(i_Bonds) % line, Tpg % imph, Tpg % atm, count_out)
      CALL Gather_Prm('Improper Torsions ',Tpg % imph, Tpg % atm, m_dihed, Prm % imph,  count_out)

      WRITE(*,*) 'Improper dihedral Parameters No. =====>',SIZE(Prm % imph)

    END SUBROUTINE Improper
    
    SUBROUTINE Find_Type(ind_in, count_out, Pot, labels)
      INTEGER, INTENT(IN) :: ind_in
      INTEGER, INTENT(INOUT) :: count_out
      TYPE(Par_Pot), DIMENSION(:), INTENT(IN) :: Pot
      CHARACTER(len=max_char), DIMENSION(:), INTENT(IN) :: labels

      INTEGER :: n,m,ma,mb,iflags
      INTEGER, DIMENSION(:), ALLOCATABLE :: ip
      LOGICAL :: ok1,ok2,ok3,ok4
      CHARACTER(len=max_char) :: lab0
      
      ALLOCATE(ip(SIZE(labels)))
      DO n=1,SIZE(labels)
         ip(n)=Type_Number(labels(n))
      END DO

      ok4=.FALSE.
      DO n=1,SIZE(Pot)
         IF(.NOT. ALLOCATED(Pot(n) % p)) CYCLE
         ok1=.TRUE.
         ok2=.TRUE.
         ok3=.FALSE.
         ma=SIZE(Pot (n) % b)
         mb=SIZE(Pot (n) % p)
         
         DO m=1,ma
            ok1=ok1 .AND. (Pot (n) % b (m) == ip(m) .OR. Pot (n) % b (m) == -1)
            ok2=ok2 .AND. (Pot (n) % b (ma-m+1) == ip(m) .OR. Pot (n) % b (ma-m+1) == -1)

            IF(.NOT. (ok1 .OR. ok2)) EXIT
            
            IF(Pot (n) % b (m) == -1) THEN
               ok3=.TRUE.
            END IF
         END DO
         IF(ok1 .OR. ok2) THEN
            IF(ALLOCATED(Params)) DEALLOCATE(Params)
            ALLOCATE(Params(mb+1))
            Params(1:mb)=Pot(n) % p
            Params(mb+1)=DBLE(ind_in)
            IF(ok3) Params(mb+1)=-DBLE(ind_in)
            CALL Add(Params,count_out)
            ok4=.TRUE.
         END IF
      END DO
      IF(.NOT. ok4) THEN
         WRITE(lab0,'(4(a7,2x))') (TRIM(labels(n)),n=1,SIZE(labels))
         errmsg_f='The following bonded parameter is undefined: '//lab0
         CALL Add_Error(-1,errmsg_f)
      END IF
    END SUBROUTINE Find_Type
    SUBROUTINE Associate_Type(lines, Tot_Tpg, atm, count_out)
      TYPE(Atom2Cnt), DIMENSION(:), INTENT(IN) :: atm
      CHARACTER(len=max_char), DIMENSION(:), INTENT(IN) :: lines
      INTEGER, DIMENSION(:,:), INTENT(IN) :: Tot_Tpg
      INTEGER, INTENT(OUT) :: count_out

      INTEGER :: xdim
      INTEGER :: nlines,n,nword,ia,m_Tpg,iflags,i1,count_A
      CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE :: labs
      LOGICAL :: ok,noerror

      xdim=SIZE(Tot_Tpg,1)
      
      ALLOCATE(labs(xdim))

      nlines=SIZE(lines)
      ALLOCATE(p_Tpg(nlines))

      DO n=1,nlines
         CALL MY_parse(lines(n), strngs)
         nword=SIZE(strngs)
         ok=.TRUE.
         DO i1=1,xdim
            IF(Type_Number(strngs(i1),noerror) < -1) ok=.FALSE.
         END DO
         IF(ok) THEN
            ALLOCATE(p_Tpg(n) % b (xdim))
            DO i1=1,xdim
               p_Tpg(n) % b(i1) = Type_Number(strngs(i1))
            END DO
            IF(nword > xdim) THEN
               ALLOCATE(p_Tpg (n) % p (nword-xdim))
               DO ia=1,SIZE(p_Tpg(n) % p)
                  CALL SP_Getnum(strngs(ia+xdim),p_Tpg(n) % p (ia), iflags)
               END DO
            END IF
         END IF
      END DO

      CALL Print_Errors()

      m_Tpg=SIZE(Tot_Tpg,2)
      DO n=1,m_Tpg
         DO i1=1,xdim
            ia=Tot_Tpg(i1,n)
            labs(i1)= atm (ia) % a % betab
         END DO
         CALL Find_Type(n, count_out, p_Tpg, labs)
      END DO
      DEALLOCATE(labs,p_Tpg)
      
    END SUBROUTINE Associate_Type
    FUNCTION Type_Number(label,noerror) RESULT(out)
      CHARACTER(len=*) :: label
      LOGICAL, OPTIONAL :: noerror
      INTEGER  :: out
      CHARACTER(len=max_char) :: lab1
      INTEGER :: i2

      IF(TRIM(label) == 'x' .OR. TRIM(label) == '*') THEN
         out = -1
         RETURN
      END IF
      DO i2=1,SIZE(Prm % Types)
         IF(TRIM(label) == TRIM(Prm % Types(i2))) THEN
            out = i2
            RETURN
         END IF
      END DO
      out = -999
      IF(PRESENT(noerror)) RETURN
      errmsg_f='Atom type '//TRIM(label)//' not on parameter list'
      CALL Add_Error(-1,errmsg_f)
    END FUNCTION Type_Number
    SUBROUTINE Gather_Prm(What_this_is,Tot_Tpg, atm, m_Tpg, Tot_Prm,  count_out)
      CHARACTER(len=*) :: What_this_Is
      TYPE(Atom2Cnt), DIMENSION(:), INTENT(IN) :: atm
      TYPE(Chain_M), DIMENSION(:), ALLOCATABLE, INTENT(OUT)  :: Tot_Prm
      INTEGER, DIMENSION(:,:), INTENT(IN) :: Tot_Tpg
      INTEGER, INTENT(IN)  :: m_Tpg,count_out

      INTEGER :: n,count_a,count_b,m1,m2,i1,ii,i
      TYPE(ChainR8), DIMENSION(:), ALLOCATABLE :: dummy
      INTEGER, DIMENSION(:), ALLOCATABLE :: ind_x, ind_y
      LOGICAL :: end_of_list, Dyn, ok
      CHARACTER(len=max_char) :: lab0, lab1

      ALLOCATE(dummy(count_out), ind_y(count_out), ind_x(m_Tpg))
      ind_x=0
      ind_y=0
      end_of_list=.FALSE.
      count_A=0
      CALL Extract(Params, end_of_list, Dyn, 0)
      DO WHILE(.NOT. end_of_list) 
         count_A=count_A+1
         m1=SIZE(Params)
         ALLOCATE(dummy(count_A) % g(m1-1))
         dummy(count_a) % g =Params(1:m1-1)
         n=INT(Params(m1))
         ind_y(count_a)=n
         ind_x(ABS(n))=ind_x(ABS(n))+1
         CALL Extract(Params, end_of_list, Dyn)
      END DO
      CALL Cleanup()
      count_a=0
      DO ii=1,m_Tpg
         ok=.TRUE.
         DO i=1,ind_x(ii)
            IF(i == 1) count_b=count_a
            count_A=count_A+1
            IF(ind_y(count_a) < 0) ok=.FALSE.
         END DO
         IF(.NOT. ok) THEN
            count_a=count_b
            DO i=1,ind_x(ii)
               count_A=count_A+1
               IF(ind_y(count_a) < 0) THEN
                  ind_y(count_a) = 0
               END IF
            END DO
         END IF
      END DO
      count_a=0
      count_b=0
      DO ii=1,m_Tpg
         IF(ind_x(ii) > 1) THEN
            DO i=1,ind_x(ii)
               count_A=count_A+1
               IF(count_a > SIZE(ind_y)) STOP
               IF(ind_y(count_A) /= 0) count_b=count_b+1
            END DO
         ELSE
            count_A=count_A+1
            count_b=count_b+1
         END IF
      END DO

      ALLOCATE(Tot_Prm (count_b))

      count_a=0
      count_b=0
      DO ii=1,m_Tpg
         IF(ind_x(ii) == 0) THEN
            lab0='!'
            DO i1=1,SIZE(Tot_Tpg,1)
               lab0=TRIM(lab0)//' '//TRIM(atm(Tot_Tpg (i1,ii)) % a % beta)
            END DO
            lab0=TRIM(lab0)//' !'
            lab1='!'
            DO i1=1,SIZE(Tot_Tpg,1)
               lab1=TRIM(lab1)//' '//TRIM(atm(Tot_Tpg (i1,ii)) % a % betab)
            END DO
            lab1=TRIM(lab1)//' !'
            errmsg_f=TRIM(What_This_Is)//TRIM(lab0)//' : '//TRIM(lab1)&
                 &//' not associated with a parameter of the list'
            CALL Add_Error(-1,errmsg_f)
         ELSE
            IF(ind_x(ii) > 1) THEN
               DO i=1,ind_x(ii)
                  count_A=count_A+1
                  IF(ind_y(count_A) /= 0) THEN
                     count_b=count_b+1
                     m1=SIZE(dummy(count_a) % g)
                     ALLOCATE(Tot_Prm (count_b) % g (m1))
                     Tot_Prm(count_b) % g = dummy(count_a) % g
                     Tot_Prm(count_b) % pt = ii
                  END IF
               END DO
            ELSE
               count_A=count_A+1
               count_b=count_b+1
               m1=SIZE(dummy(count_a) % g)
               ALLOCATE(Tot_Prm (count_b) % g (m1-1))
               Tot_Prm(count_b) % g = dummy(count_a) % g(1:m1-1)
               Tot_Prm(count_b) % pt = ii
            END IF
         END IF
      END DO
      DEALLOCATE(ind_x, ind_y, dummy)
    END SUBROUTINE Gather_Prm
  END SUBROUTINE Init
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE Class_SystemPrm
