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
MODULE PrmUtilities

!!$***********************************************************************
!!$   Time-stamp: <2007-01-11 18:12:56 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Jan 11 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*

  USE Constants
  USE AtomCnt
  USE MyParse
  USE Node
  USE TypesPrm
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f, errmsg_w
  USE STRPAK
  IMPLICIT none

  PRIVATE
  PUBLIC :: SystemPrm__Chain, SystemPrm__Chain2, PrmUtilities__Type, PrmUtilities__Gather
  TYPE :: SystemPrm__Chain
     INTEGER :: pt 
     REAL(8), ALLOCATABLE :: g(:)
  END TYPE SystemPrm__Chain
  TYPE :: SystemPrm__Chain2
     INTEGER :: pt 
     REAL(8), POINTER :: g(:)
  END TYPE SystemPrm__Chain2
  TYPE PrmUtilities__Chain
     REAL(8), DIMENSION(:), ALLOCATABLE :: g
  END TYPE PrmUtilities__Chain
  TYPE :: Par_Pot
     INTEGER, ALLOCATABLE :: b(:)
     REAL(8), ALLOCATABLE :: p(:)
  END TYPE Par_Pot
  TYPE(Par_Pot), ALLOCATABLE :: p_Tpg(:) 
  REAL(8),  POINTER  :: Params(:)
CONTAINS
  FUNCTION PrmUtilities__Type(lines, Tot_Tpg) RESULT(out)
    LOGICAL :: out
    CHARACTER(len=max_char), DIMENSION(:), INTENT(IN) :: lines
    INTEGER, DIMENSION(:,:), INTENT(IN) :: Tot_Tpg
    
    INTEGER :: xdim
    INTEGER :: nlines,n,nword,ia,m_Tpg,iflags,i1,count_A
    CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE :: labs
    LOGICAL :: ok,noerror

    out=.TRUE.
    xdim=SIZE(Tot_Tpg,1) ; ALLOCATE(labs(xdim))
    nlines=SIZE(lines); ALLOCATE(p_Tpg(nlines))
    
    DO n=1,nlines
       nword=MYparse_(lines(n))
       ok=.TRUE.
       DO i1=1,xdim
          IF(TypesPrm__Number(strngs(i1),noerror) < -1) ok=.FALSE.
       END DO
       IF(ok) THEN
          ALLOCATE(p_Tpg(n) % b (xdim))
          DO i1=1,xdim
             p_Tpg(n) % b(i1) = TypesPrm__Number(strngs(i1))
          END DO
          IF(nword > xdim) THEN
             ALLOCATE(p_Tpg (n) % p (nword-xdim))
             DO ia=1,SIZE(p_Tpg(n) % p)
                CALL SP_Getnum(strngs(ia+xdim),p_Tpg(n) % p (ia), iflags)
                IF(iflags /= 0) THEN
                   errmsg_f='Error while conventing Parameter file&
                        & lines to real numbers: '&
                        &//TRIM(strngs(ia+xdim))//' Not translated'
                   CALL Add_Errors(-1,errmsg_f)
                   out=.FALSE.
                   RETURN
                END IF
             END DO
          END IF
       END IF
    END DO
    IF(.NOT. Node_()) STOP

    m_Tpg=SIZE(Tot_Tpg,2)
    DO n=1,m_Tpg
       DO i1=1,xdim
          ia=Tot_Tpg(i1,n)
          labs(i1)=AtomCnts(ia) % betab
       END DO
       IF(.NOT. Find_Type(n, p_Tpg, labs)) THEN
          out=.FALSE.
          RETURN
       END IF
    END DO
    DEALLOCATE(labs,p_Tpg)
  CONTAINS
    FUNCTION Find_Type(ind_in, Pot, labels) RESULT(out)
      LOGICAL :: out
      INTEGER, INTENT(IN) :: ind_in
      TYPE(Par_Pot), DIMENSION(:), INTENT(IN) :: Pot
      CHARACTER(len=max_char), DIMENSION(:), INTENT(IN) :: labels

      INTEGER :: n,m,ma,mb,iflags
      INTEGER, ALLOCATABLE :: ip(:)
      LOGICAL :: ok1,ok2,ok3,ok4
      CHARACTER(len=max_char) :: lab0
      
      out=.TRUE.
      ALLOCATE(ip(SIZE(labels)))
      DO n=1,SIZE(labels)
         ip(n)=TypesPrm__Number(labels(n))
      END DO

      ok4=.FALSE.
      DO n=1,SIZE(Pot)
         IF(.NOT. ALLOCATED(Pot(n) % p)) CYCLE
         ok1=.TRUE.
         ok2=.TRUE.
         ok3=.FALSE.
         ma=SIZE(Pot (n) % b)
         mb=SIZE(Pot (n) % p)
!!$
!!$--- Pot(n) % b(m) == -1 when wild cards in the parameters list
!!$
         DO m=1,ma
            ok1=ok1 .AND. (Pot (n) % b (m) == ip(m) .OR. Pot (n) % b (m) == -1)
            ok2=ok2 .AND. (Pot (n) % b (ma-m+1) == ip(m) .OR. Pot (n) % b (ma-m+1) == -1)

            IF(.NOT. (ok1 .OR. ok2)) EXIT

!!$-- When wild cards on the parameter list make sure only one 
!!$-- parameter is counted. Only applicable for torsions and impropers

            IF(Pot (n) % b (m) == -1) THEN
               ok3=.TRUE.
            END IF
         END DO
         IF(ok1 .OR. ok2) THEN
            IF(ASSOCIATED(Params)) DEALLOCATE(Params)
            ALLOCATE(Params(mb+1))
            Params(1:mb)=Pot(n) % p
            Params(mb+1)=DBLE(ind_in)
            IF(ok3) Params(mb+1)=-DBLE(ind_in)
            CALL Node__Push(Params)
            ok4=.TRUE.
         END IF
      END DO
      IF(.NOT. ok4) THEN
         WRITE(lab0,'(4(a7,2x))') (TRIM(labels(n)),n=1,SIZE(labels))
         errmsg_f='The following bonded parameter is undefined: '//lab0
         CALL Add_Errors(-1,errmsg_f)
         out=.FALSE.
         RETURN
      END IF
      DEALLOCATE(ip)
    END FUNCTION Find_Type    
  END FUNCTION PrmUtilities__Type
  FUNCTION PrmUtilities__Gather(What_this_is,Tot_Tpg, Prm) RESULT(out)
    LOGICAL :: out
    CHARACTER(len=*) :: What_this_Is
    TYPE(SystemPrm__Chain2), POINTER :: Prm(:)
    INTEGER, INTENT(IN) :: Tot_Tpg(:,:)
    
    INTEGER :: n,count_a,count_b,m1,m2,i1,ii,i,count_out,m_Tpg
    TYPE(PrmUtilities__Chain), ALLOCATABLE :: dummy(:)
    INTEGER,  ALLOCATABLE :: ind_x(:), ind_y(:)
    LOGICAL :: end_of_list, Dyn, ok
    CHARACTER(len=max_char) :: lab0, lab1
    
    out=.TRUE.
    count_out=Node__Size() 
    IF(count_out == 0) THEN
       errmsg_w='Warning: No '//TRIM(What_This_is)//' were found: Is this correct?'
       CALL Add_Errors(1,errmsg_w)
       RETURN
    END IF
    m_Tpg=SIZE(Tot_Tpg,2)
    ALLOCATE(dummy(count_out), ind_y(count_out), ind_x(m_Tpg))
    ind_x=0
    ind_y=0
    count_A=0
    DO WHILE(Node__Pop(Params))
       count_A=count_A+1
       m1=SIZE(Params)
       ALLOCATE(dummy(count_A) % g(m1-1))
       dummy(count_a) % g =Params(1:m1-1)
       n=INT(Params(m1))
       ind_y(count_a)=n
       ind_x(ABS(n))=ind_x(ABS(n))+1
    END DO
    IF(.NOT. Node_()) STOP

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
    
    ALLOCATE(Prm(count_b))
    
    count_a=0
    count_b=0
    DO ii=1,m_Tpg
       IF(ind_x(ii) == 0) THEN
          lab0='!'
          DO i1=1,SIZE(Tot_Tpg,1)
             lab0=TRIM(lab0)//' '//TRIM(AtomCnts(Tot_Tpg (i1,ii)) % beta)
          END DO
          lab0=TRIM(lab0)//' !'
          lab1='!'
          DO i1=1,SIZE(Tot_Tpg,1)
             lab1=TRIM(lab1)//' '//TRIM(AtomCnts(Tot_Tpg (i1,ii)) % beta)
          END DO
          lab1=TRIM(lab1)//' !'
          errmsg_f=TRIM(What_This_Is)//TRIM(lab0)//' : '//TRIM(lab1)&
               &//' not associated with a parameter of the list'
          CALL Add_Errors(-1,errmsg_f)
          out=.FALSE.
       ELSE
          IF(ind_x(ii) > 1) THEN
             DO i=1,ind_x(ii)
                count_A=count_A+1
                IF(ind_y(count_A) /= 0) THEN
                   count_b=count_b+1
                   m1=SIZE(dummy(count_a) % g)
                   ALLOCATE(Prm (count_b) % g (m1))
                   Prm(count_b) % g = dummy(count_a) % g
                   Prm(count_b) % pt = ii
                END IF
             END DO
          ELSE
             count_A=count_A+1
             count_b=count_b+1
             m1=SIZE(dummy(count_a) % g)
             ALLOCATE(Prm (count_b) % g (m1-1))
             Prm(count_b) % g = dummy(count_a) % g(1:m1-1)
             Prm(count_b) % pt = ii
          END IF
       END IF
    END DO
    DEALLOCATE(ind_x, ind_y, dummy)
  END FUNCTION PrmUtilities__Gather
  
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE PrmUtilities
