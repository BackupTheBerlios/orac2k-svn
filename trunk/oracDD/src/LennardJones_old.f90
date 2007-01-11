MODULE LennardJones

!!$***********************************************************************
!!$   Time-stamp: <2007-01-11 14:28:54 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Wed Jan 10 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*

  USE Units
  USE TypesPrm
  USE Constants
  USE Resid
  USE Strings
  USE MyParse
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f
  IMPLICIT none
  PRIVATE
  PUBLIC  LennardJones_
  TYPE :: LennardJones__Type
     INTEGER :: pt 
     REAL(8), DIMENSION(:), ALLOCATABLE :: g
  END TYPE LennardJones__Type
  TYPE LennardJones__Chain
     CHARACTER(len=max_char) :: lab
     REAL(8), DIMENSION(:), ALLOCATABLE :: g
  END TYPE LennardJones__Chain
  
  TYPE(LennardJones__Type), ALLOCATABLE :: LennardJones__Par(:)
  REAL(8), ALLOCATABLE, SAVE :: LennardJones__c12(:),LennardJones__c12_(:)&
       &,LennardJones__c6(:),LennardJones__c6_(:)
CONTAINS
  FUNCTION LennardJones_() RESULT(out)
    LOGICAL :: out
    INTEGER :: n,m,ij,ntypes,i,i_bonds,count_out,count_A,iflags&
         &,nlines,nword,m1,ip
    LOGICAL :: ok,noerror
    INTEGER, DIMENSION(:), ALLOCATABLE :: ind_x
    REAL(8) :: eps,sigma,eps14,sigma14,tmpi
    TYPE(LennardJones__Chain), DIMENSION(:), ALLOCATABLE :: LJ
    CHARACTER(len=max_char) :: lab0
    REAL(8), DIMENSION(:), ALLOCATABLE  :: Params
    LOGICAL, DIMENSION(:), ALLOCATABLE  :: oks
    
    out=.FALSE.
    IF(.NOT. Valid()) THEN
       errmsg_f='Atomic Types must be defined before Lennard Jones&
            & Parameters can be set'
       CALL Add_Errors(-1,errmsg_f)
       RETURN
    END IF
    
    i_Bonds=-1
    DO n=1,SIZE(Paras)
       IF(My_Fxm('NONBOND', Paras (n) % Residue )) THEN; i_bonds=n; 
          EXIT
       END IF
    END DO
    IF(i_Bonds == -1) THEN
       errmsg_f='Can''t find atom type labels in the input file'
       CALL Add_Errors(-1,errmsg_f); out=.FALSE.; RETURN
    END IF
       
    ALLOCATE(ind_x(SIZE(Paras(i_Bonds) % line)))
    ALLOCATE(LennardJones__Par (SIZE(TypesPrm__Types)))
    ALLOCATE(LJ(SIZE(Paras(i_Bonds) % line)))
    
    count_a=0
    DO n=1,SIZE(Paras(i_Bonds) % line)
       nword=MYParse_(Paras(i_Bonds) % line (n))
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
    ALLOCATE(oks(SIZE(LennardJones__Par)))
    
    oks=.FALSE.
    DO n=1,count_a
       ip=TypesPrm__Number(LJ(n) % lab, noerror)
       IF(ip > 0) THEN
          oks(ip)=.TRUE.
          LennardJones__Par (ip) % pt = ip
          
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
                CALL Add_Errors(-1,errmsg_f)
             END IF
             ALLOCATE(LennardJones__Par (ip) % g (2))
             LennardJones__Par (ip) % g =(/sigma,eps/)
          CASE(6)
             sigma=LJ(n) % g (3)*LJ_fact
             IF(LJ(n) % g (1) == 0.0D0) THEN
                eps=-LJ(n) % g (2)
             ELSE IF(LJ(n) % g (2) == 0.0D0) THEN
                eps=LJ(n) % g (1)
             ELSE
                errmsg_f='Unrecognized Lennard-Jones format: '&
                     &//TRIM(Paras(i_Bonds) % line(ind_x(count_a)))
                CALL Add_Errors(-1,errmsg_f)
                out=.FALSE.
             END IF
             sigma14=LJ(n) % g (6)*LJ_fact
             IF(LJ(n) % g (4) == 0.0D0) THEN
                eps14=-LJ(n) % g (5)
             ELSE IF(LJ(n) % g (5) == 0.0D0) THEN
                eps14=LJ(n) % g (4)
             ELSE
                errmsg_f='Unrecognized Lennard-Jones format: '&
                     &//TRIM(Paras(i_Bonds) % line(ind_x(count_a)))
                CALL Add_Errors(-1,errmsg_f)
                out=.FALSE.
             END IF
             ALLOCATE(LennardJones__Par (ip) % g (4))
             LennardJones__Par (ip) % g =(/sigma,eps,sigma14,eps14/)
          CASE DEFAULT
             errmsg_f='Unrecognized Lennard-Jones format: '&
                  &//TRIM(Paras(i_Bonds) % line(ind_x(count_a)))
             CALL Add_Errors(-1,errmsg_f)
             out=.FALSE.
          END SELECT
       END IF
    END DO
    IF(COUNT(oks) /= SIZE(oks)) THEN
       DO n=1,SIZE(oks)
          IF(.NOT. oks(n)) lab0=TRIM(lab0)//' '//TRIM(TypesPrm__Types (n))
       END DO
       errmsg_f='Parameter types '//TRIM(lab0)//' not found in&
            & the LJ parameter list'
       CALL Add_Errors(-1,errmsg_f)
       CALL Print_Errors()
       out=.FALSE.
    END IF
    
    m=SIZE(LennardJones__Par)
    m=m*(m+1)/2
    ALLOCATE(LennardJones__c12 (m),LennardJones__c6 (m))
    ALLOCATE( LennardJones__c12_ (m),LennardJones__c6_ (m))
    
    DO n=1,SIZE(LennardJones__Par)
       DO m=n,SIZE(LennardJones__Par)
          sigma=(LennardJones__Par (n) % g(1) + LennardJones__Par (m) % g(1))*0.5D0
          eps=SQRT(LennardJones__Par (n) % g(2)*LennardJones__Par (m) % g(2))
          ij=m*(m-1)/2+n
          tmpi=1000.0D0*4.184D0/(unite*avogad)
          LennardJones__c6(ij) = tmpi*4.0D0*eps*sigma**6
          LennardJones__c12(ij) = tmpi*4.0D0*eps*sigma**12
          IF((SIZE(LennardJones__Par (n) % g) == 4)&
               & .AND. (SIZE(LennardJones__Par (m) % g) == 4)) THEN
             sigma=(LennardJones__Par (n) % g(3) + LennardJones__Par (m) % g(3))*0.5D0
             eps=SQRT(LennardJones__Par (n) % g(4)*LennardJones__Par (m) % g(4))
             ij=m*(m-1)/2+n
             tmpi=1000.0D0*4.184D0/(unite*avogad)
             LennardJones__c6_(ij) = tmpi*4.0D0*eps*sigma**6
             LennardJones__c12_(ij) = tmpi*4.0D0*eps*sigma**12
          ELSE
             LennardJones__c6_(ij) = LennardJones__c6(ij) 
             LennardJones__c12_(ij) = LennardJones__c12(ij) 
          END IF
       END DO
    END DO
    WRITE(*,*) 'Total Lennard-Jones Parameter types No. =====>',SIZE(LennardJones__C12)
    
    DEALLOCATE(LJ,oks)
  CONTAINS
    FUNCTION Valid() RESULT(out)
      LOGICAL :: out
      out=.FALSE.
      IF(ALLOCATED(TypesPrm__Types)) out=.TRUE.
    END FUNCTION Valid
    
  END FUNCTION LennardJones_
     
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*
     
END MODULE LennardJones
   
