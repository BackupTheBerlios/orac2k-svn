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
MODULE LennardJones

!!$***********************************************************************
!!$   Time-stamp: <2007-01-11 18:44:55 marchi>                           *
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

  USE Parameters
  USE Units
  USE TypesPrm
  USE Constants
  USE Resid
  USE Strings
  USE MyParse
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f
  USE Print_Defs
  IMPLICIT none
  PRIVATE
  PUBLIC  LennardJones_, LennardJones__Type, LennardJones__SigmaEps&
       &, LennardJones__Read, LennardJones__Write, LennardJones__Par
  TYPE :: LennardJones__SigmaEps
     INTEGER :: pt 
     REAL(8), DIMENSION(:), ALLOCATABLE :: g     
  END TYPE LennardJones__SigmaEps

  TYPE :: LennardJones__Type
     TYPE(LennardJones__SigmaEps), ALLOCATABLE :: Par_SE(:)
     REAL(8), ALLOCATABLE :: c12(:),c12_(:),c6(:),c6_(:)     
  END TYPE LennardJones__Type

  TYPE LennardJones__Chain
     CHARACTER(len=max_char) :: lab
     REAL(8), DIMENSION(:), ALLOCATABLE :: g
  END TYPE LennardJones__Chain

  TYPE(LennardJones__Type), SAVE, TARGET :: LennardJones__Par
  
CONTAINS
  FUNCTION LennardJones_() RESULT(out)
    TYPE(LennardJones__Type), POINTER :: out
    INTEGER :: n,m,ij,ntypes,i,i_bonds,count_out,count_A,iflags&
         &,nlines,nword,m1,ip
    LOGICAL :: ok,noerror
    INTEGER, DIMENSION(:), ALLOCATABLE :: ind_x
    REAL(8) :: eps,sigma,eps14,sigma14,tmpi
    TYPE(LennardJones__Chain), DIMENSION(:), ALLOCATABLE :: LJ
    CHARACTER(len=max_char) :: lab0
    REAL(8), DIMENSION(:), ALLOCATABLE  :: Params
    LOGICAL, DIMENSION(:), ALLOCATABLE  :: oks
    
    out=>NULL()
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
       CALL Add_Errors(-1,errmsg_f); out=>NULL(); RETURN
    END IF
       
    ALLOCATE(ind_x(SIZE(Paras(i_Bonds) % line)))
    ALLOCATE(LennardJones__Par % Par_SE (SIZE(TypesPrm__Types)))
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
    ALLOCATE(oks(SIZE(LennardJones__Par % Par_SE)))
    
    oks=.FALSE.
    DO n=1,count_a
       ip=TypesPrm__Number(LJ(n) % lab, noerror)
       IF(ip > 0) THEN
          oks(ip)=.TRUE.
          LennardJones__Par % Par_SE (ip) % pt = ip
          
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
             ALLOCATE(LennardJones__Par % Par_SE (ip) % g (2))
             LennardJones__Par % Par_SE (ip) % g =(/sigma,eps/)
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
                out=>NULL()
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
                out=>NULL()
             END IF
             ALLOCATE(LennardJones__Par % Par_SE (ip) % g (4))
             LennardJones__Par % Par_SE (ip) % g =(/sigma,eps,sigma14,eps14/)
          CASE DEFAULT
             errmsg_f='Unrecognized Lennard-Jones format: '&
                  &//TRIM(Paras(i_Bonds) % line(ind_x(count_a)))
             CALL Add_Errors(-1,errmsg_f)
             out=>NULL()                
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
       out=>NULL()
    END IF
    
    m=SIZE(LennardJones__Par % Par_SE)
    m=m*(m+1)/2
    ALLOCATE(LennardJones__Par % c12 (m),LennardJones__Par % c6 (m))
    ALLOCATE( LennardJones__Par % c12_ (m),LennardJones__Par % c6_ (m))
    
    DO n=1,SIZE(LennardJones__Par % Par_SE)
       DO m=n,SIZE(LennardJones__Par % Par_SE)
          sigma=(LennardJones__Par % Par_SE (n) % g(1) &
               &+ LennardJones__Par % Par_SE (m) % g(1))*0.5D0
          eps=SQRT(LennardJones__Par % Par_SE (n) % g(2)&
               &*LennardJones__Par % Par_SE (m) % g(2))
          ij=m*(m-1)/2+n
          tmpi=1000.0D0*4.184D0/(unite*avogad)
          LennardJones__Par % c6(ij) = tmpi*4.0D0*eps*sigma**6
          LennardJones__Par % c12(ij) = tmpi*4.0D0*eps*sigma**12
          IF((SIZE(LennardJones__Par % Par_SE (n) % g) == 4)&
               & .AND. (SIZE(LennardJones__Par % Par_SE (m) % g) == 4)) THEN
             sigma=(LennardJones__Par % Par_SE (n) % g(3) &
                  &+ LennardJones__Par % Par_SE (m) % g(3))*0.5D0
             eps=SQRT(LennardJones__Par % Par_SE (n) % g(4)&
                  &*LennardJones__Par % Par_SE (m) % g(4))
             ij=m*(m-1)/2+n
             tmpi=1000.0D0*4.184D0/(unite*avogad)
             LennardJones__Par % c6_(ij) = tmpi*4.0D0*eps*sigma**6
             LennardJones__Par % c12_(ij) = tmpi*4.0D0*eps*sigma**12
          ELSE
             LennardJones__Par % c6_(ij) = LennardJones__Par % c6(ij) 
             LennardJones__Par % c12_(ij) = LennardJones__Par % c12(ij) 
          END IF
       END DO
    END DO
    out=>LennardJones__Par
    WRITE(kprint,*) 'Total Lennard-Jones Parameter &
         &types No. =====>',SIZE(LennardJones__Par % C12)
    DEALLOCATE(LJ,oks)
  CONTAINS
    FUNCTION Valid() RESULT(out)
      LOGICAL :: out
      out=.FALSE.
      IF(ALLOCATED(TypesPrm__Types)) out=.TRUE.
    END FUNCTION Valid
    
  END FUNCTION LennardJones_
  SUBROUTINE LennardJones__Write
    INTEGER :: n,o_c12,o_c12_,o_c6,o_c6_,o_SE

    o_c12=0
    o_c12_=0
    o_c6=0
    o_c6_=0
    o_SE=0

    IF(ALLOCATED(LennardJones__Par % c12)) o_c12=SIZE(LennardJones__Par % c12)
    IF(ALLOCATED(LennardJones__Par % c12_)) o_c12_=SIZE(LennardJones__Par % c12_)
    IF(ALLOCATED(LennardJones__Par % c6)) o_c6=SIZE(LennardJones__Par % c6)
    IF(ALLOCATED(LennardJones__Par % c6_)) o_c6_=SIZE(LennardJones__Par % c6_)
    IF(ALLOCATED(LennardJones__Par % Par_SE)) o_SE=SIZE(LennardJones__Par % Par_SE)
    WRITE(kbinary) o_c12,o_c12_,o_c6,o_c6_,o_SE
    IF(o_c12 /= 0) WRITE(kbinary) LennardJones__Par % c12
    IF(o_c12_ /= 0) WRITE(kbinary) LennardJones__Par % c12_
    IF(o_c6 /= 0) WRITE(kbinary) LennardJones__Par % c6
    IF(o_c6_ /= 0) WRITE(kbinary) LennardJones__Par % c6_

    DO n=1,o_SE
       IF(ALLOCATED(LennardJones__Par % Par_SE (n) % g)) THEN
          WRITE(kbinary) SIZE(LennardJones__Par % Par_SE (n) % g)
          WRITE(kbinary) LennardJones__Par % Par_SE (n) % pt&
               &, LennardJones__Par % Par_SE (n) % g
       ELSE
          WRITE(kbinary) 0
       END IF
    END DO
  END SUBROUTINE LennardJones__Write
  FUNCTION LennardJones__Read() RESULT(out)
    TYPE(LennardJones__Type), POINTER :: out
    INTEGER :: n,m,o_c12,o_c12_,o_c6,o_c6_,o_SE

    out=>NULL()

    READ(kbinary,ERR=100,END=200) o_c12,o_c12_,o_c6,o_c6_,o_SE
    IF(o_c12 /= 0) ALLOCATE(LennardJones__Par % c12(o_c12))
    IF(o_c12_ /= 0) ALLOCATE(LennardJones__Par % c12_ (o_c12_))
    IF(o_c6 /= 0) ALLOCATE(LennardJones__Par % c6 (o_c6))
    IF(o_c6_ /= 0) ALLOCATE(LennardJones__Par % c6_ (o_c6_))

    IF(o_c12 /= 0) READ(kbinary,ERR=100,END=200) LennardJones__Par % c12
    IF(o_c12_ /= 0) READ(kbinary,ERR=100,END=200) LennardJones__Par % c12_
    IF(o_c6 /= 0) READ(kbinary,ERR=100,END=200) LennardJones__Par % c6
    IF(o_c6_ /= 0) READ(kbinary,ERR=100,END=200) LennardJones__Par % c6_

    IF(o_SE /= 0) ALLOCATE(LennardJones__Par % Par_SE (o_SE))
    DO n=1,o_SE
       READ(kbinary,ERR=100,END=200) m
       IF(m /= 0) THEN
          ALLOCATE(LennardJones__Par % Par_SE (n) % g(m))
          READ(kbinary,ERR=100,END=200) LennardJones__Par % Par_SE (n) % pt&
               &, LennardJones__Par % Par_SE (n) % g
       END IF
    END DO
    out=>LennardJones__Par
    RETURN
100 errmsg_f='Error while reading Lennard-Jones Parameters'
    CALL Add_Errors(-1,errmsg_f)
    RETURN
200 errmsg_f='End of file found while reading Lennard-Jones Parameters'
    CALL Add_Errors(-1,errmsg_f)
    RETURN    
  END FUNCTION LennardJones__Read
     
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*
     
END MODULE LennardJones
