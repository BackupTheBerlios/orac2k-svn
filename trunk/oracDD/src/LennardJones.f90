!!$/---------------------------------------------------------------------\
!!$                                                                      |
!!$  Copyright (C) 2006-2007 Massimo Marchi <Massimo.Marchi@cea.fr>      |
!!$                                                                      |
!!$      This program is free software;  you  can  redistribute  it      |
!!$      and/or modify it under the terms of the GNU General Public      |
!!$      License version 2 as published  by  the  Free  Software         |
!!$      Foundation;                                                     |
!!$                                                                      |
!!$      This program is distributed in the hope that  it  will  be      |
!!$      useful, but WITHOUT ANY WARRANTY; without even the implied      |
!!$      warranty of MERCHANTABILITY or FITNESS  FOR  A  PARTICULAR      |
!!$      PURPOSE.   See  the  GNU  General  Public License for more      |
!!$      details.                                                        |
!!$                                                                      |
!!$      You should have received a copy of the GNU General  Public      |
!!$      License along with this program; if not, write to the Free      |
!!$      Software Foundation, Inc., 59  Temple  Place,  Suite  330,      |
!!$      Boston, MA  02111-1307  USA                                     |
!!$                                                                      |
!!$\---------------------------------------------------------------------/
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

  USE Units
  USE TypesPrm
  USE Constants
  USE Resid
  USE Strings
  USE MyParse
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f
  IMPLICIT none
  PRIVATE
  PUBLIC  LennardJones_, LennardJones__Type, LennardJones__SigmaEps
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

  TYPE(LennardJones__Type), SAVE :: LennardJones__Par
  
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
    WRITE(*,*) 'Total Lennard-Jones Parameter &
         &types No. =====>',SIZE(LennardJones__Par % C12)
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
   
