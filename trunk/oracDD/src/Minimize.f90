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
MODULE Minimize
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Dec 11 2008 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*

#include "config.h"

  USE Print_Defs
  USE PI_
  USE Potential
  USE Units, ONLY: efact
  USE IndIntraBox
  USE Forces, ONLY: Force, fp=>fp_n0
  USE PI_Communicate
  USE PI_Collectives
  USE Atom
  USE Intra, ONLY: Intra_n0_,uconstr_slt,uconstr_slv
  USE PI_IntraMaps, ONLY: IntraMaps_n0_, IntraMaps_n1_
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f, errmsg_w
  USE IndBox
  IMPLICIT none
  PRIVATE
  PUBLIC Minimize__Bonds_
CONTAINS
  FUNCTION Minimize__Bonds_() RESULT(out)
    LOGICAL :: out
    REAL(8), ALLOCATABLE :: x(:),g(:)
    REAL(8), ALLOCATABLE :: d(:),Gold(:),w(:)
    REAL(8) :: tlev,Energy,eps
    INTEGER :: iflag,n,nn
    INTEGER, SAVE :: iprint(2)=(/-1,0/),method=2,irest=0,step_Max=100
    LOGICAL :: finish,ok
    INTEGER :: natom,ntot,icall,i

    WRITE(kprint,100) 
    out=.TRUE.
    eps=Rattle__Param % adjust_eps

    natom=SIZE(Atoms)
    ntot=3*natom
    ALLOCATE(x(ntot),g(ntot),d(ntot),Gold(ntot),w(ntot))

    CALL PI__ResetSecondary
    CALL IntraMaps_n0_
    CALL PI__ShiftIntra(1,_INIT_EXCHANGE_)
    IF(.NOT. IndIntraBox_n0_()) CALL Print_Errors()
    CALL Intra_n0_(_INIT_EXCHANGE_,_CONSTRS_)
    Energy=(uconstr_Slv+uconstr_Slt)

    icall=1
20  CONTINUE
    CALL PI__ShiftIntra(1,_INIT_EXCHANGE_)
    CALL Intra_n0_(_EXCHANGE_ONLY_,_CONSTRS_)
    Energy=(uconstr_Slv+uconstr_Slt)
    CALL Arrays_Get

30  CONTINUE
    CALL cgfam(ntot,x,Energy,g,d,Gold,iprint,eps,w,iflag,irest,method,finish)
    CALL Arrays_Put

    IF(iflag < 0 .OR. iCall > 1000) THEN
       errmsg_f='Stopped adjusting bond lengths. No convergence.'
       CALL Add_Errors(-1,errmsg_f)
       out=.FALSE.
       RETURN
    END IF
    IF(iflag == 1) THEN
       icall=icall+1
       IF(iCall > step_max) GOTO 50
       GOTO 20
    ELSE IF(iflag == 2) THEN
             
!!$
!!$ Termination Test.  The user may replace it by some other test. However, 
!!$ the parameter 'FINISH' must be set to 'TRUE' when the test is satisfied.
!!$

       tlev=eps*(1.0_8 + dabs(Energy))
       i=0
40     i=i+1
       IF(i > ntot) THEN
          Finish=.TRUE.
          GOTO 30
       END IF
       IF(ABS(g(i)) > Tlev) THEN
          GOTO 30
       ELSE
          GOTO 40
       END IF       
    END IF
50  CONTINUE

    CALL Arrays_Put
    IF(.NOT. Atom__Convert(_X_TO_XA_)) CALL Print_Errors()

    WRITE(kprint,200) Energy/DBLE(SIZE(Indx_Constr,2)+SIZE(Indx_Bonds,2))
100 FORMAT(10x,'=======>  Adjusting Bonds <======')
200 FORMAT('===> Bond Average RMS deviation from force field bond leng&
         &th ',e14.5,' A^2 ') 
  CONTAINS
    SUBROUTINE Arrays_Get
      REAL(8), POINTER :: xc(:),yc(:),zc(:)
      INTEGER :: n,m
      xc=>Atoms(:) % x
      yc=>Atoms(:) % y
      zc=>Atoms(:) % z
      
      CALL PI_AllGather_(xc,yc,zc)
      DO n=1,natom
         m=3*(n-1)
         x(m+1)=xc(n)
         x(m+2)=yc(n)
         x(m+3)=zc(n)
      END DO
      xc=>fp(:) % x
      yc=>fp(:) % y
      zc=>fp(:) % z
      CALL PI_AllGather_(xc,yc,zc)
      DO n=1,natom
         m=3*(n-1)
         g(m+1)=-xc(n)
         g(m+2)=-yc(n)
         g(m+3)=-zc(n)
      END DO
    END SUBROUTINE Arrays_Get
    SUBROUTINE Arrays_Put
      REAL(8), POINTER :: xc(:),yc(:),zc(:)
      INTEGER :: n,m
      xc=>Atoms(:) % x
      yc=>Atoms(:) % y
      zc=>Atoms(:) % z
      
      DO n=1,natom
         m=3*(n-1)
         xc(n)=x(m+1)
         yc(n)=x(m+2)
         zc(n)=x(m+3)
      END DO
    END SUBROUTINE Arrays_Put

  END FUNCTION Minimize__Bonds_
END MODULE Minimize
