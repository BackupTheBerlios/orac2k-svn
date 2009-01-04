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
MODULE Rattle
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Dec 12 2008 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*

  USE Potential
  USE Node
  USE SystemPrm
  USE SystemTpg
  USE IndBox
  USE IndIntraBox, ONLY: Param_Constr, Indx_Constr
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f, errmsg_w
  IMPLICIT none
  PRIVATE
  PUBLIC Init_, Parameters_,Verlet_,Correct_

  REAL(8), ALLOCATABLE, SAVE :: xp0(:),yp0(:),zp0(:),xp1(:),yp1(:)&
       &,zp1(:),vpx(:),vpy(:),vpz(:),mass0(:)
  INTEGER, SAVE :: natom,nc
  TYPE :: NList
     INTEGER, ALLOCATABLE :: indx(:)
  END type NList
  TYPE :: ClusterList0
     INTEGER :: n0
     INTEGER, ALLOCATABLE :: n1(:)
  END type ClusterList0
  TYPE :: ClusterList
     INTEGER :: n0
     INTEGER, ALLOCATABLE :: n1(:)
     REAL(8), ALLOCATABLE :: g(:)
  END type ClusterList
  TYPE :: Rattle__Type1
     INTEGER, ALLOCATABLE :: n1(:,:)
  END type Rattle__Type1
  TYPE :: Rattle__Type2
     REAL(8), ALLOCATABLE :: g(:)
  END type Rattle__Type2
  TYPE(ClusterList), ALLOCATABLE, SAVE :: cnst0(:)
  TYPE(Rattle__Type1), ALLOCATABLE, TARGET, SAVE :: cnst(:)
  TYPE(Rattle__Type2), ALLOCATABLE, TARGET, SAVE :: dss(:),coeff(:)
  LOGICAL, SAVE :: Call_Init=.FALSE.
  REAL(8), ALLOCATABLE ::   mat(:,:),gam(:),gamo(:),matx(:),maty(:)&
       &,matz(:),xc(:),yc(:),zc(:),dd(:)
  INTEGER, ALLOCATABLE :: ipiv(:)
  INTEGER, SAVE :: nmat
CONTAINS
  FUNCTION Init_(ntot) RESULT(out)
    LOGICAL :: out
!!$
!!$--- Arrange the list of constraints to a list of group constrains
!!$
    INTEGER :: ntot,ncluster
    INTEGER :: ncc,n,m,n0,nn,i1,i2,count0,q,p,count_out,count_A,ia,ib&
         &,g 
    INTEGER, SAVE :: nmax=10
    INTEGER, ALLOCATABLE :: ind0(:),ind1(:,:),cn(:)
    INTEGER, POINTER :: via(:)
    TYPE(NList), ALLOCATABLE :: cnlist(:)
    TYPE(ClusterList0), ALLOCATABLE :: cnst_d(:)
    LOGICAL, ALLOCATABLE :: oks0(:),oks(:)
    LOGICAL :: oka

    Call_Init=.TRUE.
    out=.TRUE.
    IF(.NOT. Rattle__Param % switch) RETURN
    nc=SIZE(Prm % Constr)

    ALLOCATE(ind0(ntot),ind1(nmax,ntot))

!!$
!!$-- Construct a raw list with multiple occurencies 
!!$

    ind0=0
    DO nn=1,nc
       n=Prm % Constr(nn) % pt
       i1=Tpg % Bonds(1,n)
       i2=Tpg % Bonds(2,n)
       ind0(i1)=ind0(i1)+1
       ind0(i2)=ind0(i2)+1
       ind1(ind0(i1),i1)=nn
       ind1(ind0(i2),i2)=-nn
    END DO

    ALLOCATE(cnlist(ntot),oks0(SIZE(Tpg % Bonds,2)),oks(ntot))

!!$
!!$--- Construct the final list by eliminating multiple occurencies
!!$--- each element of the list points to prm % Constr not Tpg % Bonds
!!$

    DO n=1,ntot
       m=ind0(n)
       count0=0       
       oks0(ABS(conv0_(ind1(1:m,n))))=.TRUE.
       DO q=1,m
          g=ABS(conv0_(ind1(q,n)))
          IF(oks0(g)) THEN
             count0=count0+1
             oks0(g)=.FALSE.
          END IF
       END DO
       ALLOCATE(cnlist(n) % indx(count0))
       oks0(ABS(conv0_(ind1(1:m,n))))=.TRUE.
       count0=0
       DO q=1,m
          p=ind1(q,n)
          g=ABS(conv0_(p))
          IF(oks0(g)) THEN
             count0=count0+1
             oks0(g)=.FALSE.
             cnlist(n) % indx(count0)=p
          END IF
       END DO
    END DO

!!$
!!$--- Get the number of constraint clusters
!!$

    ALLOCATE(cnst_d(ntot))
    oks=.TRUE.
    ncluster=0
    DO n=1,ntot
       IF(oks(n)) THEN
          IF(.NOT. Node_()) STOP
          CALL Next_Connection(n)
          count_out=Node__Size()
          IF(count_out == 0 ) CYCLE

          ncluster=ncluster+1
          ALLOCATE(cn(count_out))
          count_A=0
          cn=0
          DO WHILE(Node__Pop(via)) 
             oka=.TRUE.
             ia=via(1)

             DO p=count_A,1,-1
                i1=cn(p)
                IF(ia == i1) THEN
                   oka=.FALSE.
                   EXIT
                END IF
             END DO
             IF(oka) THEN
                count_A=count_A+1
                cn(count_A)=ia
             END IF
          END DO

          cnst_d(ncluster) % n0 = count_A
          ALLOCATE(cnst_d(ncluster) % n1(count_A))
          cnst_d(ncluster) % n1=cn(1:count_A)
          DEALLOCATE(cn)
       END IF
    END DO

!!$
!!$--- Collect constraints of each cluster. cnst0 (n) % n1(m) points 
!!$--- to Tpg % Bonds and cnst0 (n) % g(m) is the constraint distance
!!$

    ALLOCATE(cnst0(ncluster))
    DO n=1,ncluster
       n0=cnst_d(n) % n0
       cnst0(n) % n0 = n0
       ALLOCATE(cnst0(n) % n1(n0))
       ALLOCATE(cnst0(n) % g(n0))
       DO m=1,n0
          nn=cnst_d(n) % n1(m)
          cnst0(n) % n1(m)=Prm % Constr(nn) % pt
          cnst0(n) % g(m) = Prm % Constr(nn) % g(1)
       END DO
    END DO
    nmat=Rattle__Param % mim_Max
    ALLOCATE(mat(nmat,nmat),gam(nmat),gamo(nmat),matx(2*nmat&
         &*nmat),maty(2*nmat*nmat),matz(2*nmat*nmat),xc(nmat)&
         &,yc(nmat),zc(nmat),dd(nmat),ipiv(2*nmat))

    out=.TRUE.

  CONTAINS
    RECURSIVE SUBROUTINE Next_Connection(n)
      INTEGER :: n
      INTEGER :: m,i,ia
      INTEGER :: vi0
      LOGICAL :: oka

      IF(.NOT. oks(n)) RETURN
      oks(n)=.FALSE.

      m=SIZE(cnlist(n) % indx)

      vi0=0
      DO i=1,m
         ib=conv0_(cnlist(n) % indx(i))
         IF(ib > 0) ia=Tpg % Bonds(2,ib)
         IF(ib < 0) ia=Tpg % Bonds(1,-ib)
         CALL Next_Connection(ia)
         vi0=ABS(cnlist(n) % indx(i))
         CALL Node__Push(vi0)
      END DO
    END SUBROUTINE Next_Connection
    ELEMENTAL FUNCTION conv0_(n) RESULT(out)
      INTEGER :: out
      INTEGER, INTENT(in) :: n
      IF(n > 0) out=Prm % Constr(n) % pt
      IF(n < 0) out=-Prm % Constr(-n) % pt
    END FUNCTION conv0_
  END FUNCTION Init_

  FUNCTION Parameters_(mass,knwn) RESULT(out)
    LOGICAL :: out
    REAL(8) :: mass(:)
    INTEGER :: knwn(:)
    INTEGER :: n,n0,p,m,ia,ib,na,nb,ma,mb,nn
    TYPE(ClusterList), ALLOCATABLE:: cnst1(:)
    INTEGER, ALLOCATABLE :: indx(:)

    out=.TRUE.
    IF(.NOT. Rattle__Param % switch) RETURN
    IF(.NOT. Call_Init) THEN
       out=.FALSE.
       errmsg_f='Internal: In Rattle Parameters_() must be called afte&
            &r Init_()'
       CALL Add_Errors(-1,errmsg_f)
       RETURN
    END IF
    IF(ALLOCATED(cnst)) DEALLOCATE(cnst)
    IF(ALLOCATED(dss)) DEALLOCATE(dss)
    IF(ALLOCATED(coeff)) DEALLOCATE(coeff)
    
    ALLOCATE(indx(SIZE(cnst0)))

    nc=0
    DO n=1,SIZE(cnst0)
       n0=cnst0(n) % n0
       DO m=1,n0
          p=cnst0(n) % n1(m)
          ia=Tpg % Bonds(1,p)
          ib=Tpg % Bonds(2,p)
          IF(knwn(ia) == 1 .AND. knwn(ib) == 1) THEN
             nc=nc+1
             indx(nc)=n
             EXIT
          ELSE IF(knwn(ia) == 1 .NEQV. knwn(ib) == 1) THEN
             WRITE(*,*) 'Wow '
             STOP
          END IF
       END DO
    END DO


    ALLOCATE(cnst(nc),dss(nc),coeff(nc))

    DO nn=1,nc
       n=indx(nn)
       n0=cnst0(n) % n0
       ALLOCATE(cnst(nn) % n1(2,n0))
       ALLOCATE(dss(nn) % g(n0))
       ALLOCATE(coeff(nn) % g(n0))
       DO m=1,n0
          p=cnst0(n) % n1(m)
          na=Tpg % Bonds(1,p)
          nb=Tpg % Bonds(2,p)
          ma=BoxInd_a_p(na)
          mb=BoxInd_a_p(nb)
          cnst(nn) % n1(1,m) = ma
          cnst(nn) % n1(2,m) = mb
          dss(nn) % g(m)=cnst0(n) % g(m)**2
          coeff(nn) % g(m)=2.0_8*(mass(na)+mass(nb))/(mass(na)&
               &*mass(nb))
       END DO
    END DO

    IF(ALLOCATED(xp0)) THEN
       DEALLOCATE(xp0,yp0,zp0,xp1,yp1,zp1,vpx,vpy,vpz,mass0)
    END IF

    natom=SIZE(IndBox_a_p)

    ALLOCATE(xp0(natom),yp0(natom),zp0(natom),xp1(natom),yp1(natom)&
         &,zp1(natom),vpx(natom),vpy(natom),vpz(natom),mass0(natom))

    DO nn=1,natom
       n=IndBox_a_t(IndBox_a_p(nn))
       mass0(nn)=1.0_8/mass(n)
    END DO
    out=.TRUE.
  END FUNCTION Parameters_
  INCLUDE "Rattle__Verlet.f90"
  INCLUDE "Rattle__Correct.f90"
END MODULE Rattle
