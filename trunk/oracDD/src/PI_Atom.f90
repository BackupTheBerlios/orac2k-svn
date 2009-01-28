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
MODULE PI_Atom
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Wed Oct  1 2008 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program oracDD ----*

#include 'forces.h'
  USE Print_Defs
  USE Potential
#ifdef HAVE_MPI
  USE mpi
#endif
  USE PI_
  USE Groups
  USE Atom
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f, errmsg_w
  USE IndBox
  USE Forces, Forces_Memory=>Memory
  USE Neighbors
  USE Cell, ONLY: oc,co, Volume
  IMPLICIT none
  PRIVATE
  PUBLIC PI_Atom_,xpg,ypg,zpg,xp0,yp0,zp0,chg,gmass,Id,Slv,Grp_No,grppt,maplg&
       &,Mapnl,Maps,natom,ngroup,PI_Atom__Neigh_,List,Neigha__&
       &,Neigha_,Update_, PI_Atom_Indexes

  TYPE Neigha__
     INTEGER :: no=0
     INTEGER, DIMENSION (:), ALLOCATABLE :: nb
  END TYPE Neigha__
  TYPE Neigha_
     TYPE(Neigha__), DIMENSION(:), ALLOCATABLE :: Neigh
  END TYPE Neigha_
  TYPE(Neigha_), DIMENSION (:), ALLOCATABLE, TARGET, SAVE :: List


  REAL(8), ALLOCATABLE, SAVE :: xpg(:),ypg(:),zpg(:),xp0(:),yp0(:),zp0(:)
  REAL(8), DIMENSION(:), ALLOCATABLE, SAVE :: chg,gmass,fppx,fppy,fppz
  INTEGER, ALLOCATABLE, SAVE :: Id(:),Slv(:),grppt(:,:),grp_No(:)
  LOGICAL, ALLOCATABLE, SAVE :: maplg(:)
  TYPE :: Mapnl
     INTEGER, ALLOCATABLE :: ex(:)
  END type Mapnl
  TYPE(Mapnl), ALLOCATABLE, SAVE :: Maps(:)
  INTEGER, SAVE :: No_Calls=0,natom=0,ngroup=0
  INTEGER, SAVE :: nnx=20,nny=20,nnz=20
  
CONTAINS
  FUNCTION PI_Atom_() RESULT(out)
    LOGICAL :: out
    INTEGER, POINTER :: Index_0(:)

    IF(No_Calls /= 0) THEN
       DEALLOCATE(xpg,ypg,zpg,xp0,yp0,zp0,chg,Id,Slv,grppt,Grp_No,Maps,maplg,gmass)
    END IF

    CALL Labelling
    CALL Memory
    CALL Gather_Atoms

    No_Calls=No_Calls+1
    out=.TRUE.
  CONTAINS
    SUBROUTINE Labelling
      INTEGER :: ngroupa,nc,ncutoff
      INTEGER, ALLOCATABLE :: Ind(:)
      REAL(8) :: rcut
      INTEGER :: nx,ny,nz,numcell,l,p,q,r,numcell1,nj
      
      IF(Reordering) THEN
         WRITE(kprint, '('' --->   Reordering <---- '')')
         ncutoff=SIZE(Radii)
         rcut=Radii(ncutoff) % out + Radii(ncutoff) % update
         IF(.NOT. Neighbors_(rcut, nnx, nny, nnz)) CALL Print_Errors()
         
         IF(.NOT. Neighbors__Particles(Groupa(:) % xa, Groupa(:) % ya,&
              & Groupa(:) % za, Groupa(:) % knwn)) CALL Print_Errors()
         
         ngroupa=COUNT(Groupa(:) % knwn /= 0)
         ALLOCATE(Ind(ngroupa))
         nc=0
         DO nx=0,nnx-1
            DO ny=0,nny-1
               DO nz=0,nnz-1
                  numcell=nz+nnz*(ny+nny*nx)+1
                  l=Head_xyz(numcell)
                  DO WHILE(l > 0)
                     p=Chain_xyz(l) % i
                     q=Chain_xyz(l) % j
                     r=Chain_xyz(l) % k
                     numcell1=r+nnz*(q+nny*p)+1

                     nc=nc+1
                     Ind(nc)=l

                     l=Chain_xyz(l) % p
                  END DO
               END DO
            END DO
         END DO
         IF(.NOT. IndBoxL_(Groupa(:) % knwn,Groupa(:) % AtSt,Groupa(:) %&
              & AtEn, Ind)) CALL Print_Errors()

      ELSE
         IF(.NOT. IndBox_(Groupa(:) % knwn,Groupa(:) % AtSt,Groupa(:) %&
              & AtEn)) CALL Print_Errors()
      END IF
    END SUBROUTINE Labelling
!!$
!!$--- Get Memory
!!$
    SUBROUTINE Memory

      natom=SIZE(IndBox_a_t) ; ngroup=SIZE(IndBox_g_t)

      IF(ALLOCATED(xpg)) THEN
         DEALLOCATE(xpg,ypg,zpg,grppt)
         DEALLOCATE(xp0,yp0,zp0,chg,Id,Slv,Maps,maplg,gmass,Grp_No)
      END IF

      ALLOCATE(xpg(ngroup),ypg(ngroup),zpg(ngroup),grppt(2,ngroup))
      
      ALLOCATE(xp0(natom),yp0(natom),zp0(natom),chg(natom),Id(natom)&
           &,Slv(natom),Maps(natom),maplg(natom),gmass(natom)&
           &,Grp_No(natom))
    END SUBROUTINE Memory
!!$
!!$--- Gather Atoms to the CPU box
!!$
    SUBROUTINE Gather_Atoms
      INTEGER :: n,m,p,q,nn,count0,g1,g2
      INTEGER :: Timesa=0

      DO n=1,ngroup
         m=IndBox_g_t(n)
         xpg(n)=Groupa(m) % xa
         ypg(n)=Groupa(m) % ya
         zpg(n)=Groupa(m) % za
         grppt(1,n)=Groupa(m) % AtSt
         grppt(2,n)=Groupa(m) % AtEn
      END DO
      DO n=1,natom
         m=IndBox_a_t(n)
         xp0(n)=Atoms(m) % xa
         yp0(n)=Atoms(m) % ya
         zp0(n)=Atoms(m) % za
         chg(n)=Atoms(m) % chg
         Id(n)=Atoms(m) % Id_Type
         Slv(n)=Atoms(m) % Id_Slv
         Grp_No(n)=Atoms(m) % Grp_No
         gmass(n)=Atoms(m) % pmass
      END DO
      IF(Timesa == 1) THEN
!!$         WRITE(200+PI_Node_Cart,'(i8,3x,3e17.8)') (IndBox_g_p(n)&
!!$              &,Groupa(IndBox_g_p(n)) % xa,Groupa(IndBox_g_p(n)) % ya&
!!$              & ,Groupa(IndBox_g_p(n)) % za,n=1,SIZE(IndBox_g_p))
!!$         WRITE(200+PI_Node_Cart,'(2i8,3x,3e17.8)')&
!!$              & (Groupa(IndBox_g_t(n)) % knwn, IndBox_g_t(n)&
!!$              & ,xpg(n),ypg(n),zpg(n),n=1,ngroup)
      END IF
      Timesa=Timesa+1

      ALLOCATE(Index_0(SIZE(Atoms)))
      Index_0=-1
      DO n=1,natom
         m=IndBox_a_t(n)
         Index_0(m)=n
      END DO
!!$
!!$--- Change grppt to reset atoms known in the box
!!$
      DO n=1,ngroup
         g1=grppt(1,n)
         g2=grppt(2,n)
         IF(Index_0(g1) /= -1) THEN
            grppt(1,n)=Index_0(g1)
         ELSE
            WRITE(*,*) 'Weired!!!'
            STOP
         END IF
         IF(Index_0(g2) /= -1) THEN
            grppt(2,n)=Index_0(g2)
         ELSE
            WRITE(*,*) 'Weired!!!'
            STOP
         END IF
      END DO
      
!!$
!!$--- Change exclusion atoms 
!!$
      
      DO n=1,natom
         m=IndBox_a_t(n)
         nn=SIZE(Atoms_Tpg(m) % Ex)
         count0=0
         DO q=1,nn
            p=Atoms_Tpg(m) % Ex(q)
            IF(Index_0(p) /= -1) THEN
               count0=count0+1
            END IF
         END DO
         IF(count0 == 0) CYCLE
         ALLOCATE(Maps(n) % ex(count0))
         count0=0
         DO q=1,nn
            p=Atoms_Tpg(m) % Ex(q)
            IF(Index_0(p) /= -1) THEN
               count0=count0+1
               Maps(n) % Ex(count0)=Index_0(Atoms_Tpg(m) % Ex(q))
            END IF
         END DO
      END DO
      
    END SUBROUTINE Gather_Atoms

!!$
!!$--- Compute Forces
!!$


  END FUNCTION PI_Atom_
  FUNCTION Update_() RESULT(out)
    LOGICAL :: out
    CALL Gather_Update
    out=.TRUE.
  CONTAINS
!!$
!!$--- Gather Atoms to the CPU box
!!$
    SUBROUTINE Gather_Update
      INTEGER :: n,m,p,q
      INTEGER :: Timesa=0

      DO n=1,ngroup
         m=IndBox_g_t(n)
         IF(Groupa(m) % knwn == 0) CYCLE != no need to update atoms outside range
         xpg(n)=Groupa(m) % xa
         ypg(n)=Groupa(m) % ya
         zpg(n)=Groupa(m) % za
         DO p=grppt(1,n),grppt(2,n)
            q=IndBox_a_t(p)
            xp0(p)=Atoms(q) % xa
            yp0(p)=Atoms(q) % ya
            zp0(p)=Atoms(q) % za
         END DO
      END DO
    END SUBROUTINE Gather_Update

!!$
!!$--- Compute Forces
!!$


  END FUNCTION Update_
  SUBROUTINE PI_Atom_Indexes(nx,ny,nz)
    INTEGER :: nx,ny,nz
    nx=nnx; ny=nny; nz=nnz
  END SUBROUTINE PI_Atom_Indexes
  FUNCTION PI_Atom__Neigh_(Keep) RESULT(out)
    LOGICAL, OPTIONAL :: Keep
    LOGICAL :: out

    REAL(8) :: rcut,rcut2
    INTEGER :: ncutoff,vect0(1)
    INTEGER, POINTER :: vect(:)
    TYPE(Neigha__), DIMENSION(:), POINTER :: Neigha
    INTEGER :: ig,ii,n,p,q,r,iv,jv,kv,l,o,count0,nx,ny,nz,numcell&
         &,counter,No_Nei,knw_ig,knw_l,n0,numcell1
    REAL(8) :: xpgi,ypgi,zpgi,xpgj,ypgj,zpgj,xa,ya,za,X_PBC,Y_PBC&
         &,Z_PBC,xc,yc,zc,rsq
    INTEGER, ALLOCATABLE :: nei(:),known(:),nei0(:)
    LOGICAL, ALLOCATABLE :: okk(:)
    LOGICAL :: ok

    out=.FALSE.
    IF(ngroup == 0) THEN
       out=.TRUE.
       errmsg_f='PI_Atom_ must be called before this function can be u&
            &sed'
       CALL Add_Errors(-1,errmsg_f)       
       RETURN
    END IF

    ALLOCATE(nei(ngroup),known(ngroup))


    ncutoff=SIZE(Radii)
    rcut=Radii(ncutoff) % out + Radii(ncutoff) % update
    rcut2=rcut**2
    IF(.NOT. PRESENT(Keep)) THEN
       IF(ALLOCATED(List)) DEALLOCATE(List)
       ALLOCATE(List(ncutoff))
       DO n=1,ncutoff
          ALLOCATE(List(n) % Neigh(ngroup))
       END DO
    END IF
    Neigha=>List(ncutoff) % Neigh

    IF(.NOT. Neighbors_(rcut, nnx, nny, nnz)) CALL Print_Errors()
    IF(.NOT. Neighbors__Particles(xpg,ypg,zpg)) CALL Print_Errors()

    counter=0
    known(:)=Groupa(IndBox_g_t(:)) % knwn
    
    ALLOCATE(okk(ngroup))
    okk=.FALSE.
    DO ig=1,ngroup
       knw_ig=known(ig)
       xpgi=xpg(ig)
       ypgi=ypg(ig)
       zpgi=zpg(ig)
       p=Chain_xyz(ig) % i
       q=Chain_xyz(ig) % j
       r=Chain_xyz(ig) % k
       IF(ALLOCATED(Neigha(ig) % nb)) DEALLOCATE(Neigha(ig) % nb)
       Neigha(ig) % no=0
       count0=0
       numcell1=r+nnz*(q+nny*p)+1
       DO o=1,SIZE(Ind_xyz)
          iv=Ind_xyz(o) % i
          jv=Ind_xyz(o) % j
          kv=Ind_xyz(o) % k
          ok=.FALSE.
          IF(iv == 0 .AND. jv == 0 .AND. kv == 0) ok=.TRUE.
          nx=mod(mod(p+iv,nnx)+nnx,nnx)
          ny=mod(mod(q+jv,nny)+nny,nny)
          nz=mod(mod(r+kv,nnz)+nnz,nnz)
          numcell=nz+nnz*(ny+nny*nx)+1
          l=Head_xyz(numcell)
          DO WHILE(l > 0)
             IF(ok .AND. l < ig) THEN
                l=Chain_xyz(l) % p
                CYCLE
             END IF
             
             knw_l=known(l)+Knw_ig
             xpgj=xpg(l)
             ypgj=ypg(l)
             zpgj=zpg(l)
             xa=xpgj-xpgi
             ya=ypgj-ypgi
             za=zpgj-zpgi
             X_PBC=_PBC(xa)
             Y_PBC=_PBC(ya)
             Z_PBC=_PBC(za)
             xa=xa+X_PBC
             ya=ya+Y_PBC
             za=za+Z_PBC
             xc=co(1,1)*xa+co(1,2)*ya+co(1,3)*za
             yc=           co(2,2)*ya+co(2,3)*za
             zc=                      co(3,3)*za
             rsq=xc*xc+yc*yc+zc*zc
             IF(rsq <= rcut2 .AND. knw_l /= 4) THEN
                count0=count0+1
                nei(count0)=l
                IF(known(ig) == 2) THEN
                   okk(ig)=.TRUE.
                END IF
                IF(known(l) == 2) THEN
                   okk(l)=.TRUE.
                END IF
             END IF
             l=Chain_xyz(l) % p
          END DO
       END DO
       Neigha(ig) % no=count0
       IF(count0 /= 0) THEN
          Neigha(ig) % no=count0
          ALLOCATE(Neigha(ig) % nb(count0))
          Neigha(ig) % nb=nei(1:count0)
       END IF
       counter=counter+count0
    END DO
    WRITE(*,*) PI_Node_Cart,COUNT(okk),COUNT(Known==2)
    CALL MPI_ALLREDUCE(counter,No_Nei,1,MPI_INTEGER,MPI_SUM,PI_Comm_Cart,ierr)
    IF(PI_Node_Cart == 0) THEN
       WRITE(kprint,*) 'Neighbors (',No_Nei,')'
    END IF
  END FUNCTION PI_Atom__Neigh_
END MODULE PI_Atom
