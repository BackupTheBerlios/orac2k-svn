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
#include 'forces.h'
SUBROUTINE Lists_(NShell0)
!!$***********************************************************************
!!$   Time-stamp: <2007-01-24 10:48:13 marchi>                           *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Mon Sep 29 2008 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This routine is part of the program  oracDD ----*
!!$---- This routine is contained in the Module Direct ----*

  INTEGER :: NShell0
  INTEGER :: ii,n,p,q,r,o,iv,jv,kv,nx,ny,nz,numcell,count_g&
         &,count_gs,ngrp_j,ngrp_js,count0,i
  REAL(8) :: xpgi,ypgi,zpgi,xpgj,ypgj,zpgj,xgg,ygg,zgg,rcutb,rcutb2&
       &,rcutc,rcutc2,xa,ya,za,xc,yc,zc,rsq

  INTEGER :: AtSt_i,AtEn_i,ig,j,k,l,i1,jj,j1&
       &,count_b,count_c,p_mapb,count_cc,count_a
  TYPE(Neigha__), DIMENSION(:), POINTER :: Neigha,Neighb,Neighc
  INTEGER, ALLOCATABLE :: inda(:),indb(:)
  INTEGER :: NShell,npp


  NShell=Nshell0-2
  IF(NShell == 1) RETURN
  IF(NShell > SIZE(Radii)) RETURN
  Neigha=>List(NShell) % Neigh

  SELECT CASE(NShell)
  CASE(3)
     Neighb=>List(2) % Neigh
     rcutb=Radii(2) % out + Radii(2) % update 
     rcutb2=rcutb**2

     Neighc=>List(1) % Neigh
     rcutc=Radii(1) % out + Radii(1) % update 
     rcutc2=rcutc**2
  CASE(2)
     Neighb=>List(1) % Neigh
     rcutb=Radii(1) % out + Radii(1) % update 
     rcutb2=rcutb**2

     Neighc=>NULL()
  END SELECT

  DO ig=1,ngroup
     xpgi=xpg(ig)
     ypgi=ypg(ig)
     zpgi=zpg(ig)
     count_g=0
     count_gs=0
     IF(NShell == 3) THEN
        IF(ALLOCATED(Neighb(ig) % nb)) DEALLOCATE(Neighb(ig) % nb)
        Neighb(ig) % no = 0 
        IF(ALLOCATED(Neighc(ig) % nb)) DEALLOCATE(Neighc(ig) % nb)
        Neighc(ig) % no = 0 
     ELSE IF(Nshell == 2) THEN
        IF(ALLOCATED(Neighb(ig) % nb)) DEALLOCATE(Neighb(ig) % nb)
        Neighb(ig) % no = 0 
     END IF
     count_b=0
     count_cc=0
     DO o=1,Neigha(ig) % no
        l=Neigha(ig) % nb(o)
        xpgj=xpg(l)
        ypgj=ypg(l)
        zpgj=zpg(l)
        xa=xpgi-xpgj
        ya=ypgi-ypgj
        za=zpgi-zpgj
        xa=xa+_PBC(xa)
        ya=ya+_PBC(ya)
        za=za+_PBC(za)
        xc=co(1,1)*xa+co(1,2)*ya+co(1,3)*za
        yc=           co(2,2)*ya+co(2,3)*za
        zc=                      co(3,3)*za
        rsq=xc*xc+yc*yc+zc*zc
        IF(NShell == 3) THEN
           IF(rsq <= rcutb2) THEN
              count_b=count_b+1
              neib(count_b)=l
           END IF
           IF(rsq <= rcutc2) THEN
              count_cc=count_cc+1
              neic(count_cc)=l
           END IF
        ELSE IF(NShell == 2) THEN
           IF(rsq <= rcutb2) THEN
              count_b=count_b+1
              neib(count_b)=l
           END IF
        END IF
     END DO
     IF(count_b /= 0) THEN
        Neighb(ig) % no=count_b
        ALLOCATE(Neighb(ig) % nb(count_b))
        Neighb(ig) % nb=neib(1:count_b)
     END IF
     IF(count_cc /= 0) THEN
        Neighc(ig) % no=count_cc
        ALLOCATE(Neighc(ig) % nb(count_cc))
        Neighc(ig) % nb=neic(1:count_cc)
     END IF
  END DO
  count_a=SUM(Neigha(:) % no)
  count_b=SUM(Neighb(:) % no)
  count_c=SUM(Neighc(:) % no)
  
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,count_a,1,MPI_INTEGER4,MPI_SUM&
       &,PI_Comm_Cart,ierr)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,count_b,1,MPI_INTEGER4,MPI_SUM&
       &,PI_Comm_Cart,ierr)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,count_c,1,MPI_INTEGER4,MPI_SUM&
       &,PI_Comm_Cart,ierr)
  WRITE(kprint,*) count_a,count_b,count_c
END SUBROUTINE Lists_
