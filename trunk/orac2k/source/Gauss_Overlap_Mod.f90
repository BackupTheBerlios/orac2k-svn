MODULE Gauss_Overlap_Mod
!!$***********************************************************************
!!$   Time-stamp: <2005-02-07 11:14:40 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Feb  4 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

  USE Neighbors_Mod

  TYPE Indx
     INTEGER, DIMENSION(:), POINTER :: ind_g
     INTEGER, DIMENSION(:), POINTER :: ind_a
  END TYPE Indx

  TYPE Ind_Dip
     INTEGER :: nato
     INTEGER :: ngrp
     REAL(8), DIMENSION (:), POINTER :: pol
     REAL(8), DIMENSION (:), POINTER :: mx,my,mz
     REAL(8), DIMENSION (:), POINTER :: x,y,z
  END TYPE Ind_Dip
  TYPE Overlap
     INTEGER :: nato
     REAL(8), DIMENSION (:), POINTER :: x,y,z
  END TYPE Overlap

  TYPE(Indx), SAVE, PRIVATE :: g
  TYPE(Neighbors), DIMENSION (:), SAVE, POINTER :: nn_Pol
  TYPE(Ind_Dip), SAVE :: dip

  REAL(8), PRIVATE, SAVE :: alpha
  REAL(8), DIMENSION (:,:), POINTER, PRIVATE, SAVE :: coa,oca

CONTAINS
  SUBROUTINE Gauss_init(nato, ngrp, grppt, index, nnlpp, x, y, z, co,&
       & oc, dipoles, polar, beta)

!!$---- This subroutine is part of the program ORAC ----*

    IMPLICIT none
    
    INTEGER :: ngrp, nato, grppt(2,*)
    INTEGER :: index(ngrp),nnlpp(*)
    REAL(8), DIMENSION (3,3), TARGET :: co, oc
    REAL(8), TARGET :: polar(*)
    REAL(8), TARGET :: dipoles(3,nato),x(nato),y(nato),z(nato) 
    REAL(8) :: beta

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: i, ii, count, m
    
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    dip%mx=>dipoles(1,1:nato)
    dip%my=>dipoles(2,1:nato)
    dip%mz=>dipoles(3,1:nato)
    dip%pol=>polar(1:nato)
    dip%x=>x(1:nato)
    dip%y=>y(1:nato)
    dip%z=>z(1:nato)
    coa=>co; oca=>oc

    CALL Neighbors_init(nnlpp,ngrp)
    nn_Pol=>nn; 

    alpha=beta**2; dip%ngrp=ngrp; dip%nato=nato

    ALLOCATE(G%ind_g(ngrp))
    ALLOCATE(G%ind_a(ngrp))

    count=0
    DO ii=1,ngrp
       G%ind_g(ii)=-1
       G%ind_a(ii)=-1
       IF(index(ii) == 2) THEN
          count=count+1

          i=grppt(1,ii)
          G%ind_g(ii)=ii
          G%ind_a(ii)=i
       END IF
    END DO

  END SUBROUTINE Gauss_init

  TYPE(Overlap) FUNCTION Gauss_overlap(R_Cut)

!!$---- This subroutine is part of the program ORAC ----*

    IMPLICIT none
    REAL(8) :: R_Cut

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: i, ii, j, jj, m, o, ig, jg, im
    REAL(8) :: xpi,ypi,zpi,xg,yg,zg,xc,yc,zc,rsq,R_Cut2,U_ene,alpha2&
         &,aux
    REAL(8) :: mu_x,mu_y,mu_z,mmu_x,mmu_y,mmu_z,mmu
    REAL(8), DIMENSION (:), POINTER :: fpx,fpy,fpz
    LOGICAL, SAVE :: first_entry=.TRUE.
    
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    IF(first_entry) THEN
       NULLIFY(Gauss_Overlap%x,Gauss_Overlap%y,Gauss_Overlap%z)
       first_entry=.FALSE.
    END IF

    IF(.NOT. ASSOCIATED(Gauss_Overlap%x)) THEN
       ALLOCATE(Gauss_Overlap%x(dip%nato),Gauss_Overlap%y(dip%nato)&
            &,Gauss_Overlap%z(dip%nato))
       WRITE(*,*) 'Ippa '
    END IF

    Gauss_overlap%nato=dip%nato
    fpx=>Gauss_Overlap%x
    fpy=>Gauss_Overlap%y
    fpz=>Gauss_Overlap%z

    fpx=0.0D0;fpy=0.0D0;fpz=0.0D0

    R_Cut2=R_Cut**2
    alpha2=alpha*2.0D0
    U_ene=0.0D0
    DO ig=1,dip%ngrp
       i=G%ind_a(ig)
       ii=G%ind_g(ig)
       IF(i < 0) CYCLE
       xpi=dip%x(i)
       ypi=dip%y(i)
       zpi=dip%z(i)

       mu_x=dip%mx(i)
       mu_y=dip%my(i)
       mu_z=dip%mz(i)

       m=nn_Pol(ii)%m
       DO im=1,m
          jg=nn_Pol(i)%Neigh(im)
          j=G%ind_a(jg)
          jj=G%ind_g(jg)
          IF(j < 0) CYCLE

          mmu_x=mu_x*dip%mx(j)
          mmu_y=mu_y*dip%my(j)
          mmu_z=mu_z*dip%mz(j)
          mmu=mmu_x+mmu_y+mmu_z

          xg=xpi-dip%x(j)
          yg=ypi-dip%y(j)
          zg=zpi-dip%z(j)
          xg=xg-PBC(xg)
          yg=yg-PBC(yg)
          zg=zg-PBC(zg)
          xc=coa(1,1)*xg+coa(1,2)*yg+coa(1,3)*zg
          yc=            coa(2,2)*yg+coa(2,3)*zg
          zc=                        coa(3,3)*zg
          rsq=xc*xc+yc*yc+zc*zc
          IF(rsq < r_Cut2) THEN
             aux=mmu*DEXP(-alpha2*rsq)
             U_ene=U_ene+aux
             fpx(i)=fpx(i)-2.0D0*alpha2*xc*aux
             fpy(i)=fpy(i)-2.0D0*alpha2*yc*aux
             fpz(i)=fpz(i)-2.0D0*alpha2*zc*aux
          END IF
       END DO
    END DO
    WRITE(*,*) U_ene
  CONTAINS
    REAL(8) FUNCTION PBC(x)
      REAL(8) :: x
      PBC=DNINT(0.5D0*x)
    END FUNCTION PBC
  END FUNCTION Gauss_overlap
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE Gauss_Overlap_Mod
