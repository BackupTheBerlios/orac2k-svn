MODULE HYDRATION_Mod
!!$***********************************************************************
!!$   Time-stamp: <2008-03-10 14:44:05 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Thu Aug  4 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

  USE Module_Neighbors, ONLY: Neigh_Start=>Start, Neigh_Delete&
       &=>Delete, neigha, Neighbors
  USE xerror_mod
  INTEGER, SAVE :: nstart,nend,nlocal,ngrp,nbun,node=0,nprocs=1&
       &,nlocal_neigh
  INTEGER, SAVE :: n_neighbors=1,khydration=0,ncx=10,ncy=10,ncz=10&
       &,n_write=1,kbinary=0
  REAL(8), SAVE :: coeff=1.0D0,cutoff_max=4.0D0
  LOGICAL, SAVE :: hydration=.FALSE.
  LOGICAL, DIMENSION (:), ALLOCATABLE, SAVE :: tags_sv

  TYPE(Neighbors), DIMENSION (:), POINTER, SAVE :: neigh_s
  INTEGER, DIMENSION (:), ALLOCATABLE, SAVE :: index_st,index_sv
  INTEGER, SAVE :: indxyz,nind
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE ::  indxi,indxj,indxk
  INTEGER, DIMENSION (:,:), ALLOCATABLE, SAVE :: mres,grppt
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE ::  resg

CONTAINS

!!$---- This subroutine is part of the program ORAC ----*
  SUBROUTINE Initialize(index_a,index_b,khydration_a)
    IMPLICIT none
    
!!$----------------------------- ARGUMENTS ------------------------------*
    
    INTEGER :: index_a(*),index_b(*),khydration_a
    
!!$----------------------------- VARIABLES ------------------------------*

    INTEGER :: i,m

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    khydration=khydration_a
    m=index_a(1)
    ALLOCATE(index_st(m+1))
    index_st(1)=m
    DO i=1,m
       index_st(i+1)=index_a(i+1)
    END DO
    m=index_b(1)
    ALLOCATE(index_sv(m+1))
    index_sv(1)=m
    DO i=1,m
       index_sv(i+1)=index_b(i+1)
    END DO

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*
  END SUBROUTINE Initialize
  SUBROUTINE Initialize_P(node_a,nprocs_a,ngrp_a,nbun_a)
    IMPLICIT none
    
!!$----------------------------- ARGUMENTS ------------------------------*

    INTEGER :: node_a,nprocs_a,ngrp_a,nbun_a
    INTEGER :: ntot
    
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    nbun=nbun_a
    ngrp=ngrp_a
    nprocs=nprocs_a
    node=node_a
    ntot=index_st(1)/nprocs
    nstart=node*ntot+1
    IF(node /= nprocs-1) THEN
       nend=(node+1)*ntot
    ELSE
       nend=index_st(1)
    END IF
    nlocal=nend-nstart

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*
    
  END SUBROUTINE Initialize_P
  SUBROUTINE Initialize_Array(grppt_a,mres_a,resg_a)
    IMPLICIT none
    
!!$----------------------------- ARGUMENTS ------------------------------*

    INTEGER :: grppt_a(2,*),mres_a(2,*),resg_a(*)
    INTEGER :: ntot
    
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    ALLOCATE(grppt(2,ngrp),mres(2,nbun),resg(ngrp))
    grppt(:,1:ngrp)=grppt_a(:,1:ngrp)
    resg(1:ngrp)=resg_a(1:ngrp)
    mres(:,1:nbun)=mres_a(:,1:nbun)

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*
    
  END SUBROUTINE Initialize_Array
  SUBROUTINE Compute(xp0,yp0,zp0,co,type,sigma)
    IMPLICIT none
    
    REAL(8) :: xp0(*),yp0(*),zp0(*),co(3,3),sigma(*)
    INTEGER :: type(*)

!!$----------------------------- VARIABLES ------------------------------*

    INTEGER, SAVE :: time_of_call=0
    REAL(8), SAVE :: fact
    REAL(8) :: sqcut,sig_i,sig_j
    INTEGER :: m ,j,i,jj,ia,ja,st_nn,p_nn,count,ii,iii
    REAL(8) :: xpi,ypi,zpi,xpj,ypj,zpj,xg,yg,zg,xc,yc,zc,rsq
    INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: index_x

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*
    
    IF(time_of_call == 0) THEN
       fact=coeff*(1.0D0/(2.0D0**(1.0D0/6.0D0)))
       ALLOCATE(index_x(nbun))
       NULLIFY(neigh_s)
    END IF
    CALL Neigh_Delete(neigh_s)
    CALL Neigh_Start(neigh_s,index_st(1))

    time_of_call= time_of_call+1

    p_nn=0
    st_nn=0
    DO iii=nstart,nend
       ii=index_st(iii+1)
       index_x=0
       st_nn=st_nn+1
       DO i=mres(1,ii),mres(2,ii)
          p_nn=p_nn+1
          m=neigha(p_nn) % no
          DO jj=1,m
             j=neigha(p_nn) % nb(jj)
             DO ia=grppt(1,i),grppt(2,i)
                xpi=xp0(ia)
                ypi=yp0(ia)
                zpi=zp0(ia)
                sig_i=sigma(type(ia))
                DO ja=grppt(1,j),grppt(2,j)
                   sig_j=sigma(type(ja))
                   sqcut=(sig_j+sig_i)*fact
                   sqcut=sqcut**2
                   xpj=xp0(ja)
                   ypj=yp0(ja)
                   zpj=zp0(ja)
                   xg=xpi-xpj
                   yg=ypi-ypj
                   zg=zpi-zpj
                   xg=xg-2.0D0*PBC(xg)
                   yg=yg-2.0D0*PBC(yg)
                   zg=zg-2.0D0*PBC(zg)
                   xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
                   yc=           co(2,2)*yg+co(2,3)*zg
                   zc=                      co(3,3)*zg
                   rsq=xc*xc+yc*yc+zc*zc
                   IF(rsq < sqcut) THEN
                      index_x(resg(j))=1
                   END IF
                END DO
             END DO
          END DO
       END DO
       count=SUM(index_x)
       neigh_s(iii) % no = count
       IF(count /= 0) THEN
          ALLOCATE(neigh_s(iii) % nb(count))
          m=index_sv(1)
          count=0
          DO ii=1,m
             i=index_sv(1+ii)
             IF(index_x(i) == 1) THEN
                count=count+1
                neigh_s(iii) % nb(count)=i
             END IF
          END DO
       END IF
    END DO
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

  END SUBROUTINE Compute
  SUBROUTINE Compute_Neighbors(xpg,ypg,zpg,co)
    IMPLICIT NONE 

!!$----------------------------- ARGUMENTS ------------------------------*
    
    REAL(8) :: xpg(*),ypg(*),zpg(*),co(3,3)

!!$----------------------------- VARIABLES ------------------------------*

    INTEGER, SAVE :: time_of_call=0
    INTEGER :: i,m,nn,na,n,count
    
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    IF(time_of_call == 0) THEN
       indxyz=(2*(ncx-1)+1)*(2*(ncy-1)+1)*(2*(ncz-1)+1)
       ALLOCATE(indxi(indxyz),indxj(indxyz),indxk(indxyz))
       ALLOCATE(tags_sv(nbun))
       tags_sv=.TRUE.
       m=index_sv(1)
       DO i=1,m
          tags_sv(index_sv(i+1))=.FALSE.
       END DO
       count=0
       DO nn=nstart,nend
          na=index_st(nn+1)
          DO n=mres(1,na),mres(2,na) 
             count=count+1
          END DO
       END DO
       nlocal_neigh=count
    END IF
    time_of_call=time_of_call+1

    CALL lc_index(indxyz,ncx,ncy,ncz,nind,indxi,indxj,indxk&
         &,cutoff_max,co)
    CALL Get_neigh
  CONTAINS
    SUBROUTINE Get_Neigh
      IMPLICIT NONE 
      INTEGER ::  iret,indxyz
      INTEGER ::  nx,ny,nz
      INTEGER ::  i,j,k,l,m,n,nv,nvtot,numcell
      INTEGER ::  iv,jv,kv,nmin
      INTEGER ::  nppp,map,count,nn,na

      REAL(8) :: dx,dy,dz,co(3,3),rcut
      REAL(8) :: x1,y1,z1,x2,y2,z2,xx,yy,zz 
      REAL(8) :: sqcut,d

      INTEGER, DIMENSION (:), ALLOCATABLE :: headp,chainp,cellpi,cellpk,cellpj
      INTEGER, DIMENSION(:), ALLOCATABLE :: ind_a


      ALLOCATE(headp(ncx*ncy*ncz))
      ALLOCATE(chainp(ngrp),cellpi(ngrp),cellpj(ngrp),cellpk(ngrp))
      ALLOCATE(ind_a(ngrp))

      CALL Neigh_Delete(neigha)
      CALL Neigh_Start(neigha,nlocal_neigh)
      sqcut = cutoff_max**2
      dx=2.d0/ncx
      dy=2.d0/ncy
      dz=2.d0/ncz


      DO n=1,ncx*ncy*ncz
         headp(n)=0
      END DO

!!$*=======================================================================
!!$*     Compute chain list for system
!!$*=======================================================================

      DO n=1,ngrp
         x1=xpg(n)/dx
         y1=ypg(n)/dy
         z1=zpg(n)/dz
         nx=int(x1)+(sign(1.d0,x1-int(x1))-1.)/2
         ny=int(y1)+(sign(1.d0,y1-int(y1))-1.)/2
         nz=int(z1)+(sign(1.d0,z1-int(z1))-1.)/2
         nx=mod(mod(nx,ncx)+ncx,ncx)
         ny=mod(mod(ny,ncy)+ncy,ncy)
         nz=mod(mod(nz,ncz)+ncz,ncz)
         cellpi(n)=nx
         cellpj(n)=ny
         cellpk(n)=nz
         numcell=nz+ncz*(ny+ncy*nx)+1
         chainp(n)=headp(numcell)
         headp(numcell)=n
      END DO

!!$*=======================================================================
!!$*     Compute neighbor list nnlpp 
!!$*=======================================================================
      
      nppp=0
      nvtot=0
      count=0
      DO nn=nstart,nend
         na=index_st(nn+1)
         DO n=mres(1,na),mres(2,na) 
            count=count+1
            x1=xpg(n)
            y1=ypg(n)
            z1=zpg(n)
            nv=0         
            i=cellpi(n)
            j=cellpj(n)
            k=cellpk(n)
            DO m=1,nind
               iv=indxi(m)
               jv=indxj(m)
               kv=indxk(m)
               nx=mod(mod(i+iv,ncx)+ncx,ncx)
               ny=mod(mod(j+jv,ncy)+ncy,ncy)
               nz=mod(mod(k+kv,ncz)+ncz,ncz)
               numcell=nz+ncz*(ny+ncy*nx)+1
               l=headp(numcell)
               nmin=0
               IF(m == 1) nmin=n
               DO WHILE(l > nmin)
                  IF(l > n .AND. (.NOT. tags_sv(resg(l)))) THEN
!!$                  IF(.NOT. tags_sv(resg(l))) THEN
                     x2=x1-xpg(l)
                     y2=y1-ypg(l)
                     z2=z1-zpg(l)
                     x2=x2-2.0*pbc(x2)
                     y2=y2-2.0*pbc(y2)
                     z2=z2-2.0*pbc(z2)
                     xx=co(1,1)*x2+co(1,2)*y2+co(1,3)*z2
                     yy=           co(2,2)*y2+co(2,3)*z2
                     zz=                      co(3,3)*z2
                     d=xx**2+yy**2+zz**2
                     if(d < sqcut) then
                        nv=nv+1
                        ind_a(nv)=l
                     endif
                  endif
                  l=chainp(l)
               END DO
            END DO

            IF(nv /= 0) THEN
               ALLOCATE(neigha(count) % nb(nv))
               neigha(count) % no = nv
               neigha(count) % nb =ind_a(1:nv)
            END IF
            nvtot=nvtot+nv+1
         END DO
      END DO
#if defined PARALLEL
      IF(nprocs .GT. 1) CALL P_merge_i(nvtot)
#endif

      DEALLOCATE(headp)
      DEALLOCATE(chainp,cellpi,cellpj,cellpk,ind_a)
    END SUBROUTINE Get_Neigh
  END SUBROUTINE Compute_Neighbors
  REAL(8) FUNCTION PBC(x)
    REAL(8) :: x
    PBC=DNINT(0.5D0*x)
  END FUNCTION PBC
  SUBROUTINE Write_it(fstep)
    IMPLICIT NONE 

!!$----------------------------- ARGUMENTS ------------------------------*
    
    REAL(8) :: fstep
    
!!$----------------------------- VARIABLES ------------------------------*

    INTEGER :: ii,n,m
    TYPE(Neighbors), DIMENSION (:), POINTER :: neigh_sm

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    IF(node == 0) THEN
       IF(Kbinary == 0) THEN
          WRITE(khydration,'(''Solvation at '',f12.3,'' fs '')') fstep
       ELSE
          WRITE(khydration) fstep
       END IF
    END IF
    CALL Exchange
    IF(node == 0) THEN
       DO ii=1,SIZE(neigh_sm)
          m=neigh_sm(ii) % no
          IF(Kbinary == 0) THEN
             WRITE(khydration,'(2i6)') index_st(ii+1),m
             IF(m /=0) WRITE(khydration,'(12i6)') (neigh_sm(ii) % nb (n)&
                  &, n=1,m)
          ELSE
             WRITE(khydration) index_st(ii+1),m
             IF(m /=0) WRITE(khydration) (neigh_sm(ii) % nb (n)&
                  &, n=1,m)
          END IF
       END DO
    END IF
  CONTAINS
    SUBROUTINE Exchange
      IMPLICIT NONE 
#ifdef PARALLEL 
      INCLUDE 'mpif.h'
      include 'mpi_size.h'
#endif
      INTEGER, DIMENSION (:), ALLOCATABLE :: o,sizees,locales,displ
      INTEGER :: ii,i,n,m,size_tot,count,locale,nstarte,sizee,ierr


      size_tot=index_st(1)
      CALL Neigh_Delete(neigh_sm)
      CALL Neigh_Start(neigh_sm,size_tot)

#ifdef PARALLEL 

      ALLOCATE(locales(nprocs),displ(nprocs))
      count=0
      DO ii=nstart,nend
         m=neigh_s(ii) % no
         count=count+m+1
      END DO
      locale=count

      CALL MPI_ALLGATHER(locale,1,MPI_INTEGER4,locales,1 &
           ,MPI_INTEGER4,MPI_COMM_WORLD,ierr)

      displ(1)=0
      DO i=2,nprocs
         displ(i)=displ(i-1)+locales(i-1)
      END DO
      count=SUM(locales)

      ALLOCATE(o(count))
      count=displ(node+1)
      nstarte=count+1

      DO ii=nstart,nend
         m=neigh_s(ii) % no
         count=count+1
         o(count)=m
         DO n=1,m
            o(count+n)=neigh_s(ii) % nb(n)
         END DO
         count=count+m
      END DO

      CALL MPI_ALLGATHERV(o(nstarte),locale,MPI_INTEGER4,o,locales&
           &,displ,MPI_INTEGER4,MPI_COMM_WORLD,ierr)

      count=0
      DO ii=1,index_st(1)
         count=count+1
         m=o(count)
         neigh_sm(ii) % no = m
         ALLOCATE(neigh_sm(ii) % nb(m))
         DO n=1,m
            count=count+1
            neigh_sm(ii) % nb(n) = o(count)
         END DO
      END DO
      DEALLOCATE(o)
#else
      DO ii=1,SIZE(neigh_s)
         m=neigh_s(ii) % no
         neigh_sm(ii) % no =m
         ALLOCATE(neigh_sm(ii) % nb (m))
         DO n=1,m
            neigh_sm(ii) % nb (n)=neigh_s(ii) % nb(n)
         END DO
      END DO
#endif 
    END SUBROUTINE Exchange

  END SUBROUTINE Write_it
  
END MODULE HYDRATION_Mod
