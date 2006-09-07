MODULE VORONOI_Mod

!!$***********************************************************************
!!$   Time-stamp: <2006-08-11 15:39:11 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Wed Apr  5 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This MODULE is part of the program ORAC ----*

#include "config.h"

  USE Xerror_Mod

  USE INPUT_Mod, ONLY: Read_String, Parser, err_open,err_end,err_unr&
       &,err_fnf,err_args

  INTEGER :: max_neigh,pig_nnl,pnnlpp_vor,maxpla,maxver
  PARAMETER ( &
       & max_neigh=_MAX_NEIGH_, &
       & pnnlpp_vor=MAXPLA+1, &
       & maxpla=MAXPLA, &
       & maxver=MAXVER)

!!$======================== Input Parameters=============================*

  LOGICAL, SAVE :: voronoi=.FALSE.,heavy=.FALSE.,access=.FALSE.&
       &,volume=.FALSE.,fluct=.FALSE.,neighbor=.FALSE.,rewind_vor&
       &=.FALSE.,noprint=.FALSE.,compress=.FALSE.,dynamics=.FALSE.
  INTEGER, SAVE :: nvoronoi=1,kvoronoi=0,nfluct=2,bmax,counter=0&
       &,nprint_press=1,ncx=-1,ncy=-1,ncz=-1,ntype=3,kdynamics=0
  REAL(8), SAVE :: cutoff,bin=1.00D0,vol_tot=0.0D0,unit_fluct
  CHARACTER(80), SAVE :: filename
  TYPE RESIDUE
     INTEGER :: n0
     CHARACTER(8), DIMENSION(:), POINTER :: label
  END TYPE RESIDUE
  INTEGER, DIMENSION (:), ALLOCATABLE :: Res_Class
  TYPE(RESIDUE), DIMENSION(:), ALLOCATABLE, SAVE :: Res_type

!!$======================== Input Parameters end =========================*

!!$============================= Arrays ==================================*

  INTEGER, SAVE :: ntap,nbun,nunits,node
  INTEGER, DIMENSION (:,:), ALLOCATABLE, SAVE :: nnlpp_vor
  INTEGER, DIMENSION (:,:), ALLOCATABLE, SAVE :: res
  REAL(8), DIMENSION (:,:), ALLOCATABLE, SAVE :: area_vor
  REAL(8), DIMENSION (:), ALLOCATABLE, SAVE :: volume_vor,dist_min&
       &,sigmas,N1,N2,VOL1
  REAL(8), DIMENSION (:,:), ALLOCATABLE, SAVE :: MM1,MM2,VVOL1
  REAL(8), SAVE :: vol_slv
  REAL(8), ALLOCATABLE, SAVE :: vol_res(:),vol_slt(:),area_slt(:)&
       &,area_slv(:),area_frac(:),area_tot(:),vol_type(:)&
       &,surface_res(:),ratio_res(:)
  REAL(4), DIMENSION (:,:), ALLOCATABLE, SAVE :: dynas
  INTEGER, ALLOCATABLE, SAVE :: index(:),vol_typep(:),N_N1(:),M_MM1(:&
       &,:),dist_id(:)
  LOGICAL, DIMENSION (:), ALLOCATABLE, SAVE :: solvent

  
!!$============================= Arrays End ==============================*

!!$============================= Fluctations =============================*

  REAL(8), DIMENSION(:,:), ALLOCATABLE, SAVE :: dvij
  REAL(8), DIMENSION(:,:), ALLOCATABLE, SAVE :: dist

  REAL(8), DIMENSION(:), ALLOCATABLE, SAVE :: dvi
  INTEGER, SAVE :: fluct_calls=0
  
!!$============================= Fluctations End =========================*


CONTAINS
  SUBROUTINE Init(nodea,mma,natom,nbuna,mres,grppt,ss_index,prnb&
       &,nbtype,rcut,prsymb,mend,mtraj)

!!$======================== DECLARATIONS ================================*

    IMPLICIT none
    
    INTEGER :: nodea,mma,natom,nbuna,mtraj,mres(2,*), grppt(2,*)&
         &,ss_index(*),nbtype(*),mend(*)
    REAL(8) :: prnb(*),rcut
    CHARACTER(8) :: prsymb(*)

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INCLUDE 'unit.h'
    INTEGER :: n,m,i,count,nslv,ierr,ntraj
    LOGICAL :: ok
    REAL(8) :: facta,t=300.0D0,kt=300.0D0*8.3144/1000.0D0

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    ntraj=mtraj/nvoronoi
    unit_fluct=efact/unitp/kt/1000.0D0
    WRITE(*,*) 'Unit_Fluct = ',unit_fluct
    bmax=INT(rcut/bin)+1
    ALLOCATE(N1(-bmax:bmax),N2(-bmax:bmax),N_N1(-bmax:bmax),VOL1(&
         &-bmax:bmax)) 
    ALLOCATE(MM1(-bmax:bmax,ntype),MM2(-bmax:bmax,ntype),M_MM1(-bmax:bmax&
         &,ntype),VVOL1(-bmax:bmax,ntype))
    N1=0.0D0
    N2=0.0D0
    N_N1=0
    MM1=0.0D0
    MM2=0.0D0
    M_MM1=0
    VOL1=0.0D0
    VVOL1=0.0D0

    nbun=nbuna
    ntap=natom
    node=nodea
    ALLOCATE(solvent(ntap),dist_min(ntap),dist_id(ntap),sigmas(ntap))

    solvent=.FALSE.

    facta=1.0D0/(2.0D0**(1.0D0/6.0D0))
    DO n=1,ntap
       IF(ss_index(n) == 2) solvent(n)=.TRUE.
       sigmas(n)=prnb(nbtype(n))*facta
    END DO

    ALLOCATE(nnlpp_vor(pnnlpp_vor,mma))
    ALLOCATE(area_vor(pnnlpp_vor,mma),volume_vor(ntap)) 
    count=0
    DO n=1,nbun
       ok=.FALSE.
       DO m=mres(1,n),mres(2,n)
          DO i=grppt(1,m),grppt(2,m)
             IF(ss_index(i) == 1) ok=.TRUE.               
          END DO
       END DO
       IF(ok) THEN
          count=count+1
       END IF
    END DO
    ALLOCATE(res(2,nbun))
    ALLOCATE(vol_res(nbun),vol_slt(nbun),area_slt(nbun)&
         &,area_slv(nbun),area_frac(nbun),area_tot(nbun)&
         &,vol_type(nbun),surface_res(nbun),ratio_res(nbun)&
         &,index(nbun),vol_typep(nbun))
    ALLOCATE(dvij(count+1,count+1),dvi(count+1),dist(count+1,count+1))
    dvij=0.0D0
    dvi=0.0D0
    surface_res=0.0D0
    ratio_res=0.0D0
    nunits=count
    count=0
    DO n=1,nbun
       res(1,n)=grppt(1,mres(1,n))
       res(2,n)=grppt(2,mres(2,n))       
    END DO

    IF(fluct) THEN
       access=.TRUE.
       volume=.TRUE.
    END IF
    CALL Define_Residues
    CALL Assign_Residues
    ALLOCATE(dynas(nbun-nunits,ntraj))
    dynas=-1.0D0
  CONTAINS
    SUBROUTINE Define_Residues
      IMPLICIT NONE 
      ALLOCATE(Res_Class(ntap))
      Res_Class=-1

      ALLOCATE(Res_Type(5))

!!$      
!!$ ----- Aliphatic
!!$

      Res_Type(1)%n0 = 1
      ALLOCATE(Res_Type(1)%label(4))
      Res_Type(1)%label(1)='ala'
      Res_Type(1)%label(2)='val'
      Res_Type(1)%label(3)='leu'
      Res_Type(1)%label(4)='ile'

!!$      
!!$ ----- Non Polar
!!$

      Res_Type(2)%n0 = 1
      ALLOCATE(Res_Type(2)%label(4))
      Res_Type(2)%label(1)='gly'
      Res_Type(2)%label(2)='pro'
      Res_Type(2)%label(3)='cys'
      Res_Type(2)%label(4)='met'

!!$      
!!$ ----- Aromatic
!!$

      Res_Type(3)%n0 = 1
      ALLOCATE(Res_Type(3)%label(5))
      Res_Type(3)%label(1)='his'
      Res_Type(3)%label(2)='hsd'
      Res_Type(3)%label(3)='phe'
      Res_Type(3)%label(4)='tyr'
      Res_Type(3)%label(5)='trp'
      
!!$      
!!$ ----- Polar
!!$

      Res_Type(4)%n0 = 2
      ALLOCATE(Res_Type(4)%label(5))
      Res_Type(4)%label(1)='asn'
      Res_Type(4)%label(2)='gln'
      Res_Type(4)%label(3)='ser'
      Res_Type(4)%label(4)='thr'
      Res_Type(4)%label(5)='ace'

!!$      
!!$ ----- Charged
!!$

      Res_Type(5)%n0 = 3
      ALLOCATE(Res_Type(5)%label(4))
      Res_Type(5)%label(1)='lys'
      Res_Type(5)%label(2)='arg'
      Res_Type(5)%label(3)='asp'
      Res_Type(5)%label(4)='glu'

    END SUBROUTINE Define_Residues
    SUBROUTINE Assign_Residues
      IMPLICIT NONE 
      INTEGER :: n,m,o,r,p
      CHARACTER(8) :: label

      DO m=1,nbun
         label=prsymb(mend(m)) 
         DO n=1,5
            r=SIZE(Res_Type(n)%label)
            ok=.FALSE.
            DO o=1,r
               IF(Res_Type(n)%label(o)(1:3) == label(1:3)) THEN
                  ok=.TRUE.
                  DO p=res(1,m),res(2,m)
                     Res_Class(p)=Res_Type(n)%n0
                  END DO
                  EXIT
               END IF
            END DO
            IF(ok) EXIT
         END DO

      END DO
      
!!$      DO m=1,nbun
!!$         label=prsymb(mend(m))
!!$         DO p=res(1,m),res(2,m)
!!$            WRITE(*,'(i6,i3,a5)') m,Res_Class(p),label
!!$         END DO
!!$      END DO
    END SUBROUTINE Assign_Residues
  END SUBROUTINE Init

  SUBROUTINE Collect_Dynamics

!!$----------------------------- ARGUMENTS ------------------------------*

    IMPLICIT NONE 

    REAL(8) :: rmin
    INTEGER, SAVE :: counter=0
    INTEGER :: n,i

!!$------------------------- LOCAL VARIABLES ----------------------------*


    
    counter=counter+1
    IF(counter > SIZE(dynas,2)) RETURN
    DO n=nunits+1,nbun
       rmin=1.0D9
       DO i=res(1,n),res(2,n)
          IF(dist_min(i) < rmin) THEN
             rmin=dist_min(i)
          END IF
       END DO
       dynas(n-nunits,counter)=rmin
    END DO
    WRITE(*,*) SIZE(dynas,2),counter
    IF(SIZE(dynas,2) == counter) THEN
       REWIND(kdynamics)

       WRITE(*,*) '-----> Write Dynamic file at the end of run. Step =',counter

       WRITE(kdynamics) nbun-nunits,counter
       DO n=1,nbun-nunits
          WRITE(kdynamics) dynas(n,:)
       END DO
    ELSE IF(MOD(counter,500) == 0) THEN
       REWIND(kdynamics)

       WRITE(*,*) '-----> Write Dynamic file. Step =',counter

       WRITE(kdynamics) nbun-nunits,counter
       DO n=1,nbun-nunits
          WRITE(kdynamics) dynas(n,1:counter)
       END DO
    END IF

  END SUBROUTINE Collect_Dynamics
  SUBROUTINE density(volume)

!!$----------------------------- ARGUMENTS ------------------------------*

    IMPLICIT NONE 

    REAL(8) :: volume

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: i,n,ih,o,i_min
    REAL(8) :: rmin
    REAL(8), DIMENSION (:), ALLOCATABLE :: NW,N_NW
    REAL(8), DIMENSION (:,:), ALLOCATABLE :: MMW,M_MMW

    ALLOCATE(NW(-bmax:bmax),N_NW(-bmax:bmax),MMW(-bmax:bmax,ntype),M_MMW(&
         &-bmax:bmax,ntype))

    counter=counter+1
    
    NW=0.0D0
    N_NW=0.0D0
    MMW=0.0D0
    M_MMW=0.0D0
    Vol_tot=Vol_Tot+volume
    WRITE(*,'(''Cell Volume ='',f14.5,'' Dev. = '',f14.5,'' Ang. '',e1&
         &4.6,''  % '')') volume,volume-SUM(vol_res),100.0D0*(volume&
         &-SUM(vol_res))/volume 
    DO n=nunits+1,nbun
       rmin=1.0D9
       DO i=res(1,n),res(2,n)
          IF(dist_min(i) < rmin) THEN
             rmin=dist_min(i)
             i_min=dist_id(i)
          END IF
       END DO
       IF(rmin < 1.0D9) THEN
          ih=INT(rmin/bin)+1
          NW(ih)=NW(ih)+vol_res(n)
          N_N1(ih)=N_N1(ih)+1
          N_NW(ih)=N_NW(ih)+1.0D0
          
          o=Res_Class(i_min)
          IF(o > 0) THEN
             MMW(ih,o)=MMW(ih,o)+vol_res(n)
             M_MMW(ih,o)=M_MMW(ih,o)+1.0D0
             M_MM1(ih,o)=M_MM1(ih,o)+1
          END IF
       END IF
    END DO

    DO ih=-bmax,bmax
       IF(N_NW(ih) /= 0.0D0) THEN
          N1(ih)=N1(ih)+NW(ih)/N_NW(ih)
          N2(ih)=N2(ih)+volume*NW(ih)/N_NW(ih)
       END IF
    END DO

    DO o=1,ntype
       DO ih=-bmax,bmax
          IF(M_MMW(ih,o) /= 0.0D0) THEN
             MM1(ih,o)=MM1(ih,o)+MMW(ih,o)/M_MMW(ih,o)
             MM2(ih,o)=MM2(ih,o)+volume*MMW(ih,o)/M_MMW(ih,o)
          END IF
       END DO
    END DO

    DEALLOCATE(NW,MMW,N_NW,M_MMW)
  END SUBROUTINE Density

  SUBROUTINE print_density

!!$----------------------------- ARGUMENTS ------------------------------*

    IMPLICIT NONE 

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: i,n,ih,o
    REAL(8) :: n00,n0,nn,vol,kk,r,aux,auxx,dens

    IF(MOD(counter,nprint_press) == 0) THEN
       DO ih=-bmax,bmax
          IF(N_N1(ih) /= 0) THEN
             r=DBLE(ih-1)*bin
             aux=DBLE(N_N1(ih))/DBLE(counter)

             auxx=VOL1(ih)/DBLE(counter)
             dens=aux/auxx

             kk=0.0D0

             n00=N1(ih)/DBLE(counter)

             n0=n00*aux
             nn=aux*N2(ih)/DBLE(counter)
             vol=vol_tot/DBLE(counter)

             IF(n0 /= 0.0D0) kk=(nn-n0*vol)/n0
             WRITE(*,'(f12.3,5e14.6)') r,n00,aux,dens,kk*unit_fluct
          END IF
       END DO
       WRITE(*,*) '#  '
       WRITE(*,*) '#  '

       DO o=1,ntype
          WRITE(*,*) '# Class ',o
          DO ih=-bmax,bmax
             IF(M_MM1(ih,o) /= 0) THEN
                r=DBLE(ih-1)*bin
                aux=DBLE(M_MM1(ih,o))/DBLE(counter)

                auxx=VVOL1(ih,o)/DBLE(counter)
                dens=aux/auxx


                kk=0.0D0
                n00=MM1(ih,o)/DBLE(counter)
                n0=n00*aux
                nn=aux*MM2(ih,o)/DBLE(counter)
                vol=vol_tot/DBLE(counter)
                IF(n0 /= 0.0D0) kk=(nn-n0*vol)/n0
                WRITE(*,'(f12.3,4e14.6)') r,n00,aux,dens,kk*unit_fluct
             END IF
          END DO
       END DO
    END IF
    
  END SUBROUTINE Print_density
  SUBROUTINE Volume_Shell(volume,co,xp0,yp0,zp0)
    IMPLICIT NONE
    REAL(8) :: volume,co(3,3),xp0(*),yp0(*),zp0(*)



    REAL(8), SAVE :: one=1.0D0,ncxi,ncyi,nczi,Dvolume
    INTEGER :: i,ii,k,l,m,i_min,o,ih
    INTEGER, SAVE :: times=0
    REAL(8) :: xb,xd,xc,yb,yd,yc,zb,zd,zc,rmin,rs,rsp
    LOGICAL :: ok

    IF(ncx < 1) RETURN

    Dvolume=volume/DBLE(ncx*ncy*ncz)

    IF(Times == 0) THEN
       Times=1
       WRITE(*,*) 'Dvolume =',Dvolume
    END IF

    ncxi=2.0D0/DBLE(ncx)
    ncyi=2.0D0/DBLE(ncy)
    nczi=2.0D0/DBLE(ncz)

    DO k=1,ncx
       xb=ncxi*DBLE(k)-one
       DO l=1,ncy
          yb=ncyi*DBLE(l)-one
          DO m=1,ncz
             zb=nczi*DBLE(m)-one
             rmin=1.0D9
             i_min=-100
             ok=.FALSE.
             DO ii=1,nunits
                DO i=res(1,ii),res(2,ii)
                   xd=xb-xp0(i)
                   yd=yb-yp0(i)
                   zd=zb-zp0(i)
                   xd=xd-2.0D0*PBC(xd)
                   yd=yd-2.0D0*PBC(yd)
                   zd=zd-2.0D0*PBC(zd)
                   xc=co(1,1)*xd+co(1,2)*yd+co(1,3)*zd
                   yc=co(2,1)*xd+co(2,2)*yd+co(2,3)*zd
                   zc=co(3,1)*xd+co(3,2)*yd+co(3,3)*zd
                   rs=xc*xc+yc*yc+zc*zc
                   rsp=SQRT(rs)-sigmas(i)
                   IF(rmin .GT. rsp) THEN
                      rmin=rsp
                      i_min=i
                      ok=.TRUE.
                   END IF
                END DO
             END DO
             IF(ok) THEN
                ih=INT(rmin/bin)+1
                VOL1(ih)=VOL1(ih)+Dvolume
                o=Res_Class(i_min)
                VVOL1(ih,o)=VVOL1(ih,o)+Dvolume
             END IF
          END DO
       END DO
    END DO

  END SUBROUTINE Volume_Shell

  SUBROUTINE Fluctuations(kprint,xp0,yp0,zp0)

!!$----------------------------- ARGUMENTS ------------------------------*

    IMPLICIT NONE 

    REAL(8) :: xp0(*),yp0(*),zp0(*)
    INTEGER :: kprint

!!$----------------------- VARIABLES IN COMMON --------------------------*

!!$------------------------- LOCAL VARIABLES ----------------------------*

    REAL(8), DIMENSION(:,:), ALLOCATABLE, SAVE :: Kij,Kij_orig
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: U,VT,sigma,Multiply
    REAL(8), DIMENSION(:), ALLOCATABLE, SAVE :: vols,work,sing
    
    INTEGER(8) :: lda,info,lwork
    INTEGER :: n_U,i,j
    REAL(8) :: f_Calls,xpi,xpj

    IF(rewind_vor) REWIND(kvoronoi)

    n_U=nunits+1
    
    ALLOCATE(vols(n_U))
         

    vols(1:nunits)=vol_res(1:nunits)
    vols(n_U)=vol_slv

    fluct_calls=fluct_calls+1
    f_Calls=1.0D0/DBLE(fluct_calls)

    DO i=1,n_U
       xpi=vols(i)
       dvi(i)=dvi(i)+xpi
       DO j=i,n_U
          xpj=vols(j)
          dvij(i,j)=dvij(i,j)+xpi*xpj
       END DO
    END DO
    DO i=1,n_U-1
       surface_res(i)=surface_res(i)+area_slv(i)
       ratio_res(i)=ratio_res(i)+area_frac(i)
    END DO

    CALL Get_Dist

    DEALLOCATE(vols)

    IF(MOD(fluct_calls,nfluct) == 0) THEN
       ALLOCATE(U(n_U,n_U),VT(n_U,n_U),sigma(n_U,n_U),sing(n_U)&
            &,kij_orig(n_U,n_U),kij(n_U,n_U),Multiply(n_U,n_U))
       U=0.0D0
       VT=0.0D0
       sigma=0.0D0
       sing=0.0D0
       Multiply=0.0D0
       Kij=0.0D0
       DO i=1,n_U
          xpi=dvi(i)*f_Calls
          DO j=i,n_U
             xpj=dvi(j)*f_Calls
             Kij(i,j)=dvij(i,j)*f_Calls-xpi*xpj
             Kij_orig(i,j)=Kij(i,j)

             Kij(j,i)=Kij(i,j)
             Kij_orig(j,i)=Kij_orig(i,j)
          END DO
       END DO

       IF(node ==0) THEN
          WRITE(kvoronoi,'(''Average Volumes '',i7)') fluct_calls
          DO i=1,n_U
             IF(i /= n_U) THEN
                WRITE(kvoronoi,'(i7,2f15.4,f15.5)') i,dvi(i)*f_Calls&
                     &,surface_res(i)*f_Calls,ratio_res(i)*f_Calls
             ELSE
                WRITE(kvoronoi,'(i7,2f15.5)') i,dvi(i)*f_Calls
             END IF
          END DO

          WRITE(kvoronoi,'(''Direct Fluctuations at step '',i7)') fluct_calls
          
          DO i=1,n_U
             DO j=1,n_U
                WRITE(kvoronoi,'(2i7,e15.7)') i,j,Kij(i,j)             
             END DO
          END DO
       END IF

       lda=SIZE(Kij,1)

       lwork=10

       ALLOCATE(work(lwork))

       lwork=-1

       CALL DGESVD('A','A',lda,lda,Kij,lda,sing,U,lda,vt,lda,work,lwork,info)
       lwork=INT(work(1))
       DEALLOCATE(work)

       ALLOCATE(work(lwork))
       CALL DGESVD('A','A',lda,lda,Kij,lda,sing,U,lda,vt,lda,work,lwork,info)
       DEALLOCATE(work)

       sigma=0.0D0
       DO i=1,n_U
          IF(DABS(sing(i)) > 1.0D-6) THEN
             sigma(i,i)=1.0D0/sing(i)
          END IF
       END DO

       vt=TRANSPOSE(vt)
       U=TRANSPOSE(U)
       U=MATMUL(sigma,U)
       Kij=MATMUL(vt,U)

       IF(node ==0) THEN
          WRITE(kvoronoi,'(''Inverse Fluctuations at step '',i7)') fluct_calls
          
          DO i=1,n_U
             DO j=i,n_U
                WRITE(kvoronoi,'(2i7,e15.7,2x,f15.7)') i,j,Kij(i,j)&
                     &,dist(i,j)*f_Calls
             END DO
          END DO
       END IF
       DEALLOCATE(U,VT,sigma,sing,kij_orig,Kij,Multiply)
    END IF

  CONTAINS
    SUBROUTINE Get_Dist
      IMPLICIT NONE 
      INTEGER :: n,m
      REAL(8), DIMENSION(:), ALLOCATABLE :: rx,ry,rz
      REAL(8) :: fcount,xpi,ypi,zpi,xc,yc,zc

      ALLOCATE(rx(nunits),ry(nunits),rz(nunits))
      
      DO n=1,nunits
         rx(n)=0.0D0
         ry(n)=0.0D0
         rz(n)=0.0D0
         DO m=res(1,n),res(2,n)
            rx(n)=rx(n)+xp0(m)
            ry(n)=ry(n)+yp0(m)
            rz(n)=rz(n)+zp0(m)
         END DO
         fcount=1.0D0/DBLE(res(2,n)-res(1,n)+1)
         rx(n)=rx(n)*fcount
         ry(n)=ry(n)*fcount
         rz(n)=rz(n)*fcount
      END DO
      DO n=1,nunits
         xpi=rx(n)
         ypi=ry(n)
         zpi=rz(n)
         DO m=n,nunits
            xc=xpi-rx(m)
            yc=ypi-ry(m)
            zc=zpi-rz(m)
            dist(n,m)=dist(n,m)+SQRT(xc**2+yc**2+zc**2)
         END DO
      END DO
      DEALLOCATE(rx,ry,rz)

    END SUBROUTINE Get_Dist
  END SUBROUTINE Fluctuations

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*
  REAL(8) FUNCTION PBC(x)
    REAL(8) :: x
    PBC=DNINT(0.5D0*x)
  END FUNCTION PBC
  INCLUDE 'VORONOI_Read.f90'
END MODULE VORONOI_Mod
