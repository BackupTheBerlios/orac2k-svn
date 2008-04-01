MODULE FIELD_Mod

!!$***********************************************************************
!!$   Time-stamp: <2008-03-14 12:01:39 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Mon Nov 19 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program orac2k ----*


!!$======================== DECLARATIONS ================================*

  USE DENSITY_Mod, ONLY: DEN_write=>n_write,Density_Calc&
       &,DEN_Initialize=>Initialize,DEN_Initialize_Ar=>Initialize_Ar&
       &,DEN_compute=>Compute,DEN_write_it=>write_it&
       &,DEN_Initialize_Par=>Initialize_Par

  USE Node
  USE INPUT_Mod, ONLY: Read_String, Parser, err_open,err_end,err_unr&
       &,err_fnf,err_args
  USE PDB, ONLY: PDB_out=>Write_it, PDB_
  USE REDUCE, ONLY: REDUCE_Phi, REDUCE_Forces, REDUCE_Init=>Init
  IMPLICIT none
  PRIVATE
  PUBLIC :: Compute, Init,Read_it, Electric, n_write

  REAL(8), SAVE :: a=0.0D0,b=0.0D0,c=0.0D0,alph,bet,gamm,co(3,3),oc(3&
       &,3),volume,xcm,ycm,zcm,Dvolume,rho=0.0D0,unitfluc,UnitElec&
       &,dx_Pr=0.1D0,max_Pr=15.0D0
  INTEGER, SAVE :: natom,n_write=1,nccx=0,nccy=0,nccz=0,ntot=1,nats=0&
       &,natoms_Tot=0,nmax_Pr
  INTEGER, SAVE :: kelectric=0,kelectric2=0,kpdb=0,kpr=0

  INTEGER, SAVE :: nstart,nend,nstart_a,nend_a,nlocal_a,noeds=0,nprocs=1&
       &,ncube,nodex,nodey,nodez,nstart_h,nend_h,nlocal_h &
       &,nstart_ah,nend_ah,nlocal_ah,nstart_1,nend_1,nlocal_1&
       &,nstart_ex0,nend_ex0,nlocal_ex0,nstart_ex,nend_ex,nlocal_ex&
       &,nstart_2,nend_2,nlocal_2,fudgec,nbyte=4,ncpu_h

  CHARACTER(1) ::  P_shell
  INTEGER, ALLOCATABLE, SAVE :: mapnl(:),worka(:),tag_bndg(:)
  REAL(8), ALLOCATABLE, SAVE :: cpu_h(:)
  CHARACTER(80), SAVE :: filename_ele,filename_ele2,filename_pdb&
       &,filename_Pr
  CHARACTER(8), SAVE :: file_format='xplor   ',Target_Res='all'
  LOGICAL, SAVE :: Electric=.FALSE.,do_Pr=.FALSE.
  INTEGER, DIMENSION (:), ALLOCATABLE, SAVE :: atoms
  INTEGER, DIMENSION (:), ALLOCATABLE, SAVE :: indxi,indxj,indxk
  INTEGER, SAVE :: indxyz,nind
  INTEGER, ALLOCATABLE :: mapdn(:,:),nmapdn(:)
  REAL(8), ALLOCATABLE :: fpx_p(:),fpy_p(:),fpz_p(:),stressd_p(:,:)
  REAL(8), ALLOCATABLE :: fpx_h(:),fpy_h(:),fpz_h(:),stressd_h(:,:)&
       &,stressr_h(:,:),phi(:),phid(:)
  REAL(8), ALLOCATABLE, SAVE :: Pr(:)
  INTEGER, ALLOCATABLE, SAVE :: nPr(:)
CONTAINS
!!$
!!$-------- COMPUTE 
!!$
  SUBROUTINE Compute(xp0,yp0,zp0,xpa,ypa,zpa,xpga,ypga,zpga,xpcma&
       &,ypcma,zpcma,co,oc,volume,fstep)
    REAL(8) :: xpa(*),ypa(*),zpa(*),xpga(*),ypga(*),zpga(*),xpcma(*)&
         &,ypcma(*),zpcma(*),co(3,3),oc(3,3),xp0(*),yp0(*),zp0(*)&
         &,volume,fstep
    INTEGER, SAVE :: No_of_Calls=0

    No_of_Calls=No_of_Calls+1
    CALL Update
    CALL Direct
    IF(noeds == 0) THEN
       IF(do_Pr) CALL P_of_r
       CALL DEN_Compute(xp0,yp0,zp0,volume,phi)
       IF(MOD(No_of_Calls,DEN_write) == 0) THEN
          CALL DEN_write_it(fstep)
          IF(do_Pr) CALL Pr_Write(fstep)
       END IF
    END IF
  CONTAINS
    SUBROUTINE P_of_r
      USE PBC_Mod
      IMPLICIT NONE 
      REAL(8) :: xa,ya,za,xb,yb,zb,xi,yi,zi,xj,yj,zj,rs,bin
      INTEGER :: n,m
      LOGICAL :: exist

      IF(No_of_calls == 1) THEN
         n=INT(max_Pr/dx_Pr)+1
         nmax_Pr=n
         ALLOCATE(Pr(n),nPr(n))
         Pr=0.0D0
         nPr=0
         filename_Pr='ELECTRIC_Pr.dat'
         INQUIRE(FILE=filename_Pr,EXIST=exist)
         IF(exist) THEN
            CALL openf(kpr,filename_Pr,'FORMATTED','OLD',0)
         ELSE
            CALL openf(kpr,filename_Pr,'FORMATTED','NEW',0)
         END IF
         RETURN
      END IF

      xi=xpa(1)
      yi=ypa(1)
      zi=zpa(1)

      DO n=1,natom
         xj=xpa(n)
         yj=ypa(n)
         zj=zpa(n)
         xa=xj-xi
         ya=yj-yi
         za=zj-zi
         xa=xa-2.0D0*PBC(xa)
         ya=ya-2.0D0*PBC(ya)
         za=za-2.0D0*PBC(za)
         xb=co(1,1)*xa+co(1,2)*ya+co(1,3)*za
         yb=co(2,1)*xa+co(2,2)*ya+co(2,3)*za
         zb=co(3,1)*xa+co(3,2)*ya+co(3,3)*za
         rs=SQRT(xb*xb+yb*yb+zb*zb)
         bin=INT(rs/dx_Pr)
         IF(bin <= nmax_Pr) THEN
            Pr(bin)=Pr(bin)+phi(n)
            NPr(bin)=NPr(bin)+1
         END IF
      END DO
    END SUBROUTINE P_of_r
    SUBROUTINE Pr_Write(fstep)
      IMPLICIT NONE 
      REAL(8) :: fstep
      REAL(8) :: rs
      INTEGER :: n,m,nb
      REAL(8), POINTER :: Pn(:)

      WRITE(kpr,'(''P of r at '',f12.3,'' fs '')') fstep
      DO n=0,nmax_Pr
         rs=dx_Pr*n
         IF(nPr(n) /= 0) THEN
            WRITE(kpr,'(f12.3,2x,e15.8)') rs,Pr(n)/nPr(n)
         ELSE
            WRITE(kpr,'(f12.3,2x,e15.8)') rs,0.0D0
         END IF
      END DO
      CLOSE(kpr)
      CALL openf(kpr,filename_Pr,'FORMATTED','OLD',0)
      REWIND(kpr)
    END SUBROUTINE Pr_Write
    SUBROUTINE Direct
      IMPLICIT NONE 
      INCLUDE 'parst.h'
      INCLUDE 'cpropar.h'
      INCLUDE 'pme.h'
      INCLUDE 'parallel.h'
      INCLUDE 'unit.h'
      CHARACTER(1) :: rshell,rshk
      REAL(8) :: ucns_h,ucos_h,virs_h,virsp_h,ucnp_h,ucop_h,ucnsp_h&
           &,ucosp_h
      REAL(8) :: fsin14,fsbend,fsbond,fscnstr_slt,fscnstr_slv&
           &,coul_bnd_slt,coul_bnd_slv,eer_h,urcsp_h,urcs_h,urcp_h&
           &,virp_h
      REAL(8) ::   unb14_slt,cnb14_slt,ungrp_slt,cngrp_slt,uptors_slt&
           &,unb14_slv,cnb14_slv,ungrp_slv,cngrp_slv,uptors_slv
      INTEGER :: i
      

      rshell='h'
      rshk='h'
      fpx_h=0.0D0
      fpy_h=0.0D0
      fpz_h=0.0D0
      ucns_h=0.0D0
      ucos_h=0.0D0
      virs_h=0.0D0
      virsp_h=0.0D0
      ucnp_h=0.0D0
      ucop_h=0.0D0
      ucnsp_h=0.0D0
      ucosp_h=0.0D0
      phid(1:ntap)=0.0D0
      CALL mts_forces(rshell,xpa,ypa,zpa,xpga,ypga,zpga,xpcma,ypcma&
           &,zpcma,mapnl,mapdn,nmapdn,ucns_h,ucos_h,virs_h,virsp_h&
           &,ucnp_h,ucop_h,ucnsp_h,ucosp_h,fpx_h,fpy_h,fpz_h&
           &,phid,stressd_h,worka,cpu_h,ncpu_h,nstart_h,nend_h,nstart_ah&
           &,nend_ah,nlocal_ah,noeds,nprocs,ncube,P_dyn_update_shell)
#if defined PARALLEL
      CALL REDUCE_Phi(phid)
#endif      

      CALL fnb14(ss_index,xp0,yp0,zp0,chrge,ntap,ecc1412,ecc146&
           &,cutoff,clewld,alphal,int14,int14p,type14,int14_x,fudge&
           &,lj_fudge,fpx_h,fpy_h,fpz_h,phid,unb14_slt,unb14_slv&
           &,cnb14_slt,cnb14_slv)

      CALL fnbgrp(ss_index,xp0,yp0,zp0,chrge,nbtype,ecc12,ecc6,clewld&
           &,alphal,ingrp,ingrpp,ingrp_x,fpx_h,fpy_h,fpz_h,phid&
           &,ungrp_slt,ungrp_slv,cngrp_slt,cngrp_slv)

#if defined PARALLEL
      CALL P_expand_r8(Phid,nstart_1,nend_1,nlocal_1,noeds,nprocs)
#endif      

      phi(1:ntap)=0.0D0
      CALL mts_furier(noeds,nodex,nodey,nodez,ictxt,npy,npz,descQ&
           &,nprocs,ncube,nstart_1,nend_1,nlocal_1,nstart_2,nend_2&
           &,nlocal_2,xp0,yp0,zp0,xpa,ypa,zpa,xpcma,ypcma,zpcma&
           &,urcsp_h,urcs_h,urcp_h,virsp_h,virs_h,virp_h,fpx_h,fpy_h&
           &,fpz_h,phi,fsin14,fsbend,fsbond,fscnstr_slt,fscnstr_slv&
           &,coul_bnd_slt,coul_bnd_slv,rshell,rshk,eer_h,stressr_h&
           &,fudgec,tag_bndg) 
      
#if defined PARALLEL
      CALL REDUCE_Phi(phi)
      CALL P_expand_r8(Phi,nstart_2,nend_2,nlocal_2,noeds,nprocs)
#endif

      CALL cself_phi(1,ntap,alphal,chrge,volume,phi)
      phi=(phi+phid)*efact/1000.0D0/DSQRT(unitc)
#if defined PARALLEL
      CALL P_Barrier
#endif
    END SUBROUTINE Direct
    SUBROUTINE Update
      IMPLICIT NONE 
      INCLUDE 'parst.h'
      INCLUDE 'cpropar.h'
      INCLUDE 'pme.h'
      INCLUDE 'parallel.h'
      INCLUDE 'unit.h'
      REAL(8) :: ucns_p,ucos_p,virs_p,virsp_p,ucnp_p,ucop_p,ucnsp_p&
           &,ucosp_p,aux
      

      IF( .NOT. linked_cell) THEN
!!$*--      update shell h neighbor list
         CALL mts_forces('u',xpa,ypa,zpa,xpga,ypga,zpga,xpcma,ypcma&
              &,zpcma,mapnl,mapdn,nmapdn,ucns_p,ucos_p,virs_p&
              &,virsp_p,ucnp_p,ucop_p,ucnsp_p,ucosp_p,fpx_p,fpy_p&
              &,fpz_p,phi,stressd_p,worka,cpu_h,ncpu_h,nstart_h,nend_h&
              &,nstart_ah,nend_ah,nlocal_ah,node,nprocs,ncube&
              &,P_dyn_update_shell) 
      ELSE
         aux= rcuth+rtolh+rneih
         CALL lc_index(indxyz,ncx,ncy,ncz,nind,indxi,indxj,indxk,aux&
              &,co)
         CALL lc_list(ncx,ncy,ncz,nind,indxi,indxj,indxk,aux,co,xpga&
              &,ypga,zpga,ngrp,nstart_h,nend_h,node,nprocs,ncube&
              &,worka,kprint,.TRUE.) 
      END IF
    END SUBROUTINE Update
  END SUBROUTINE Compute
!!$
!!$-------- Initialize 
!!$
  SUBROUTINE Init(xpa,ypa,zpa,xpga,ypga,zpga,xpcma,ypcma,zpcma,co,oc&
       &,nodea,nprocsa,prsymb,beta,res1,res2,mass,natom0)

    IMPLICIT NONE
    CHARACTER(8) :: prsymb(*)
    CHARACTER(7) :: beta(*)
    INTEGER :: res1(*),res2(*),natom0
    REAL(8) :: mass(*)

    INTEGER :: nodea,nprocsa,ntapa,ngrpa
    REAL(8) :: xpa(*),ypa(*),zpa(*),xpga(*),ypga(*),zpga(*),xpcma(*)&
         &,ypcma(*),zpcma(*),co(3,3),oc(3,3)

    REAL(8) :: ucns_p,ucos_p,virs_p,virsp_p,ucnp_p,ucop_p,ucnsp_p&
         &,ucosp_p

    noeds=nodea; nprocs=nprocsa
    

    IF(noeds == 0) WRITE(*,*) '------------> No. of Nodes for the run  = '&
         &,nprocs,' <----------------'
    CALL Map
    ncpu_h=nprocs
    ALLOCATE(cpu_h(nprocs))

    CALL SetUp

    IF(natoms_Tot /=0) filename_pdb='ELECTRIC_FILE_p.pdb'
    SELECT CASE(file_format)
    CASE DEFAULT
       filename_ele='ELECTRIC_FILE_e.cube'
       filename_ele2='ELECTRIC2_FILE.cube'
    CASE('xplor')
       filename_ele='ELECTRIC_FILE_e.xplor'
       filename_ele2='ELECTRIC2_FILE.xplor'
    END SELECT

    IF(noeds == 0) THEN
       CALL DEN_Initialize_Par(a,b,c,nccx,nccy,nccz,n_write&
            &,file_format,atoms,nats,natoms_Tot,filename_ele&
            &,filename_ele2,filename_pdb,Target_Res)
       CALL DEN_Initialize(noeds,prsymb,beta,res1,res2,mass,natom0)
       CALL DEN_Initialize_Ar
    END IF
   CONTAINS
     SUBROUTINE Map
       IMPLICIT NONE 
       INCLUDE 'parst.h'
       INCLUDE 'cpropar.h'
       LOGICAL, POINTER :: lwork(:)
       INTEGER :: i,j,j1,j2,l,n,m,noff,count_a,count_out
       LOGICAL :: ok

       
       ALLOCATE(fpx_p(1),fpy_p(1),fpz_p(1),stressd_p(3,3),mapdn(2,1)&
            &,nmapdn(1))
       ALLOCATE(fpx_h(ntap),fpy_h(ntap),fpz_h(ntap),stressd_h(3,3)&
            &,stressr_h(3,3),phi(ntap),phid(ntap))

       ALLOCATE(lwork(ntap))

       n=0
       IF(.NOT. Node_()) STOP
       DO i=1,ntap
          lwork(1:ntap-i)=.TRUE.
          DO l=1,lbond
             ok=.FALSE.
             IF(lbnd(1,l) .EQ. i) THEN
                j1=lbnd(2,l)
                ok=.TRUE.
             END IF
             IF(lbnd(2,l) .EQ. i) THEN
                j1=lbnd(1,l)
                ok=.TRUE.
             END IF
             IF(ok) THEN
                IF(j1 .GT. i) THEN
                   lwork(j1-i)=.FALSE.
                END IF
             END IF
          END DO
          
          
          DO l=1,lbend
             ok=.FALSE.
             IF(lbndg(1,l).EQ.i) THEN
                j1=lbndg(2,l)
                j2=lbndg(3,l)
                ok=.TRUE.
             END IF
             IF(lbndg(2,l).EQ.i) THEN
                j1=lbndg(3,l)
                j2=lbndg(1,l)
                ok=.TRUE.
             END IF
             IF(lbndg(3,l).EQ.i) THEN
                j1=lbndg(1,l)
                j2=lbndg(2,l)
                ok=.TRUE.
             END IF
             IF(ok) THEN
                IF(j1.GT.i) THEN
                   lwork(j1-i)=.FALSE.
                END IF
                IF(j2.GT.i) THEN
                   lwork(j2-i)=.FALSE.
                END IF
             END IF
          END DO
          
          IF(int14p.NE.0) THEN
             DO l=1,int14p
                IF(int14(1,l).EQ.i) THEN
                   j1=int14(2,l)
                   IF(j1.GT.i) THEN
                      lwork(j1-i)=.FALSE.
                   END IF
                ELSE IF(int14(2,l).EQ.i) THEN
                   j1=int14(1,l)
                   IF(j1.GT.i) THEN
                      lwork(j1-i)=.FALSE.
                   END IF
                END IF
             END DO
          END IF
          m=0
          DO j=i+1,ntap
             IF(.NOT.lwork(j-i)) THEN
                m=m+1
             END IF
          END DO
          CALL Node__Push(m)
          m=0
          DO j=i+1,ntap
             IF(.NOT.lwork(j-i)) THEN
                m=m+1
                CALL Node__Push(j)
             END IF
          END DO
          noff=m+1
          n=n+noff
       END DO
       count_out=Node__Size()

       ALLOCATE(mapnl(count_out))
       count_a=0
       DO WHILE(Node__Pop(j))
          count_A=count_A+1
          mapnl(count_a)=j
       END DO
       CALL Node__Delete()

     END SUBROUTINE Map
     SUBROUTINE SetUp
       INCLUDE 'parst.h'
       INCLUDE 'cpropar.h'
       INCLUDE 'pme.h'
       INCLUDE 'parallel.h'

       INTEGER :: offset,abmd_dir,cnstpp,cnst_protp,cnst_protp_1&
            &,cnst_protp_ex,ffwork(2)
       INTEGER, POINTER :: cnstp(:,:), cnst_protl(:),cnst_protl_1(:)&
            &,cnst_protl_ex(:)
       INTEGER :: iret,i,count,i1,m,j
       CHARACTER(80) :: errmsg
       CHARACTER(1)  :: rshk
       REAL(8), POINTER ::  work(:)

!!$       chrge(2:ntap)=0.0D0
!!$       chrge(2:4)=0.0D0
       iret=0
       ALLOCATE(work(mspline))
       natom=ntap
       pme =.TRUE.                                          
       grpcut=.FALSE.
       clewld=.TRUE.                                           

       cnstpp=0
       cnst_protp=0
       
       ALLOCATE(worka(ntap))
#if !defined PARALLEL
       nstart_h=1
       nend_h=ngrp
       nlocal_h=ngrp
       nstart_ah=1
       nend_ah=ntap
       nlocal_ah=ntap
       
#else
       
!!$*----- Fake split and setup
       
       CALL P_split_scalar0(noeds,nprocs,nstart_h,nend_h,nlocal_h&
            &,nstart_ah,nend_ah,nlocal_ah,ngrp,grppt) 
!!$*----- Compute neighbors

       P_shell='l'
       CALL mts_forces('z',xpa,ypa,zpa,xpga,ypga,zpga,xpcma,ypcma&
            &,zpcma,mapnl,mapdn,nmapdn,ucns_p,ucos_p,virs_p,virsp_p&
            &,ucnp_p,ucop_p,ucnsp_p,ucosp_p,fpx_p,fpy_p,fpz_p&
            &,phi,stressd_p,worka,cpu_h,ncpu_h,nstart_h,nend_h,nstart_ah&
            &,nend_ah,nlocal_ah,noeds,nprocs,ncube,P_shell)
       
!!$*----- Now can split and set up correctly (nstart_h)

       cpu_h(1:nprocs)=0.0D0
       CALL P_split_scalar(noeds,nprocs,ncube,nbyte,nstart_h,nend_h&
            &,nlocal_h,nstart_ah,nend_ah,nlocal_ah,ngrp,grppt,worka&
            &,cpu_h,ncpu_h)

#endif

!!$*----- Split bonded interaction and constraints (nstart_1, nstart_ex)
       nprot=ntap
       count=0
       DO i=1,ntap
          protl(1+count)=1
          protl(2+count)=i
          count=count+2
       END DO
       count=0
       DO i=1,nprot
          m=protl(count+1)
          DO i1=1,m
             j=protl(count+1+i1)
             atomp(j)=i
          END DO
          count=count+m+1
       END DO
       CALL P_atoms_split_intra(noeds,nprocs,ncube,nbyte,nstart_1&
            &,nend_1,nlocal_1,nstart_ex0,nend_ex0,nlocal_ex0&
            &,nstart_ex,nend_ex,nlocal_ex,nprot,protl,cnstp&
            &,cnst_protp,cnst_protl,cnst_protp_1,cnst_protl_1&
            &,cnst_protp_ex,cnst_protl_ex,ntap,lstrtch,lstretch,lbndg&
            &,lbend,atomp,iret,errmsg) 

!!$*----- Split intramolecular interaction arrays

      CALL P_split_intra(noeds,nprocs,ncube,nstart_1,nend_1,worka,iret&
           &,errmsg) 

      IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)

!!$*----- Split nonbonded interaction loop and do set up (nstart_2) 

      nstart_2=nstart_1
      nend_2=nend_1
      nlocal_2=nlocal_1

      IF(nfft1 == 0 .OR. nfft2 == 0 .OR. nfft3 == 0 ) THEN
         errmsg='Cannot run with zero grid length'
         CALL xerror(errmsg,80,1,2)
      END IF
      IF(pme_order == 0) THEN
         errmsg='Cannot run with zero spline order'
         CALL xerror(errmsg,80,1,2)
      END IF
      CALL fft_pme_init(ntap,nfft1,nfft2,nfft3,pme_order,sizfftab&
           &,sizffwrk,siztheta,siz_Q,sizheap,sizstack,bsp_mod1&
           &,bsp_mod2,bsp_mod3,fftable,ffwork)
      rshk='l'
      IF(erf_corr) THEN
         CALL erf_corr_cutoff(oc,delew,rkcut,nfft1,nfft2,nfft3)
         CALL int_corr_erf_spline(rlew,ruew,nbinew,alphal,rkcut&
              &,erf_arr_corr,work)
      END IF
      CALL Pme_init(noeds,nprocs,nodex,nodey,nodez,npy,npz,ictxt,descQ&
           &,fftable,nfft1,nfft2,nfft3,nfft3_start,nfft3_local&
           &,nfft2_start,nfft2_local,iret,errmsg)
      IF(iret .EQ. 1) CALL xerror(errmsg,80,1,2)
      IF(linked_cell) THEN
         indxyz=(2*(ncx-1)+1)*(2*(ncy-1)+1)*(2*(ncz-1)+1)
         ALLOCATE(indxi(indxyz),indxj(indxyz),indxk(indxyz))
      END IF
      ALLOCATE(tag_bndg(lbend))
      CALL gen_stress_tag_bnd(lbend,3,lbndg,lbndg_x,atomp,tag_bndg)

#ifdef PARALLEL
      CALL REDUCE_Init(nstart_2,nend_2,nlocal_2,noeds,nprocs)
#endif
      CALL mapnl_divide(node,nstart_h,nend_h,grppt,mapnl)

    END SUBROUTINE SetUp
  END SUBROUTINE Init
  INCLUDE 'FIELD__Read.f90'
END MODULE FIELD_Mod
