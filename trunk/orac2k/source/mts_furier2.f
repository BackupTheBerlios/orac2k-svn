      SUBROUTINE mts_furier2(node,nodex,nodey,nodez,ictxt,npy,npz,descQ
     &     ,nprocs,ncube,nstart_1,nend_1,nlocal_1,nstart_2,nend_2
     &     ,nlocal_2,xpb,ypb,zpb,xp0,yp0,zp0,xpcm,ypcm,zpcm,urcsp,urcs
     &     ,urcp,fpx,fpy,fpz,Ex_rec,Ey_rec,Ez_rec,Ex_Cor,Ey_Cor,Ez_Cor
     &     ,ene,eer,fudgec,tag_bndg,ene_ferrf,charges,dipole)

************************************************************************
*                                                                      *
*     Transit routine : calls MTS_FURIWW, MTS_FURIPP, MTS_FURIPW,      *
*     FERFF, WATSELF, PME_FFT                                          *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     Written by Piero Procacci, CECAM-ENS Lyon 1995                   *
*     Modified for parallel machines by                                *
*              Massimo Marchi                                          *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed Feb 10 1999 -                                     *
*                                                                      *
*                                                                      *
************************************************************************


*======================= DECLARATIONS ==================================

      IMPLICIT none
      INTERFACE
         SUBROUTINE P_Reduce_Forces(x,y,z,nstart,nend,nlocal,node,nprocs
     &        )
         REAL(8) :: x(*),y(*),z(*)
         INTEGER, OPTIONAL :: nstart,nend,nlocal,node,nprocs
         END SUBROUTINE P_Reduce_Forces
      END INTERFACE

*----------------------- ARGUMENTS -------------------------------------

      INTEGER node,nodex,nodey,nodez,ictxt,npy,npz,ncube,nprocs,descQ(*)
      INTEGER nstart_2,nend_2,nlocal_2,nstart_1,nend_1,nlocal_1
      REAL*8  xp0(*),yp0(*),zp0(*),fpx(*),fpy(*),fpz(*),ene(*)
      REAL*8  Ex_rec(*),Ey_rec(*),Ez_rec(*),dipole(3,*),charges(*)
      REAL*8  xpb(*),ypb(*),zpb(*),urcsp,urcs,urcp,fudgec
      REAL*8 Ex_Cor(*),Ey_Cor(*),Ez_Cor(*)
      REAL*8  xpcm(*),ypcm(*),zpcm(*)
      REAL*8  eer,ene_ferrf
      INTEGER tag_bndg(*)

c---  stuff for b-spline interpolation and FFT 

*---- VARIABLES IN INCLUDE --------------------------------------------*

      INCLUDE 'parst.h'
      INCLUDE 'cpropar.h'
      INCLUDE 'pme.h'
      REAL*8  fxx(m1,2),fyy(m1,2),fzz(m1,2),Ex(m1),Ey(m1),Ez(m1)
      SAVE fxx,fyy,fzz

      REAL*8  fsbond_slt,fsbend_slt,fsin14_slt,fsbond_slv,fsbend_slv
     &     ,fsin14_slv,coul_bnd_slt,coul_bnd_slv,cnb14_slt,cnb14_slv
     &     ,ene_dd,ene_cc,Utotal
      LOGICAL self

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER  i,j
#if defined T3E | defined _CRAY_
#define _AXPY_  saxpy
      INTEGER*8 n_loc,one
#else
#define _AXPY_  daxpy
      INTEGER n_loc,one
#endif
      DATA one/1/

*==== EXECUTABLE STATEMENTS: ==========================================*

      self=.TRUE.
      n_loc=nlocal_1

*=======================================================================
*      subtract "intramolecular term" in the ZERO cell:   BONDS     ----
*=======================================================================

      Ex_Cor(1:ntap)=0.0D0
      Ey_Cor(1:ntap)=0.0D0
      Ez_Cor(1:ntap)=0.0D0
      ene_cc=0.0D0
      ene_dd=0.0D0
      ene_ferrf=0.0D0
      CALL zeroa(fxx(nstart_1,1),fyy(nstart_1,1),fzz(nstart_1,1)
     &     ,nlocal_1,1)
      CALL zeroa(fxx(nstart_1,2),fyy(nstart_1,2),fzz(nstart_1,2)
     &     ,nlocal_1,1)

      IF(lstretch .NE. 0) THEN
         CALL ferrf(ss_index,alphal,charges,1.0D0,xpb,ypb,zpb,1,lstrtch
     &        ,lstretch,lbnd_x,fsbond_slt,fsbond_slv,fxx(1,2),fyy(1,2)
     &        ,fzz(1,2),erf_corr,erf_arr_corr,delew,rlew)
         CALL Correct_Ewald(alphal,charges,dipole,xpb,ypb,zpb,lcnstr
     &        ,lconstr,lconstr_x,Utotal,Ex_Cor,Ey_Cor,Ez_Cor,fpx,fpy,fpz
     &        )
         ene_ferrf=ene_ferrf+fsbond_slt+fsbond_slv+Utotal
      END IF
      IF(lconstr .NE. 0) THEN
c$$$         CALL ferrf(ss_index,alphal,charges,1.0D0,xpb,ypb,zpb,1,lcnstr
c$$$     &        ,lconstr,lconstr_x,fsbond_slt,fsbond_slv,Ex_Cor,Ey_Cor
c$$$     &        ,Ez_Cor,erf_corr,erf_arr_corr,delew,rlew)
c$$$         Ex_Cor(1:ntap)=Ex_Cor(1:ntap)/charges(1:ntap)
c$$$         Ey_Cor(1:ntap)=Ey_Cor(1:ntap)/charges(1:ntap)
c$$$         Ez_Cor(1:ntap)=Ez_Cor(1:ntap)/charges(1:ntap)
c$$$
c$$$         CALL Correct_Ewald(alphal,charges,dipole,xpb,ypb,zpb,lcnstr
c$$$     &        ,lconstr,lconstr_x,Utotal,fpx,fpy,fpz)

         CALL Correct_Ewald(alphal,charges,dipole,xpb,ypb,zpb,lcnstr
     &        ,lconstr,lconstr_x,Utotal,Ex_Cor,Ey_Cor,Ez_Cor,fpx,fpy
     &        ,fpz)
         ene_ferrf=ene_ferrf+Utotal
      END IF
*=======================================================================
*----- Subtract "intramolecular term" in the ZERO cell:   BENDS     ----
*=======================================================================
      
      IF(int13p .NE. 0) THEN
         CALL ferrf_tag(ss_index,alphal,chrge,xpb,ypb,zpb,int13
     &        ,int13p,int13_x,fsbend_slt,fsbend_slv,fxx,fyy,fzz,m1
     &        ,tag_bndg,erf_corr,erf_arr_corr,delew,rlew)
         CALL Correct_Ewald(alphal,charges,dipole,xpb,ypb,zpb,int13
     &        ,int13p,int13_x,Utotal,Ex_Cor,Ey_Cor,Ez_Cor,fxx(1,2),fyy(1
     &        ,2),fzz(1,2))
         ene_ferrf=ene_ferrf+fsbend_slt+fsbend_slv+Utotal
      END IF

*=======================================================================
*----- Subtract "intramolecular term" in the ZERO cell:   1-4 Fudged ---
*=======================================================================
      
      IF(int14p .NE. 0) THEN
         CALL ferrf(ss_index,alphal,charges,fudgec,xpb,ypb,zpb,1,int14
     &        ,int14p,int14_x,fsin14_slt,fsin14_slv,fxx(1,2),fyy(1,2)
     &        ,fzz(1,2),erf_corr,erf_arr_corr,delew,rlew)
         
         CALL fnb14_2(ss_index,xpb,ypb,zpb,charges,ntap,alphal
     &        ,int14,int14p,type14,int14_x,fudge
     &        ,fxx(1,2),fyy(1,2),fzz(1,2),cnb14_slt,cnb14_slv)
         ene_ferrf=ene_ferrf+fsin14_slv+cnb14_slv+fsin14_slt+cnb14_slt
      END IF
*=======================================================================
*---- Add up all contributions from tagged ferrf and then to         ---
*---- full force                                                     ---
*=======================================================================

      CALL _AXPY_(n_loc,1.0D0,fxx(nstart_1,2),one,fpx(nstart_1),one)
      CALL _AXPY_(n_loc,1.0D0,fyy(nstart_1,2),one,fpy(nstart_1),one)
      CALL _AXPY_(n_loc,1.0D0,fzz(nstart_1,2),one,fpz(nstart_1),one)
      CALL _AXPY_(n_loc,1.0D0,fxx(nstart_1,1),one,fpx(nstart_1),one)
      CALL _AXPY_(n_loc,1.0D0,fyy(nstart_1,1),one,fpy(nstart_1),one)
      CALL _AXPY_(n_loc,1.0D0,fzz(nstart_1,1),one,fpz(nstart_1),one)
      
      Ex(nstart_1:nend_1)=fxx(nstart_1:nend_1,1)+fxx(nstart_1:nend_1,2)
      Ey(nstart_1:nend_1)=fyy(nstart_1:nend_1,1)+fyy(nstart_1:nend_1,2)
      Ez(nstart_1:nend_1)=fzz(nstart_1:nend_1,1)+fzz(nstart_1:nend_1,2)


*=======================================================================
*---- Merge energies and expands forces in parallel runs             ---
*=======================================================================

*=======================================================================
*--- Call fft_pme                                                   ----
*=======================================================================

      CALL fft_pme_dipole(node,nodex,nodey,nodez,ictxt,npy,npz,descQ
     &     ,nprocs,ncube,nbyte,rbyte,nstart_2,nend_2,nlocal_2,ntap
     &     ,xp0,yp0,zp0,xpcm,ypcm,zpcm,charges,dipole,co,oc,volume
     &     ,alphal,pme_order,nfft1,nfft2,nfft3,nfft3_start,nfft3_local
     &     ,nfft2_start,nfft2_local,eer,fpx,fpy,fpz,ene,Ex_rec,Ey_rec
     &     ,Ez_rec,rkcut)
      
#ifdef PARALLEL
c$$$
c$$$--- Not yet parallel-ready!!!
c$$$--- Not yet parallel-ready!!!
c$$$--- Not yet parallel-ready!!!
c$$$

      CALL P_Reduce_Forces(fpx,fpy,fpz)
      IF(nprocs .GT. 1) THEN
         CALL P_merge_r8(coul_bnd_slt,node,nprocs,ncube,rbyte)
         CALL P_merge_r8(coul_bnd_slv,node,nprocs,ncube,rbyte)
      END IF
#endif         
      Ex_Cor(1:ntap)=Ex_Cor(1:ntap)+Ex(1:ntap)/charges(1:ntap)
      Ey_Cor(1:ntap)=Ey_Cor(1:ntap)+Ey(1:ntap)/charges(1:ntap)
      Ez_Cor(1:ntap)=Ez_Cor(1:ntap)+Ez(1:ntap)/charges(1:ntap)

      

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
