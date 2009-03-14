MODULE ElecPotential_Mod

!!$***********************************************************************
!!$   Time-stamp: <04/12/16 14:56:07 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Mar 30 2001 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program ORAC ----*

  REAL(8), DIMENSION (:,:,:), ALLOCATABLE, SAVE :: Phi,Rho,Phi0,Rho0 &
       ,HardSphere

  INTEGER, SAVE :: ntap,ngrp,nstart_h,nend_h,nlocal_h &
       ,nstart_ah,nend_ah,nlocal_ah,nprot,nprocs,node,nfft1,nfft2,nfft3 &
       ,nodex,nodey,nodez,npy,npz,nay,naz,order,kstart,kend,jstart,jend &
       ,ncube,rbyte,nbyte,CallsToRoutine=0,NatoSlt,nf1,nf2,nf3,UnitPhi  &
       ,UnitPhiDat,RecordLength,UnitPhiCoord
  REAL(8), SAVE :: ewald_coeff,rkcut,efact,unitc &
       ,RadiusToSigma=1.78179743628067860D0,Bin1,Bin2,Bin3              &
       ,SigmaIonInit,oc(3,3),kT,volume,recip(3,3),HardSphereCutoff      &
       ,SigmaIon,ChargeIon,a=0.0D0,b=0.0D0,c=0.0D0,alfa=0.0D0           &
       ,beta=0.0D0,gamma=0.0D0,SmoothFactor,EnergyVacuum=0.0D0

  REAL(8), DIMENSION (:), ALLOCATABLE, SAVE :: charge,ChargesWithIon,pmass  &
       ,mass,pnbd1
  INTEGER, DIMENSION (:,:), ALLOCATABLE, SAVE :: grppt
  REAL(8), DIMENSION (:,:), ALLOCATABLE, SAVE :: bsp_mod1 &
       ,bsp_mod2,bsp_mod3
  INTEGER, DIMENSION (:), ALLOCATABLE, SAVE :: ss_index,protl,nbtype
  INTEGER, DIMENSION (:), POINTER, SAVE :: IndX,IndY,IndZ

  REAL(8), DIMENSION (:), ALLOCATABLE, SAVE :: xpa,ypa,zpa,xpg,ypg,zpg &
       ,xpga,ypga,zpga,xpcm,ypcm,zpcm,xpcma,ypcma,zpcma,fr1,fr2,fr3,Sigma &
       ,XpAvg,YpAvg,ZpAvg
  REAL(8), DIMENSION (:,:), ALLOCATABLE, SAVE :: theta1,theta2,theta3,dtheta1 &
       ,dtheta2,dtheta3
  INTEGER(2), DIMENSION (:,:), ALLOCATABLE, SAVE :: indk1,indk2,indj1,indj2 &
       ,mk,mj
  CHARACTER(80), SAVE :: FileNamePhi,FileNamePhiDat,FileNamePhiCoord
  CHARACTER(7), DIMENSION(:), ALLOCATABLE, SAVE :: Label
  LOGICAL :: DoFreeEnergy=.FALSE.

END MODULE ElecPotential_Mod
