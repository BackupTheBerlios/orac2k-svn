!!$***********************************************************************
!!$   Time-stamp: <2005-02-25 14:09:35 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Mon Feb 14 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program ORAC ----*
MODULE Module_Fourier
  USE xerror_mod
  TYPE Input
     INTEGER :: label=-1
     INTEGER, POINTER :: node,nodex,nodey,nodez,ictxt,npy,npz&
          &,descQ(:),nprocs,ncube,nbyte,rbyte,nstart_2,nend_2&
          &,nlocal_2,ntap, order,nfft1,nfft2,nfft3,nfft3_start&
          &,nfft3_local,nfft2_start,nfft2_local,atomp(:),grppt(:,:)
     REAL(8), POINTER :: xp(:),yp(:),zp(:),xpcm(:),ypcm(:),zpcm(:)&
          &,pmechg(:),co(:,:),oc(:,:),volume,alphal,eer,fpx(:),fpy(:)&
          &,fpz(:),phi(:),stressc(:,:),rkcut
     LOGICAL, POINTER :: pressure
  END TYPE Input
  TYPE Fourier_Out
     REAL(8), POINTER  :: energy
     REAL(8), POINTER  :: fx(:),fy(:),fz(:)
     REAL(8), POINTER  :: phi(:)
  END TYPE Fourier_Out
  TYPE (Input), PRIVATE, SAVE, TARGET :: s
  PRIVATE :: Input
CONTAINS
  SUBROUTINE Init(node,nodex,nodey,nodez,ictxt,npy,npz,descQ,nprocs&
       &,ncube,nbyte,rbyte,nstart_2,nend_2,nlocal_2,ntap,xp,yp,zp&
       &,xpcm,ypcm,zpcm,pmechg,co,oc,volume,alphal,order,nfft1,nfft2&
       &,nfft3,nfft3_start,nfft3_local,nfft2_start,nfft2_local,eer&
       &,fpx,fpy,fpz,phi,stressc,atomp,grppt,pressure,rkcut)

!!$======================== DECLARATIONS ================================*

    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

    INTEGER, TARGET :: NFFT1,NFFT2,NFFT3,nfft3_start,nfft3_local,nfft2_start&
         &,nfft2_local,NUMATOMS,ORDER,node,nprocs,ncube,nbyte,rbyte&
         &,nodex,nodey,nodez,ictxt,npy,npz,descQ(:)
    INTEGER, TARGET :: ntap,atomp(:),grppt(:,:),count,mia,ia,nstart_2,nend_2,nlocal_2
    REAL(8), TARGET :: fpx(:),fpy(:),fpz(:),phi(:)
    REAL(8), TARGET :: xp(:),yp(:),zp(:),xpcm(:),ypcm(:),zpcm(:),pmechg(:)
    REAL(8), TARGET :: co(3,3),oc(3,3),volume,alphal,stressc(3,3),eer,rkcut
    LOGICAL, TARGET :: pressure
    
!!$----------------------- VARIABLES IN COMMON --------------------------*
    
!!$------------------------- LOCAL VARIABLES ----------------------------*
    
    INTEGER :: n,i
    
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*
    
    IF(s % label == -1) THEN
       s % label = 0
    END IF
    s % node => node
    s % nodex => nodex
    s % nodey => nodey
    s % nodez => nodez
    s % ictxt => ictxt
    s % npy => npy
    s % npz => npz
    s % descQ => descQ
    s % nprocs => nprocs
    s % ncube => ncube
    s % nbyte => nbyte
    s % rbyte => rbyte
    s % nstart_2 => nstart_2
    s % nend_2 => nend_2
    s % nlocal_2 => nlocal_2
    s % ntap => ntap
    s % xp => xp
    s % yp => yp
    s % zp => zp
    s % xpcm => xpcm
    s % ypcm => ypcm
    s % zpcm => zpcm
    s % pmechg => pmechg
    s % co => co
    s % oc => oc
    s % volume => volume
    s % alphal => alphal
    s % order => order
    s % nfft1 => nfft1
    s % nfft2 => nfft2
    s % nfft3 => nfft3
    s % nfft3_start => nfft3_start
    s % nfft3_local => nfft3_local
    s % nfft2_start => nfft2_start
    s % nfft2_local => nfft2_local
    s % eer => eer
    s % fpx => fpx
    s % fpy => fpy
    s % fpz => fpz
    s % phi => phi
    s % stressc => stressc
    s % atomp => atomp
    s % grppt => grppt
    s % pressure => pressure
    s % rkcut => rkcut
  END SUBROUTINE Init
  FUNCTION Do_Fourier() RESULT (out)

!!$======================== DECLARATIONS ================================*

    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

    TYPE(Fourier_Out) :: out
    
!!$------------------------- LOCAL VARIABLES ----------------------------*
    
    INTEGER :: n,i
    REAL(8) :: a

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    IF(s % label == -1) CALL abort_now('Gaussian Input type not initia&
         &lized') 
    CALL fft_pme(s % node,s % nodex,s % nodey,s % nodez,s % ictxt,s %&
         & npy,s % npz,s % descQ,s % nprocs,s % ncube,s % nbyte,s %&
         & rbyte,s % nstart_2,s % nend_2,s % nlocal_2,s % ntap,s % xp&
         &,s % yp,s % zp,s % xpcm,s % ypcm,s % zpcm,s % pmechg,s % co&
         &,s % oc,s % volume,s % alphal,s % order,s % nfft1,s % nfft2&
         &,s % nfft3,s % nfft3_start,s % nfft3_local,s % nfft2_start&
         &,s % nfft2_local,s % eer,s % fpx,s % fpy,s % fpz,s % phi,s &
         &% stressc,s % atomp,s % grppt,s % pressure,s % rkcut)

    out % fx     =>s % fpx
    out % fy     =>s % fpy
    out % fz     =>s % fpz
    out % phi    =>s % phi
    out % energy =>s % eer
  END FUNCTION Do_Fourier
  SUBROUTINE Modify_Charges(chg,ind) 
!!$======================== DECLARATIONS ================================*
    IMPLICIT none
!!$----------------------------- ARGUMENTS ------------------------------*
    REAL(8) :: chg (:)
    INTEGER :: ind(:)
!!$------------------------- LOCAL VARIABLES ----------------------------*
    INTEGER :: n,i,ii
    REAL(8) :: test_a
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    IF(s % label == -1) CALL abort_now('Gaussian Input type not initia&
         &lized') 

    test_a=s % pmechg(1)
    n=SIZE(ind)
!!$    s % pmechg = 0.0D0
!!$    s % pmechg(1)=test_a
    DO ii=1,n
       i=ind(ii)
       s % pmechg(i) = chg(ii)
    END DO

  END SUBROUTINE Modify_Charges
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE Module_Fourier
