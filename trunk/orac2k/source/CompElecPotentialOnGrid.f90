SUBROUTINE CompElecPotentialOnGrid(co,xp0,yp0,zp0,AChargeIon,ASigmaIon   &
     ,AKElecPot,ASmoothFactor,ADoFreeEnergy,ALabel,Apmass,Ass_index      &
     ,Amass,Anfft1,Anfft2,Anfft3,Anprocs,Anode &
     ,Ancube,Arbyte,Anbyte,Antap,Angrp,Agrppt,Acharge,Anstart_h,Anend_h  &
     ,Anlocal_h,Anstart_ah,Anend_ah,Anlocal_ah,Anprot,Aprotl,Aorder      &
     ,Aewald_coeff,Absp_mod1,Absp_mod2,Absp_mod3,Arkcut,Aefact,Aunitc    &
     ,AnatoSlt,Anbtype,Apnbd1)

!!$***********************************************************************
!!$   Time-stamp: <01/04/13 10:06:06 marchi>                           *
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

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*

  USE ElecPotential_Mod
  IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

  INTEGER, OPTIONAL :: Antap,Angrp,Agrppt(2,*),Anstart_h,Anend_h,Anlocal_h  &
       ,Anstart_ah,Anend_ah,Anlocal_ah,Anprot,Aprotl(*),Ass_index(*),Anfft1 &
       ,Anfft2,Anfft3,Anprocs,Anode,Ancube,Arbyte,Anbyte,Aorder,AnatoSlt    &
       ,Anbtype(*),AKElecPot
  REAL(8), OPTIONAL :: AChargeIon,ASigmaIon,Aewald_coeff,Acharge(*)         &
       ,Amass(*),Apmass(*),Absp_mod1(*),Absp_mod2(*),Absp_mod3(*),Arkcut    &
       ,Aefact,Aunitc,Apnbd1(*),ASmoothFactor
  CHARACTER(7), OPTIONAL :: ALabel(*)
  LOGICAL, OPTIONAL :: ADoFreeEnergy
  REAL(8) :: co(3,3),xp0(*),yp0(*),zp0(*)

  INTERFACE
     SUBROUTINE GetHardSphere(co,HardSPhere)
       REAL(8), DIMENSION (:,:,:) :: HardSphere
       REAL(8) :: co(3,3)
     END SUBROUTINE GetHardSphere
     SUBROUTINE ScalarSum(Q,Rho_k,eer)
       REAL(8), DIMENSION (:,:,:,:) :: Q,Rho_k
       REAL(8) :: eer       
     END SUBROUTINE ScalarSum
     SUBROUTINE do_pfft3d(isign,z,na1,na2,na3,nda1,nda2,nda3,iproca,nproca)
       INTEGER, INTENT(in) :: isign
       INTEGER, INTENT(in), OPTIONAL :: na1,na2,na3,iproca,nproca
       REAL(8), DIMENSION (*), INTENT(inout) :: z
       INTEGER, INTENT(out), OPTIONAL :: nda1,nda2,nda3
     END SUBROUTINE do_pfft3d
     SUBROUTINE LinkeCellIndex(ncx,ncy,ncz,IndX,IndY,IndZ,ctoff,co)
       INTEGER :: ncx,ncy,ncz,nind
       INTEGER, DIMENSION (:), POINTER :: IndX,IndY,IndZ
       REAL(8)  :: ctoff,co(3,3)
     END SUBROUTINE LinkeCellIndex
  END INTERFACE
!!$------------------------- LOCAL VARIABLES ----------------------------*

  INTEGER :: count,j,i,m,Natom,ierr,k,k0,ith1,ith2,ith3,j0,i0,IMax,JMax &
       ,KMax,IIon,JIon,KIon,nka,nia,nja,ka,ka0,ia,ia0,ja,ja0
  REAL(8), SAVE :: pi=3.14159265358979323846D0,CFact
  REAL(8), SAVE :: virial(3,3) &
       ,U_eps,U_eps0,FreeEnergy,DeltaG &
       ,Dist,xc,yc,zc,Rmax,TempExp,TempPhi,TempHS &
       ,ctoff,a0,b0,c0,alfa0,beta0,gamma0,TempFact,Tmass,xpcm0,ypcm0,zpcm0

  REAL(8), DIMENSION (:,:,:,:), ALLOCATABLE, SAVE :: Q,Rho_k
  REAL(8), DIMENSION (:,:,:), ALLOCATABLE, SAVE :: TempHardSphere
  REAL(8) :: dummy1,dummy2,time_a,time_b,dv

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*


  IF(present(AChargeIon)) THEN
     CALL InitializeRoutineVariables
     RETURN
  END IF

  CallsToRoutine=CallsToRoutine+1

!!$*=======================================================================
!!$*--- Compute the volume of the system ----------------------------------
!!$*=======================================================================

  CALL matinv(3,3,co,oc,volume)
  volume=volume*2.0D0**3
  CALL rotb(a0,b0,c0,alfa0,beta0,gamma0,co)
  a=a+a0
  b=b+b0
  c=c+c0
  alfa=alfa+alfa0
  beta=beta+beta0
  gamma=gamma+gamma0

!!$*=======================================================================
!!$*-------- Calculate group position  ------------------------------------
!!$*=======================================================================
                  
  CALL appbou(xp0,yp0,zp0,xpg,ypg,zpg,pmass,1,ngrp,grppt)

!!$*=======================================================================
!!$*-------- Accumulate coordinates for average ---------------------------
!!$*=======================================================================
  
  XpAvg=XpAvg+xp0(1:NatoSlt)
  YpAvg=YpAvg+yp0(1:NatoSlt)
  ZpAvg=ZpAvg+zp0(1:NatoSlt)
  Tmass=SUM(mass(1:NatoSlt))
  xpcm0=SUM(xp0(1:NatoSlt)*mass(1:NatoSlt))/Tmass
  ypcm0=SUM(yp0(1:NatoSlt)*mass(1:NatoSlt))/Tmass
  zpcm0=SUM(zp0(1:NatoSlt)*mass(1:NatoSlt))/Tmass

!!$*=======================================================================
!!$*-------- Calculate solute center of mass ------------------------------
!!$*=======================================================================
                  
  CALL inicmp(ss_index,xp0,yp0,zp0,xpcm,ypcm,zpcm,mass,nprot,protl)
                  
!!$*=======================================================================
!!$*--- Change frame to get xpa, ypa, zpa etc in box fractions ------------
!!$*=======================================================================
                  
  CALL change_frame(co,oc,-1,ntap,xp0,yp0,zp0,xpa,ypa,zpa)
  CALL change_frame(co,oc,-1,ngrp,xpg,ypg,zpg,xpga,ypga,zpga)
  CALL change_frame(co,oc,-1,nprot,xpcm,ypcm,zpcm,xpcma,ypcma,zpcma)

!!$======================================================================
!!$--- Get cell indeces for the step ------------------------------------
!!$======================================================================

  CALL LinkeCellIndex(nfft1,nfft2,nfft3,IndX,IndY,IndZ,HardSphereCutoff,co)

!!$======================================================================
!!$--- Get accessible space ---------------------------------------------
!!$======================================================================

  ALLOCATE(TempHardSphere(nfft1,nfft2,naz))

  CALL GetHardSphere(co,TempHardSphere)
  HardSphere(:,:,kstart:kend)=HardSphere(:,:,kstart:kend) &
       +TempHardSphere(:,:,1:naz)
  DEALLOCATE(TempHardSphere)

!!$======================================================================
!!$--- PME contribution on the reciprocal space -------------------------
!!$======================================================================

  recip=oc*0.5D0
  CALL GetScaledFractions

  ALLOCATE(Q(2,nfft1,nfft2,naz))
  ALLOCATE(Rho_k(2,nfft1,nfft2,naz))

  Rho_k=0.0D0
  Q=0.0D0

!!$======================================================================
!!$--- PME contribution on the reciprocal space -------------------------
!!$======================================================================

  charge(1:NatoSlt)=0.0D0
  CALL GetEnergyPhiAndRho(ntap,U_eps)
  charge(1:NatoSlt)=ChargesWithIon(1:NatoSlt)

!!$======================================================================
!!$--- Smooth electrostatic potential with Gaussian ---------------------
!!$======================================================================

  CALL GetSmoothingGaussian(Q)
  CALL do_pfft3d(-1,Q)
  TempFact=Cfact*efact/1000.0D0
  Phi(:,:,kstart:kend)=Phi(:,:,kstart:kend)+TempFact*Q(1,:,:,1:naz)

!!$======================================================================
!!$--- Smooth electrostatic charges with Gaussian -----------------------
!!$======================================================================

  CALL GetSmoothingGaussian(Rho_k)
  CALL do_pfft3d(-1,Rho_k)
  TempFact=1.0D0/Cfact
  Rho(:,:,kstart:kend)=Rho(:,:,kstart:kend)+TempFact*Rho_k(1,:,:,1:naz)

!!$========================================================================
!!$--- Compute Energy in vacuum -------------------------------------------
!!$========================================================================

  Rho_k=0.0D0
  Q=0.0D0

  CALL GetEnergyPhiAndRho(NatoSlt,U_eps0)

  CALL do_pfft3d(-1,Rho_k)
  CALL do_pfft3d(-1,Q)

  TempFact=Cfact*efact/1000.0D0
  Phi(:,:,kstart:kend)=Phi(:,:,kstart:kend)+TempFact*Q(1,:,:,1:naz)
  Phi0(:,:,kstart:kend)=Phi0(:,:,kstart:kend)+TempFact*Q(1,:,:,1:naz)
  TempFact=1.0D0/Cfact
  Rho(:,:,kstart:kend)=Rho(:,:,kstart:kend)+TempFact*Rho_k(1,:,:,1:naz)
  Rho0(:,:,kstart:kend)=Rho0(:,:,kstart:kend)+TempFact*Rho_k(1,:,:,1:naz)

!!$========================================================================
!!$--- Evaluate Free Energy -----------------------------------------------
!!$========================================================================

  IF(DoFreeEnergy) THEN
     CALL GetFreeEnergy(FreeEnergy,EnergyVacuum)
     DeltaG=FreeEnergy-EnergyVacuum
     IF(node == 0) WRITE(*,'(''DeltaG '',3f14.5)') DeltaG,U_eps,U_eps0
  END IF


!!$========================================================================
!!$--- Save Result --------------------------------------------------------
!!$========================================================================

  IF(MOD(CallsToRoutine,10) == 0) THEN
     OPEN(unit=UnitPhi,file=FileNamePhi,access='DIRECT' &
          ,form='UNFORMATTED',status='UNKNOWN',recl=RecordLength)
     TempFact=1.0D0/DFLOAT(CallsToRoutine)
     DO k0=1,naz
        k=node*naz+k0
        WRITE(UnitPhi,rec=k) HardSphere(:,:,k)*Tempfact &
             ,Phi0(:,:,k)*TempFact,Rho0(:,:,k)*TempFact &
             , Phi(:,:,k)*TempFact, Rho(:,:,k)*TempFact 
     END DO
     IF(node == 0) THEN
        OPEN(unit=UnitPhiCoord,file=FileNamePhiCoord,form='UNFORMATTED' &
             ,status='UNKNOWN')
        WRITE(UnitPhiCoord) a*TempFact,b*TempFact,c*Tempfact
        WRITE(UnitPhiCoord) alfa*TempFact,beta*TempFact,gamma*TempFact
        WRITE(UnitPhiCoord) XpAvg*TempFact,YpAvg*TempFact,ZpAvg*TempFact
        WRITE(UnitPhiCoord) Label(1:NatoSlt),Charge(1:NatoSlt)*DSQRT(unitc) &
             ,Sigma(1:NatoSlt)
     END IF
     CLOSE(UnitPhi)
     CLOSE(UnitPhiCoord)
  END IF
  DEALLOCATE(Q)
  DEALLOCATE(Rho_k)

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*
CONTAINS
  SUBROUTINE GetEnergyPhiAndRho(Natom,Energy)
    IMPLICIT none
    INTEGER :: Natom
    REAL(8) :: Energy,eer


    CALL get_bspline_coeffs(Natom,fr1,fr2,fr3,order,theta1,theta2,theta3 &
         ,dtheta1,dtheta2,dtheta3,nfft1,nfft2,nfft3,kstart,kend,jstart,jend &
         ,indk1,indk2,indj1,indj2,mk,mj)

    CALL fill_charge_grid(node,nodey,nodez,nay,naz,Natom,charge,theta1 &
         ,theta2,theta3,fr1,fr2,fr3,order,nfft1,nfft2,nfft3,nfft1,nfft2   &
         ,naz,Q,indk1,indk2,indj1,indj2,mk,mj)

    CALL do_pfft3d(1,Q)

    CALL ScalarSum(Q,Rho_k,eer)

    IF(nprocs > 1) THEN
       CALL P_merge_r8(eer)
    END IF
!!$  eer=eer-ewald_coeff*SUM(charge(1:Natom)**2)/DSQRT(pi)
    Energy=eer*efact/1000.0D0

  END SUBROUTINE GetEnergyPhiAndRho

  SUBROUTINE GetFreeEnergy(FreeEnergy,EnergyInVacuum)
    IMPLICIT NONE 
    REAL(8) :: FreeEnergy,EnergyInVacuum
    REAL(8) :: FreeEnergyOnNode,EnergyInVacuumOnNode,dv

    dv=volume/DFLOAT(nfft1*nfft2*nfft3)
    FreeEnergyOnNode=SUM(Rho0(:,:,kstart:kend)*Phi(:,:,kstart:kend)) &
         /DFLOAT(CallsToRoutine**2)
    EnergyInVacuumOnNode=SUM(Rho0(:,:,kstart:kend)*Phi0(:,:,kstart:kend)) &
         /DFLOAT(CallsToRoutine**2)
    IF(nprocs > 1) THEN
       CALL P_merge_r8(FreeEnergyOnNode)
       CALL P_merge_r8(EnergyInVacuumOnNode)
    END IF
    FreeEnergy=FreeEnergyOnNode*dv*0.5D0
    EnergyInVacuum=EnergyInVacuumOnNode*dv*0.5D0
  END SUBROUTINE GetFreeEnergy
  SUBROUTINE GetScaledFractions
    IMPLICIT NONE 
    INTEGER :: n,k,k0
    REAL(8) :: w1,w2,w3
    
    DO n=1,ntap
        w1 = xp0(n)*recip(1,1)+yp0(n)*recip(2,1)+zp0(n)*recip(3,1)
        w2 = xp0(n)*recip(1,2)+yp0(n)*recip(2,2)+zp0(n)*recip(3,2)
        w3 = xp0(n)*recip(1,3)+yp0(n)*recip(2,3)+zp0(n)*recip(3,3)
        fr1(n) = nfft1*(w1 - anint(w1)+0.5D0)
        fr2(n) = nfft2*(w2 - anint(w2)+0.5D0)
        fr3(n) = nfft3*(w3 - anint(w3)+0.5D0)
    END DO
  END SUBROUTINE GetScaledFractions
  SUBROUTINE GetSmoothingGaussian(Phi)
    REAL(8), DIMENSION (:,:,:,:) :: Phi
    REAL(8) :: fac,Rc
    INTEGER  ::  k,k1,k2,k3,m1,m2,m3
    INTEGER  ::  ka3,ka2
    REAL(8)  ::  mhat1,mhat2,mhat3,msq

    Rc=SmoothFactor
    fac = 2.0D0*Rc*pi**2
    DO ka3=1,naz
       k3=nodez*naz+ka3
       m3 = k3 - 1
       if ( k3 .gt. nf3 )m3 = k3 - 1 - nfft3

       DO ka2=1,nay
          k2=nodey*nay+ka2
          m2 = k2 - 1
          if ( k2 .gt. nf2 )m2 = k2 - 1 - nfft2
          
          DO k1=1,nfft1
             m1 = k1 - 1
             if ( k1 .gt. nf1 )m1 = k1 - 1 - nfft1

             mhat1 = recip(1,1)*m1+recip(1,2)*m2+recip(1,3)*m3
             mhat2 = recip(2,1)*m1+recip(2,2)*m2+recip(2,3)*m3
             mhat3 = recip(3,1)*m1+recip(3,2)*m2+recip(3,3)*m3
             msq = mhat1*mhat1+mhat2*mhat2+mhat3*mhat3
             Phi(1,k1,ka2,ka3)=Phi(1,k1,ka2,ka3)*EXP(-fac*msq)
             Phi(2,k1,ka2,ka3)=Phi(2,k1,ka2,ka3)*EXP(-fac*msq)
          END DO
       END DO
    END DO
  END SUBROUTINE GetSmoothingGaussian
  SUBROUTINE InitializeRoutineVariables

    INTEGER :: iceil
    order=Aorder
    ewald_coeff=Aewald_coeff
    efact=Aefact
    unitc=Aunitc
    CFact=1.0D0/DSQRT(unitc)
    SigmaIon=ASigmaIon
    ChargeIon=AChargeIon
    SmoothFactor=ASmoothFactor
    DoFreeEnergy=ADoFreeEnergy
    rkcut=Arkcut
    nfft1=Anfft1
    nfft2=Anfft2
    nfft3=Anfft3
    nf1 = nfft1/2
    if ( 2*nf1 .lt. nfft1 )nf1 = nf1+1
    nf2 = nfft2/2
    if ( 2*nf2 .lt. nfft2 )nf2 = nf2+1
    nf3 = nfft3/2
    if ( 2*nf3 .lt. nfft3 )nf3 = nf3+1
    
    Bin1=2.0D0/nfft1
    Bin2=2.0D0/nfft2
    Bin3=2.0D0/nfft3
    NatoSlt=ANatoslt
    nprocs=Anprocs
    node=Anode
    ncube=Ancube
    nbyte=Anbyte
    rbyte=Arbyte
    ntap=Antap; ngrp=Angrp
    WRITE(*,*) ntap,NatoSlt
    nstart_h=Anstart_h; nend_h=Anend_h; nlocal_h=Anlocal_h
    nstart_ah=Anstart_ah; nend_ah=Anend_ah; nlocal_ah=Anlocal_ah
    nprot=Anprot
    npy=1
    npz=nprocs
    nodex=0
    nodey=0
    nodez=node

    kT=8.3143D0*300/1000.0D0
    nay=iceil(nfft2,npy)
    naz=iceil(nfft3,npz)
    kstart=nodez*naz+1
    kend  =(nodez+1)*naz
    jstart=nodey*nay+1
    jend  =(nodey+1)*nay
    ALLOCATE(nbtype(ntap))
    ALLOCATE(Sigma(NatoSlt))
    
    nbtype(:)=Anbtype(1:ntap)
    Sigma(1:NatoSlt)=Apnbd1(nbtype(1:NatoSlt))
    Sigma=Sigma*RadiusToSigma
    HardSphereCutoff=-1.0D9
    DO i=1,NatoSlt
       IF(0.5D0*(Sigma(i)+SigmaIon) > HardSphereCutoff)  &
            HardSphereCutoff=(Sigma(i)+SigmaIon)*0.5D0
    END DO
    
    
    ALLOCATE(grppt(2,ngrp),charge(ntap),ChargesWithIon(ntap))
    ALLOCATE(mass(ntap),pmass(ntap),ss_index(ntap),Label(ntap))
    ALLOCATE(xpa(ntap),ypa(ntap),zpa(ntap))
    ALLOCATE(XpAvg(NatoSlt),YpAvg(NatoSlt),ZpAvg(NatoSlt))
    ALLOCATE(fr1(ntap),fr2(ntap),fr3(ntap))
    ALLOCATE(xpg(ngrp),ypg(ngrp),zpg(ngrp))
    ALLOCATE(xpga(ngrp),ypga(ngrp),zpga(ngrp))
    ALLOCATE(xpcm(nprot),ypcm(nprot),zpcm(nprot))
    ALLOCATE(xpcma(nprot),ypcma(nprot),zpcma(nprot))
    ALLOCATE(bsp_mod1(2,nfft1),bsp_mod2(2,nfft2),bsp_mod3(2,nfft3))
    ALLOCATE(Phi(nfft1,nfft2,nfft3))
    ALLOCATE(Rho(nfft1,nfft2,nfft3))
    ALLOCATE(Rho0(nfft1,nfft2,nfft3))
    ALLOCATE(Phi0(nfft1,nfft2,nfft3))
    ALLOCATE(HardSphere(nfft1,nfft2,nfft3))

    RecordLength=nfft1*nfft2*10

    READ(AKElecPot,'(a80)',ERR=2000) FileNamePhiDat
    READ(AKElecPot,'(a80)',ERR=2000) FileNamePhiCoord
    READ(AKElecPot,'(a80)',ERR=2000) FileNamePhi


    CALL openf(UnitPhi,FileNamePhi,'U','U',RecordLength)
    CLOSE(UnitPhi)

    CALL openf(UnitPhiDat,FileNamePhidat,'F','U',0)
    WRITE(UnitPhiDat,*) nfft1,nfft2,nfft3,RecordLength,NatoSlt

    CALL openf(UnitPhiCoord,FileNamePhiCoord,'U','U',0)

    Phi=0.0D0
    Rho=0.0D0
    Phi0=0.0D0
    Rho0=0.0D0
    HardSphere=0.0D0
    XpAvg=0.0D0
    YpAvg=0.0D0
    ZpAvg=0.0D0
    bsp_mod1=0.0D0
    bsp_mod2=0.0D0
    bsp_mod3=0.0D0
    CALL LoadBspModuli
    grppt(:,:)=Agrppt(1:2,1:ngrp)
    Charge(:)=Acharge(1:ntap)
    ChargesWithIon(:)=Acharge(1:ntap)
    mass(:)=Amass(1:ntap)
    pmass(:)=Apmass(1:ntap)
    Label(:)=ALabel(1:ntap)
    ss_index(:)=Ass_index(1:ntap)
    count=0
    DO j=1,nprot
       m=Aprotl(1+count)
       count=count+m+1
    END DO
    ALLOCATE(protl(count+1))
    count=0
    DO j=1,nprot
       m=Aprotl(1+count)
       protl(1+count)=Aprotl(1+count)
       DO i=1,m
          protl(1+count+i)=Aprotl(1+count+i)
       END DO
       count=count+m+1
    END DO
    
!!$*=======================================================================
!!$*--- Allocate Memory for PME -------------------------------------------
!!$*=======================================================================
    
    ALLOCATE(theta1(order,ntap),theta2(order,ntap),theta3(order,ntap))
    ALLOCATE(dtheta1(order,ntap),dtheta2(order,ntap),dtheta3(order,ntap))
    ALLOCATE(indk1(order,ntap),indk2(order,ntap) &
         ,indj1(order,ntap),indj2(order,ntap))
    ALLOCATE(mk(order,ntap),mj(order,ntap))
    RETURN
2000 WRITE(*,'(''Error Reading Electrostatic Potential auxiliary file. Filename missing.'')') 
    STOP
    
  END SUBROUTINE InitializeRoutineVariables
END SUBROUTINE CompElecPotentialOnGrid
