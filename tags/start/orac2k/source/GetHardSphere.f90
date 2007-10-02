SUBROUTINE GetHardSphere(co,TempHardSPhere)

!!$***********************************************************************
!!$   Time-stamp: <01/04/13 17:16:22 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Tue Apr  3 2001 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*


!!$======================== DECLARATIONS ================================*

  USE ElecPotential_Mod  
  IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

  REAL(8), DIMENSION (:,:,:) :: TempHardSphere
  REAL(8) :: co(3,3)

!!$----------------------- VARIABLES IN COMMON --------------------------*

!!$------------------------- LOCAL VARIABLES ----------------------------*

  REAL(8), DIMENSION (:), ALLOCATABLE :: Xnode,Ynode,Znode
  INTEGER, DIMENSION (:), ALLOCATABLE :: AtomInode
  INTEGER :: TotalAtoms,ii,i,j,k,iv,jv,kv,BorderRight,BorderLeft,k0
  REAL(8) :: dx,dy,dz,x1,y1,z1,x2,y2,z2,TempCutoff2,xc,yc,zc,Dist,xpn   &
       ,ypn,zpn
  INTEGER :: DeltaZ,nzPlus,nzMinus,nz,nind,nx,ny,nz0,n,m

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

  xpa=xpa+1.0D0
  ypa=ypa+1.0D0
  zpa=zpa+1.0D0
  nind=SIZE(IndX)
  dx=2.0D0/nfft1
  dy=2.0D0/nfft2
  dz=2.0D0/nfft3

  DeltaZ=2.0D0*SigmaIon*(oc(3,1)+oc(3,2)+oc(3,3))/dz
  BorderLeft=nodez*naz+1-DeltaZ
  BorderRight=nodez*naz+naz+DeltaZ
  TotalAtoms=0
  DO i=1,NatoSlt
     zc=zpa(i)/dz
     nz=INT(zc)+(SIGN(1.0D0,zc-INT(zc))-1)/2
     nzPlus=nz+nfft3
     nzMinus=nz-nfft3
     IF(  (BorderLeft < nz      .AND. nz      < BorderRight) .OR. &
          (BorderLeft < nzPlus  .AND. nzPlus  < BorderRight) .OR. &
          (BorderLeft < nzMinus .AND. nzMinus < BorderRight)) THEN
        TotalAtoms=TotalAtoms+1
     END IF
  END DO

  ALLOCATE(Xnode(TotalAtoms),Ynode(TotalAtoms),Znode(TotalAtoms))
  ALLOCATE(AtomInode(TotalAtoms))
  TotalAtoms=0
  DO i=1,NatoSlt
     zc=zpa(i)/dz
     nz=INT(zc)+(SIGN(1.0D0,zc-INT(zc))-1)/2
     nzPlus=nz+nfft3
     nzMinus=nz-nfft3
     IF(  (BorderLeft < nz      .AND. nz      < BorderRight) .OR. &
          (BorderLeft < nzPlus  .AND. nzPlus  < BorderRight) .OR. &
          (BorderLeft < nzMinus .AND. nzMinus < BorderRight)) THEN
        TotalAtoms=TotalAtoms+1
        Xnode(TotalAtoms)=xpa(i)
        Ynode(TotalAtoms)=ypa(i)
        Znode(TotalAtoms)=zpa(i)
        AtomInode(TotalAtoms)=i
     END IF
  END DO

  TempHardSphere=1.0D0
  DO ii=1,TotalAtoms
     n=AtomInode(ii)
     x1=xpa(n)/dx
     y1=ypa(n)/dy
     z1=zpa(n)/dz
     xpn=xpa(n)
     ypn=ypa(n)
     zpn=zpa(n)
     i=INT(x1)+(SIGN(1.0D0,x1-int(x1))-1)/2
     j=INT(y1)+(SIGN(1.0D0,y1-int(y1))-1)/2
     k=INT(z1)+(SIGN(1.0D0,z1-int(z1))-1)/2
     i=MOD(mod(i,nfft1)+nfft1,nfft1)
     j=MOD(mod(j,nfft2)+nfft2,nfft2)
     k=MOD(mod(k,nfft3)+nfft3,nfft3)
     TempCutoff2=0.25D0*(Sigma(n)+SigmaIon)**2
     DO m=1,nind
        iv=IndX(m)
        jv=IndY(m)
        kv=IndZ(m)
        nx=MOD(MOD(i+iv,nfft1)+nfft1,nfft1)
        ny=MOD(MOD(j+jv,nfft2)+nfft2,nfft2)
        nz=MOD(MOD(k+kv,nfft3)+nfft3,nfft3)
        nz0=nz-nodez*naz
        IF(nz0 >= 0 .AND. nz0 < naz) THEN
           x2=nx*dx
           y2=ny*dy
           z2=nz*dz
           x2=xpn-x2
           y2=ypn-y2
           z2=zpn-z2
           x2=x2-2.0D0*PBC(x2)
           y2=y2-2.0D0*PBC(y2)
           z2=z2-2.0D0*PBC(z2)
           xc=co(1,1)*x2+co(1,2)*y2+co(1,3)*z2
           yc=co(2,1)*x2+co(2,2)*y2+co(2,3)*z2
           zc=co(3,1)*x2+co(3,2)*y2+co(3,3)*z2
           Dist=xc*xc+yc*yc+zc*zc
           IF(Dist <= TempCutoff2) THEN
              TempHardSphere(nx+1,ny+1,nz0+1)=0.0D0
           END IF
        END IF
     END DO
  END DO
  xpa=xpa-1.0D0
  ypa=ypa-1.0D0
  zpa=zpa-1.0D0
  DEALLOCATE(Xnode,Ynode,Znode)
  DEALLOCATE(AtomInode)
  
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

CONTAINS
  FUNCTION PBC(x)
    IMPLICIT NONE 
    REAL(8) :: x,PBC

    PBC=DNINT(0.5D0*x)
  END FUNCTION PBC
END SUBROUTINE GetHardSphere
