      SUBROUTINE Get_virtual_selfenergy(xp0,yp0,zp0,virtual_atoms,ecc6
     &     ,ecc12,nbtype,type,ma,virtual_energy)
      IMPLICIT none
      REAL*8 xp0(*),yp0(*),zp0(*),ecc6(*),ecc12(*),virtual_energy
      INTEGER ma,virtual_atoms(*)
      INTEGER nbtype(*),type(ma,ma)

      INTEGER i,j,typei,typej,ii,jj,m,lij
      REAL*8 conf,xpi,ypi,zpi,rsq,xc,yc,zc,rsqi,r6,r12

      m=virtual_atoms(1)
      virtual_energy=0.0D0
      DO ii=1,m
         i=virtual_atoms(1+ii)
         xpi=xp0(i)
         ypi=yp0(i)
         zpi=zp0(i)
         DO jj=ii+1,m
            j=virtual_atoms(1+jj)
            lij=type(nbtype(i),nbtype(j))
            xc=xpi-xp0(j)
            yc=ypi-yp0(j)
            zc=zpi-zp0(j)
            rsq=xc*xc+yc*yc+zc*zc
            rsqi=1.0D0/rsq
            r6=rsqi*rsqi*rsqi
            r12=r6*r6
            conf=ecc12(lij)*r12-ecc6(lij)*r6
            virtual_energy=virtual_energy+conf
         END DO
      END DO
      virtual_energy=-virtual_energy
      RETURN
      END
