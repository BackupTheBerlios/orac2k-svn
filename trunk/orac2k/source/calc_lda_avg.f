      SUBROUTINE calc_lda_avg(krmin,keend,
     &     nlda_mol,nlda_atm,nlda_zero,wca,xp_ini,
     &     yp_ini,zp_ini,qt,xp0,yp0,zp0,nato,nstep,ninst_lda,
     &     fstep)

************************************************************************
*                                                                      *
*              Author:  Matteo Ceccarelli                              *
*              CECAM/ENS Lyon, FRANCE                                  *
*                                                                      *
*              - May 1999 -                                            *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none
      INCLUDE  'parst.h'

      integer maxnpoint
      parameter(maxnpoint=100)

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8  wca(*),xp0(*),yp0(*),zp0(*),xp_ini(*),yp_ini(*),zp_ini(*)
      real*8  qt(4),fstep
      INTEGER nlda_mol,nlda_atm,nlda_zero,nato,nstep,ninst_lda
      INTEGER krmin,keend

*----------------------- VARIABLES IN COMMON --------------------------*

      REAL*8  xyz(3,m1),xyz0(3,m1),xyzfit(3,m1),wca2(m1),work(m1)
      COMMON /rag1/ xyz,xyz0,xyzfit,work,wca2

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8 err,dist,diff,dpoint,fact,sum,tpoint
      REAL*8 a1,a1_min
      REAL*8 dd,dd_min(12),gr_c_prot(12,maxnpoint)
      REAL*8 vcp(12,maxnpoint)
      REAL*8 gr_a1(maxnpoint)
      real*8 gr_cc(maxnpoint),v4(maxnpoint),eend
      real*8 gr_tt(maxnpoint),v5(maxnpoint),zend
      real*8 gr_hm(maxnpoint),v6(maxnpoint)
      real*8 gr_xy(maxnpoint),v7(maxnpoint),xyend,xend,yend
      REAL*8 va1(maxnpoint)
      integer cdx
      INTEGER n,i,j,k,idx_rc,idx,npoint,idx_min,iter_lda,idx2,bin,nfreq
      DATA iter_lda/0/

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      if( iter_lda .eq. 0 ) then
          CALL zero0(gr_c_prot,12*maxnpoint)
          CALL zero0(gr_a1,maxnpoint)
          CALL zero0(gr_xy,maxnpoint)
          CALL zeroa(gr_tt,gr_cc,gr_hm,maxnpoint,1)
      end if

      nfreq = 2
      npoint = 100
      dpoint = 0.20d0
      tpoint = 180.0d0/npoint

      sum=0.0D0
      DO i=1,nato
         xyz(1,i)=xp0(i)
         xyz(2,i)=yp0(i)
         xyz(3,i)=zp0(i)
         xyz0(1,i)=xp_ini(i)
         xyz0(2,i)=yp_ini(i)
         xyz0(3,i)=zp_ini(i)
         wca2(i)=wca(i)
         sum = sum + wca2(i)
      END DO

      CALL normal(wca2,nato)

      IF(DABS(sum) .GT. 1.0D-3) THEN
         CALL xfit(xyz0,xyz,xyzfit,qt,wca2,work,nato,err)
        
 
         DO k=1,nlda_mol
            a1_min  = 100.0d0
            do cdx = 1,12
               dd_min(cdx) = 100.0d0
            end do
            idx = nlda_zero + (k-1)*nlda_atm

c
c        RMIN
c
            do i=1,nlda_zero
               a1  = 0.0d0

               do j=1,3
                  a1 = a1 + (xyzfit(j,idx+1)-xyzfit(j,i))**2
               end do
               a1 = sqrt(a1)
               if ( a1 .lt. a1_min ) a1_min = a1

               
               do cdx = 1,12
                  dd = 0.0d0
                  do j=1,3
            dd = dd + (xyzfit(j,idx+11+3*(cdx-1))-xyzfit(j,i))**2
                  end do
                  dd = dsqrt(dd)
                  if(dd .lt. dd_min(cdx)) then
                     dd_min(cdx) = dd
                  end if
               end do

            end do

            bin = int(a1_min/dpoint) +1
            if ( bin .le. npoint ) then
               gr_a1(bin) = gr_a1(bin) + 1.0
            end if

            do cdx=1,12
               bin = int(dd_min(cdx)/dpoint) +1
               if ( bin .le. npoint ) then
                  gr_c_prot(cdx,bin) = gr_c_prot(cdx,bin) + 1.0
               end if
            end do

c
c        C1-C12 and (C1-C12)_Z
c
            xend = xyzfit(1,idx+11)-xyzfit(1,idx+44)
            yend = xyzfit(2,idx+11)-xyzfit(2,idx+44)
            zend = xyzfit(3,idx+11)-xyzfit(3,idx+44)
            xyend = xend**2 + yend**2
            eend = xyend + zend**2
            eend = dsqrt(eend)
            zend = acos(zend/eend)*180.0/3.14
            xyend = dsqrt(xyend)
            xyend = acos(xend/xyend)*180.0/3.14

            bin = int(eend/dpoint) + 1
            if ( bin .le. npoint ) then
               gr_cc(bin) = gr_cc(bin) + 1.0
            end if

            bin = int(zend/tpoint) + 1
            if ( bin .le. npoint ) then
               gr_tt(bin) = gr_tt(bin) + 1.0
            end if

            bin = int(xyend/tpoint) + 1
            if ( bin .le. npoint ) then
               gr_xy(bin) = gr_xy(bin) + 1.0
            end if


            do i=k+1,nlda_mol
               idx2 = nlda_zero + (i-1)*nlda_atm
               dist = 0.0d0
               do j=1,3
                  dist = dist + (xyzfit(j,idx+1)-xyzfit(j,idx2+1))**2
               end do
               dist=dsqrt(dist)
               bin = int(dist/dpoint) + 1
               if ( bin .le. npoint ) then
                  gr_hm(bin) = gr_hm(bin) + 2.0
               end if
            end do

         END DO

         iter_lda = iter_lda + 1

         IF ( mod(iter_lda,nfreq).EQ.0) THEN
            CALL dcopy(npoint,gr_a1,1,va1,1)
            CALL dcopy(npoint*12,gr_c_prot,1,vcp,1)
            CALL dcopy(npoint,gr_cc,1,v4,1)
            CALL dcopy(npoint,gr_tt,1,v5,1)
            CALL dcopy(npoint,gr_hm,1,v6,1)
            CALL dcopy(npoint,gr_xy,1,v7,1)
            fact=1.0D0/DBLE(iter_lda*nlda_mol)
            CALL dscal(npoint,fact,va1,1)
            CALL dscal(npoint*12,fact,vcp,1)
            CALL dscal(npoint,fact,v4,1)
            CALL dscal(npoint,fact,v5,1)
            CALL dscal(npoint,fact,v6,1)
            CALL dscal(npoint,fact,v7,1)

c            REWIND krmin
c            REWIND keend
            write(krmin,'(a7,2x,f12.3)') 'Tstep =',fstep 
            write(keend,'(a7,2x,f12.3)') 'Tstep =',fstep 
            do j=1,npoint
                write(krmin,'(f5.2,13(1x,f5.3))') j*dpoint,va1(j)
     &          ,(vcp(cdx,j),cdx=1,12)
                write(keend,'(6(2x,f7.3))') j*dpoint,v4(j)
     &          ,v6(j)/(j*dpoint**2),j*tpoint,v5(j),v7(j)
            end do
            iter_lda = 0.0d0
          END IF

      ELSE
         write(krmin,*) 'No atoms to fit'
      END IF

         
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END

