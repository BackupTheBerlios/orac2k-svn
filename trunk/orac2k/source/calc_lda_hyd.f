
      SUBROUTINE calc_lda_hyd(kout,ninst,nstep,fstep,nmol,natm,nzero,co
     &                       ,xpa,ypa,zpa,nato_slt,nato_slv,nmol_slv)
*=======================================================================
      implicit none
 

      INTEGER kout,ninst,nstep,nmol,natm,nzero
      INTEGER nmol_slv,nato_slv,nato_slt
      REAL*8 fstep,xpa(*),ypa(*),zpa(*),co(3,3)

c---- LOCAL VARIABLES

      REAL*8 nhyd(13),norma,maxdist,dist_NO,dist_CO,dist
      REAL*8 xv,yv,zv,xl,yl,zl,xg,yg,zg,xc,yc,zc,xmap0,ymap0,zmap0
      INTEGER i,j,k,lv,ll,idx,first,sumt,idx_slv
      DATA first/1/

      INCLUDE 'pbc.h'

*=======================================================================

      IF ( first .EQ. 1) THEN
         first = first + 1
c------- O atom in water, number=1
         idx_slv = nato_slt + 1

c------- N atom in lda, number=1
         idx = nzero + 1
c
         maxdist = 25.0d0
         dist_NO = 5.0d0
         dist_CO = 3.50d0
         sumt = 0
         CALL zero0(nhyd,13)
      END IF

*=======================================================================


      DO i=1,nmol_slv
         lv = idx_slv + (i-1)*nato_slv
         xv = xpa(lv)
         yv = ypa(lv)
         zv = zpa(lv)

         DO j=1,nmol

c---------- first C atom in the tail, number=11
            ll = idx + (j-1)*natm + 10
            xl = xpa(ll)
            yl = ypa(ll)
            zl = zpa(ll)

            xg = xv - xl
            yg = yv - yl
            zg = zv - zl

            xmap0 = 2.0d0*PBC(xg)
            ymap0 = 2.0d0*PBC(yg)
            zmap0 = 2.0d0*PBC(zg)
 
            xg = xg - xmap0
            yg = yg - ymap0
            zg = zg - zmap0

            xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
            yc=           co(2,2)*yg+co(2,3)*zg
            zc=                      co(3,3)*zg

            dist = xc**2+yc**2+zc**2
            dist = DSQRT(dist)
 
            IF( dist .GE. maxdist ) GO TO 20

            IF( dist. LE. dist_CO ) THEN
               nhyd(2) = nhyd(2) + 1.0d0
c               write(80,*) i,j,'C1',dist
            END IF

            DO k=2,12

c---------- next C atom in the tail, + 3
               ll = ll + 3
               xl = xpa(ll)
               yl = ypa(ll)
               zl = zpa(ll)

               xg = xv - xl
               yg = yv - yl
               zg = zv - zl

               xmap0 = 2.0d0*PBC(xg)
               ymap0 = 2.0d0*PBC(yg)
               zmap0 = 2.0d0*PBC(zg)
 
               xg = xg - xmap0
               yg = yg - ymap0
               zg = zg - zmap0

               xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
               yc=           co(2,2)*yg+co(2,3)*zg
               zc=                      co(3,3)*zg

               dist = xc**2+yc**2+zc**2

               dist = DSQRT(dist)

               IF( dist. LE. dist_CO ) THEN
                  nhyd(1+k) = nhyd(1+k) + 1.0d0
c               write(82,*) i,j,k,'Ck',dist
               END IF

 30         END DO

c---------- now the head
            ll = idx + (j-1)*natm
            xl = xpa(ll)
            yl = ypa(ll)
            zl = zpa(ll)

            xg = xv - xl
            yg = yv - yl
            zg = zv - zl

            xmap0 = 2.0d0*PBC(xg)
            ymap0 = 2.0d0*PBC(yg)
            zmap0 = 2.0d0*PBC(zg)
 
            xg = xg - xmap0
            yg = yg - ymap0
            zg = zg - zmap0

            xc=co(1,1)*xg+co(1,2)*yg+co(1,3)*zg
            yc=           co(2,2)*yg+co(2,3)*zg
            zc=                      co(3,3)*zg

            dist = xc**2+yc**2+zc**2
            dist = DSQRT(dist)

            IF( dist. LE. dist_NO ) THEN
               nhyd(1) = nhyd(1) + 1.0d0
c               write(81,*) i,j,'N',dist
            END IF

 20      END DO

 10   END DO

c      write(80,*)
c      write(81,*)
c      write(82,*)

      sumt = sumt + 1

*=======================================================================
            
      IF( MOD(nstep,ninst) .EQ. 0) THEN
          norma = 1.0d0/(sumt*nmol)
          WRITE(kout,'('' Tstep = '',f12.2)') fstep
          WRITE(kout,'('' Atom      Hydration number'')')
          DO j=1,13
             WRITE(kout,100) j, nhyd(j)*norma
          END DO
          WRITE(kout,*)
      END IF

*=======================================================================

      RETURN
100   FORMAT(2x,i4,2x,f12.6)
      END
