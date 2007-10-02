
      SUBROUTINE calc_lda_flu(kout,ninst,nstep,fstep,nmol,natm,nzero
     &                       ,xp0,yp0,zp0)
*=======================================================================
      implicit none

      INTEGER kout,ninst,nstep,nmol,natm,nzero
      REAL*8 fstep,xp0(*),yp0(*),zp0(*)

c---- LOCAL VARIABLES

      REAL*8 trans(9),gauche(9),norma,rate(9),freq
      REAL*8 xr1,xr2,xr3,xr4,yr1,yr2,yr3,yr4,zr1,zr2,zr3,zr4,x21,x32
     &     ,x43,y21,y32,y43,z21,z32,z43,rsq21,rsq32,rsq43,rsp21,rsp32
     &     ,rsp43,cb1,cb2,cb3,sb1,sb2,sb3,aux,coa,quasi_zero
      INTEGER i,j,idx,first,sumt,idx2,nmax
      PARAMETER(nmax=1500)
      INTEGER*4 old_gauche(nmax)
      DATA quasi_zero/1.0D-10/
      DATA first/1/

*=======================================================================

      IF ( first .EQ. 1) THEN
          first = first + 1
          sumt = 0
          freq=nstep/fstep
          CALL zero0(trans,9)
          CALL zero0(gauche,9)
          CALL zero0(rate,9)
          do i=1,nmax
             old_gauche(i) = 0
          end do
          IF(nmax .lt. nmol*9) THEN
            write(*,*) 
            write(*,*) 'WARNING: calc_lda_flu.f, increase nmax'
            write(*,*) 
            STOP
          END IF
      END IF

*=======================================================================

      idx = nzero+11
c      nmol=1
      DO i=1,nmol
         do j=1,9
             xr1=xp0(idx)
             yr1=yp0(idx)
             zr1=zp0(idx)
             xr2=xp0(idx+3)
             yr2=yp0(idx+3)
             zr2=zp0(idx+3)
             xr3=xp0(idx+6)
             yr3=yp0(idx+6)
             zr3=zp0(idx+6)
             xr4=xp0(idx+9)
             yr4=yp0(idx+9)
             zr4=zp0(idx+9)
             x21=xr2-xr1
             y21=yr2-yr1
             z21=zr2-zr1
             x32=xr3-xr2
             y32=yr3-yr2
             z32=zr3-zr2
             x43=xr4-xr3
             y43=yr4-yr3
             z43=zr4-zr3
             rsq21=x21**2+y21**2+z21**2
             rsq32=x32**2+y32**2+z32**2
             rsq43=x43**2+y43**2+z43**2
             rsp21=DSQRT(rsq21)
             rsp32=DSQRT(rsq32)
             rsp43=DSQRT(rsq43)
             cb1=(x21*x32+y21*y32+z21*z32)/(rsp21*rsp32)
             cb2=(x43*x32+y43*y32+z43*z32)/(rsp43*rsp32)
             cb3=(x21*x43+y21*y43+z21*z43)/(rsp21*rsp43)
             sb1=DSQRT(DABS(1.0d0-cb1**2))
             sb2=DSQRT(DABS(1.0d0-cb2**2))
             sb3=DSQRT(DABS(1.0d0-cb3**2))
             IF(sb1 .EQ. 0.0D0) sb1=quasi_zero
             IF(sb2 .EQ. 0.0D0) sb2=quasi_zero
             IF(sb3 .EQ. 0.0D0) sb3=quasi_zero
             aux=sb1*sb2
             coa=(cb1*cb2-cb3)/aux

             idx2=(i-1)*9+j
             IF(coa .GE. 0.0d0 ) THEN
                gauche(j) = gauche(j) + 1.0d0

                if ( old_gauche(idx2) .EQ. 0 ) rate(j)=rate(j)+1.0d0
                old_gauche(idx2) = 1

             ELSE

                trans(j) = trans(j) + 1.0d0
                old_gauche(idx2) = 0

             END IF

             idx = idx + 3
          end do

          idx = idx + natm
      END DO

      sumt = sumt + 1
c      write(80,*) sumt,old_gauche(2),rate(2),coa

*=======================================================================
            
      IF( MOD(nstep,ninst) .EQ. 0) THEN
          norma = 1.0d0/(sumt*nmol)
          WRITE(kout,'('' Time/fs = '',f12.2)') (sumt/freq)
          WRITE(kout,'('' Dihedral   P-trans       P-gauche       k-rate
     &-ns'')')
          DO j=1,9
c             write(kout,*) rate(j),norma,freq
             WRITE(kout,100) j, trans(j)*norma,gauche(j)*norma
     &                       ,(rate(j)*norma*freq*1.0d+6)
          END DO
          WRITE(kout,*)
      END IF

*=======================================================================

      RETURN
100   FORMAT(2x,i4,2(2x,f12.6),2x,f12.3)
      END
