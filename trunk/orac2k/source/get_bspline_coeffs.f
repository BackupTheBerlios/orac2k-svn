      SUBROUTINE get_bspline_coeffs(numatoms,fr1,fr2,fr3,order,theta1
     &     ,theta2,theta3,dtheta1,dtheta2,dtheta3,nfft1,nfft2,nfft3
     &     ,kstart,kend,jstart,jend,indk1,indk2,indj1,indj2,mk,mj)

************************************************************************
*   Time-stamp: <2009-03-12 12:28:10 marchi>                             *
*                                                                      *
c---------------------------------------------------------------------
c INPUT:
c      numatoms: number of atoms
c      fr1,fr2,fr3 the scaled and shifted fractional coords
c      order: the order of spline interpolation
c OUTPUT
c      theta1,theta2,theta3: the spline coeff arrays
c      dtheta1,dtheta2,dtheta3: the 1st deriv of spline coeff arrays
c---------------------------------------------------------------------
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author: Tom Darden NIH                                  *
*              Modified by: Massimo Marchi                             *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Feb  9 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER numatoms,order,nfft1,nfft2,nfft3,kstart,kend
     &     ,jstart,jend
      INTEGER indk1(order,*),indk2(order,*),indj1(order,*),indj2(order
     &     ,*),mk(*),mj(*)
      REAL*8  fr1(*),fr2(*),fr3(*)
      REAL*8  theta1(order,*),theta2(order,*),theta3(order,*)
     &     ,dtheta1(order,*),dtheta2(order,*),dtheta3(order,*)

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8  w
      INTEGER n,k0,j0,k,j,ith2,ith3,mjj,mkk

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO n = 1,numatoms
         k0 = int(fr3(n)) - order
         mkk=0
         DO ith3 = 1,order
            k0 = k0 + 1
            k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
            IF(k .GE. kstart .AND. k .LE. kend) THEN
               mkk=mkk+1
               indk1(mkk,n)=ith3
               indk2(mkk,n)=k
            END IF
         END DO
         mk(n)=mkk
         IF(mkk .NE. 0) THEN
            j0 = int(fr2(n)) - order
            mjj=0
            DO ith2 = 1,order
               j0 = j0 + 1
               j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2
               IF(j .GE. jstart .AND. j .LE. jend) THEN
                  mjj=mjj+1
                  indj1(mjj,n)=ith2
                  indj2(mjj,n)=j
               END IF
            END DO
            mj(n)=mjj
            IF(mjj .NE. 0) THEN
               w = fr1(n)-int(fr1(n))
               call fill_bspline(w,order,theta1(1,n),dtheta1(1,n))
               w = fr2(n)-int(fr2(n))
               call fill_bspline(w,order,theta2(1,n),dtheta2(1,n))
               w = fr3(n)-int(fr3(n))
               call fill_bspline(w,order,theta3(1,n),dtheta3(1,n))
            END IF
         END IF
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
