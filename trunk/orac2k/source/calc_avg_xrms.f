      SUBROUTINE calc_avg_xrms(avg_ca,avg_he,avg_bc,fstep,kout,xp_ini
     &     ,yp_ini,zp_ini,xp0,yp0,zp0,wca,whe,wbc,protl
     &     ,anprot,annpro,anpoint,ntap)

************************************************************************
*   Time-stamp: <97/09/12 12:40:07 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Apr 20 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none
*----------------------- ARGUMENTS -------------------------------------

      INTEGER ntap,kout,protl(*),annpro,anpoint(2,*)
      REAL*8  xp0(*),yp0(*),zp0(*),xp_ini(*),yp_ini(*),zp_ini(*),wca(*)
     &     ,whe(*),wbc(*),fstep
      LOGICAL avg_ca,avg_he,avg_bc,anprot

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,n,m,k,nprot,count,atoms,initial,n_ca,n_he,n_bc
      REAL*8  xb,yb,zb,dx,dy,dz,dr,sum_ca,sum_he,sum_bc
      DATA  initial/0/

*==================== EXECUTABLE STATEMENTS ============================

      IF(initial .EQ. 0) THEN
         IF(avg_ca) WRITE(kout,3) 
         IF(avg_he) WRITE(kout,4) 
         IF(avg_bc) WRITE(kout,5) 
         WRITE(kout,1) 
         initial=initial+1
      END IF


      IF ( anprot ) THEN

      do n = 1,annpro
         sum_ca=0.0D0
         sum_he=0.0D0
         sum_bc=0.0D0
         n_ca=0
         n_he=0
         n_bc=0
         m=anpoint(2,n)-anpoint(1,n)
         DO k=1,m
            i=anpoint(1,n)-1+k
            xb=xp0(i)
            yb=yp0(i)
            zb=zp0(i)
            dx=xb-xp_ini(i)
            dy=yb-yp_ini(i)
            dz=zb-zp_ini(i)
            dr=DSQRT(dx**2+dy**2+dz**2)
            IF(avg_ca .AND. wca(i) .GT. 1.0D-6) THEN
               sum_ca=sum_ca+dr
               n_ca=n_ca+1
            END IF
            IF(avg_he .AND. whe(i) .GT. 1.0D-6) THEN
               sum_he=sum_he+dr
               n_he=n_he+1
            END IF
            IF(avg_bc .AND. wbc(i) .GT. 1.0D-6) THEN
               sum_bc=sum_bc+dr
               n_bc=n_bc+1
            END IF
         END DO

         IF(n_ca .EQ. 0) THEN
            sum_ca=0.0D0
         ELSE
            sum_ca=sum_ca/DBLE(n_ca)
         END IF
         IF(n_he .EQ. 0) THEN
            sum_he=0.0D0
         ELSE
            sum_he=sum_he/DBLE(n_he)
         END IF
         IF(n_bc .EQ. 0) THEN
            sum_bc=0.0D0
         ELSE 
            sum_bc=sum_bc/DBLE(n_bc)
         END IF
         IF(DABS(sum_he+sum_bc+sum_ca) .GT. 1.0D-6) THEN

            WRITE(kout,2) fstep,n,sum_ca,sum_he,sum_bc

         END IF
      END DO

      ELSE

      count=0
      atoms=0
      nprot=0
      DO WHILE(atoms .LT. ntap)
         m=protl(1+count)
         sum_ca=0.0D0
         sum_he=0.0D0
         sum_bc=0.0D0
         n_ca=0
         n_he=0
         n_bc=0
         DO k=1,m
            i=protl(1+count+k)
            xb=xp0(i)
            yb=yp0(i)
            zb=zp0(i)
            dx=xb-xp_ini(i)
            dy=yb-yp_ini(i)
            dz=zb-zp_ini(i)
            dr=DSQRT(dx**2+dy**2+dz**2)
            IF(avg_ca .AND. wca(i) .GT. 1.0D-6) THEN
               sum_ca=sum_ca+dr
               n_ca=n_ca+1
            END IF
            IF(avg_he .AND. whe(i) .GT. 1.0D-6) THEN
               sum_he=sum_he+dr
               n_he=n_he+1
            END IF
            IF(avg_bc .AND. wbc(i) .GT. 1.0D-6) THEN
               sum_bc=sum_bc+dr
               n_bc=n_bc+1
            END IF
         END DO
         IF(n_ca .EQ. 0) THEN
            sum_ca=0.0D0
         ELSE
            sum_ca=sum_ca/DBLE(n_ca)
         END IF
         IF(n_he .EQ. 0) THEN
            sum_he=0.0D0
         ELSE
            sum_he=sum_he/DBLE(n_he)
         END IF
         IF(n_bc .EQ. 0) THEN
            sum_bc=0.0D0
         ELSE 
            sum_bc=sum_bc/DBLE(n_bc)
         END IF
         nprot=nprot+1
         count=count+m+1
         atoms=atoms+m
c
c    modified by matteo   24/04/98
c
         IF(DABS(sum_he+sum_bc+sum_ca) .GT. 1.0D-6) THEN

            WRITE(kout,2) fstep,nprot,sum_ca,sum_he,sum_bc

         END IF
      END DO

      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

1     FORMAT('    Tstep     Solute No.    CA           HE',
     &'             BC'/)
2     FORMAT(f12.3,2x,i6,2x,f12.5,2x,f12.5,2x,f12.5)
3     FORMAT('***   Rigid body fit on CA atoms       ***')
4     FORMAT('***   Rigid body fit on heavy atoms    ***')
5     FORMAT('***   Rigid body fit on backbone atoms ***')

      RETURN
      END
