      SUBROUTINE write_electric_field(xp0,yp0,zp0,xpa,ypa,zpa,xpt,ypt  
     &           ,zpt,wca,fstep,fpx,fpy,fpz,Edx,Edy,Edz,dipx          
     &           ,dipy,dipz,ene,enedip,ucoul,Udd,Ued,Uind,node)

************************************************************************
*                                                                      *
*              Author:  Matteo Ceccarelli                              *
*              CECAM/ENS Lyon, FRANCE, may 99                          *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


!======================== DECLARATIONS ================================*

      IMPLICIT none

!----------------------------- ARGUMENTS ------------------------------*

      REAL*8  xpa(*),ypa(*),zpa(*),xp0(*),yp0(*),zp0(*)
      REAL*8  xpt(*),ypt(*),zpt(*),fpx(*),fpy(*),fpz(*),ene(*),enedip(*)
      REAL*8  wca(*),Edx(*),Edy(*),Edz(*),dipx(*),dipy(*),dipz(*)
      REAL*8  fstep,ucoul,Udd,Ued,Uind
      INTEGER node
!----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'parst.h'
      INCLUDE 'cpropar.h'
      INCLUDE 'unit.h'

      REAL*8  xyz(3,m1),xyz0(3,m1),xyzfit(3,m1),wca2(m1),work(m1),qt(4)

      COMMON /rag1/ xyz,xyz0,xyzfit,work,wca2

!------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8 sum,rot(3,3)
      REAL*8 xpi,ypi,zpi,furpar,chrgi,r_x,r_y,r_z,rr,error
      REAL*8 fnx(m1),fny(m1),fnz(m1)
      INTEGER aux1
      INTEGER i,j,k,i1,j1
      real*8 crij(3),r3ij,r5ij,arg,fact,cexp,TDi(3),TDj(3),Tij(3,3)
      real*8 auxi,auxj,auxd
      integer lij,l,m,nbti

      if(node.ne.0) RETURN

!=======================================================================
!     get template coordinate system by xfit 
!=======================================================================

      do i=1,3
         do j=1,3
            rot(i,j)=0.0d0
         end do
         rot(i,i)=1.0d0
      end do

      sum=0.0D0

      DO i=1,ntap
         xyz(1,i)=xp0(i)
         xyz(2,i)=yp0(i)
         xyz(3,i)=zp0(i)
         xyz0(1,i)=xpt(i)
         xyz0(2,i)=ypt(i)
         xyz0(3,i)=zpt(i)
         wca2(i)=wca(i)    
         sum = sum + wca2(i)
      END DO
      CALL normal(wca2,ntap)
      IF(DABS(sum) .GT. 1.0D-3) THEN
         CALL xfit(xyz0,xyz,xyzfit,qt,wca2,work,ntap,error)
   
!---- ROTATION MATRIX IN TERMS OF QUATERNIONS: ----------------------*

         rot(1,1)=-2.d0*qt(3)**2-2.d0*qt(4)**2+1.d0
         rot(1,2)=2.d0*(-qt(1)*qt(4)+qt(2)*qt(3))
         rot(1,3)=2.d0*(qt(1)*qt(3)+qt(2)*qt(4))
         rot(2,1)=2.d0*(qt(1)*qt(4)+qt(2)*qt(3))
         rot(2,2)=-2.d0*qt(2)**2-2.d0*qt(4)**2+1.d0
         rot(2,3)=2.d0*(-qt(1)*qt(2)+qt(3)*qt(4))
         rot(3,1)=2.d0*(-qt(1)*qt(3)+qt(2)*qt(4))
         rot(3,2)=2.d0*(qt(1)*qt(2)+qt(3)*qt(4))
         rot(3,3)=-2.d0*qt(2)**2-2.d0*qt(3)**2+1.d0
      END IF

      write(kout_ef,'(a3,2x,f12.2)') ' T ', fstep

      do k=1,nmol_ef
         aux1 = pmol_ef(1,k)
   
         write(kout_ef,'(a5,2x,i6,i6)') ' CRO ', k,aux1
         do i1=1,aux1
            i=pmol_ef(i1+1,k)
            chrgi=chrge(i)
            nbti=nbtype(i)
            xpi=xp0(i)
            ypi=yp0(i)
            zpi=zp0(i)
            do j1=i1+1,aux1
               j=pmol_ef(j1+1,k)
               r_x = xpi-xp0(j)
               r_y = ypi-yp0(j)
               r_z = zpi-zp0(j)
               rr  = r_x*r_x+r_y*r_y+r_z*r_z
               rr  = dsqrt(rr)
               furpar = chrgi*chrge(j)/rr
   
               if(polar .and. rr .le. rspoff) then
   
                  lij=type_table(nbti,nbtype(j))
                  arg = rr/plrzbij(lij)
                  fact = 1.0d0 + arg + 0.50d0*arg**2
                  cexp = exp(-arg)

                  crij(1) = r_x        
                  crij(2) = r_y        
                  crij(3) = r_z
                  r3ij = rr**3
                  r5ij = rr**5        
   
                  do m=1,3
                     TDi(m) = 0.0d0
                     TDj(m) = 0.0d0
                     do l=1,3
                        Tij(m,l) = 3.0d0*crij(m)*crij(l)/r5ij
                  Tij(m,l) = Tij(m,l)*(1.0d0-(fact+arg**3/6.0d0)*cexp)
                     end do
                     Tij(m,m) = Tij(m,m) - (1.0d0-fact*cexp)/r3ij

                     TDj(m) = TDj(m) + Tij(m,1)*dipx(j)
                     TDj(m) = TDj(m) + Tij(m,2)*dipy(j)
                     TDj(m) = TDj(m) + Tij(m,3)*dipz(j)
                     TDi(m) = TDi(m) + Tij(m,1)*dipx(i)
                     TDi(m) = TDi(m) + Tij(m,2)*dipy(i)
                     TDi(m) = TDi(m) + Tij(m,3)*dipz(i)
   
                  end do

                  Edx(i) = Edx(i) - TDj(1)
                  Edy(i) = Edy(i) - TDj(2)
                  Edz(i) = Edz(i) - TDj(3)
                  Edx(j) = Edx(j) - TDi(1)
                  Edy(j) = Edy(j) - TDi(2)
                  Edz(j) = Edz(j) - TDi(3)
   
           auxd = dipx(i)*TDj(1)+dipy(i)*TDj(2)+dipz(i)*TDj(3)
           Udd  = Udd  + auxd
   
                  auxi = dipx(j)*r_x+dipy(j)*r_y+dipz(j)*r_z
                  auxj = dipx(i)*r_x+dipy(i)*r_y+dipz(i)*r_z
                  enedip(i) = enedip(i) - auxi/r3ij
                  enedip(j) = enedip(j) + auxj/r3ij
   
               end if
!----- ENE
               ucoul  = ucoul  - furpar
               ene(i) = ene(i) - furpar
               ene(j) = ene(j) - furpar
!----- FORCES
               furpar = furpar / (rr*rr)
               fpx(i) = fpx(i) - furpar*r_x
               fpy(i) = fpy(i) - furpar*r_y
               fpz(i) = fpz(i) - furpar*r_z
               fpx(j) = fpx(j) + furpar*r_x
               fpy(j) = fpy(j) + furpar*r_y
               fpz(j) = fpz(j) + furpar*r_z
            end do
         END DO

!   ELECTRIC POTENTIAL  unit:  ElectPot=Hartree/e-
!   ELECTRIC FIELD      unit:  Elec=Hartree/e-/ao
         DO i1=1,aux1
            i=pmol_ef(i1+1,k)
            chrgi=chrge(i)
            fpx(i) = (Edx(i) + fpx(i)/chrgi) * unitefield
            fpy(i) = (Edy(i) + fpy(i)/chrgi) * unitefield
            fpz(i) = (Edz(i) + fpz(i)/chrgi) * unitefield
   
            fnx(i) = rot(1,1)*fpx(i) + rot(1,2)*fpy(i) + rot(1,3)*fpz(i)
            fny(i) = rot(2,1)*fpx(i) + rot(2,2)*fpy(i) + rot(2,3)*fpz(i)
            fnz(i) = rot(3,1)*fpx(i) + rot(3,2)*fpy(i) + rot(3,3)*fpz(i)

            ene(i) = ene(i) - alphal/DSQRT(pi)*chrgi**2
            ene(i) = (enedip(i) + ene(i)/chrgi) * unitepot

            write(kout_ef,101) 
     &            beta(i)(1:4),i,fnx(i),fny(i),fnz(i),ene(i)
         END DO
      end do

      RETURN
101   FORMAT(a4,1x,i7,4(1x,e15.9))
      END
