!!$***********************************************************************
!!$   Time-stamp: <2005-03-06 22:23:36 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Feb 25 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This MODULE is part of the program ORAC ----*
MODULE Module_Extra_Forces
  TYPE extra_stretch
     INTEGER :: ind_a,ind_b
     REAL(8) :: K, r0
  END TYPE extra_stretch
  TYPE(extra_stretch), DIMENSION(:), POINTER, SAVE :: pop1, pop2, pop3
  TYPE extra_input
     CHARACTER(7) :: type_a,type_b
     REAL(8) :: K=0.0D0,r0=0.0D0
  END TYPE extra_input
  TYPE(extra_input), SAVE :: inp_1,inp_2,inp_3
  LOGICAL, SAVE :: Extra_Force=.FALSE.
CONTAINS

  SUBROUTINE Stretch_Init(beta,nprot,protl,nato,inp,pop)
!!$======================== DECLARATIONS ================================*
    IMPLICIT none
!!$----------------------------- ARGUMENTS ------------------------------*
    INTEGER :: nato,nprot,protl(*)
    CHARACTER(7) :: beta(*)
    TYPE(extra_stretch), DIMENSION(:), POINTER :: pop
    TYPE(extra_input) :: inp
!!$------------------------- LOCAL VARIABLES ----------------------------*

    INCLUDE 'unit.h'
    INTEGER :: i,count,p,n,ii,m,ip_a,ip_b
    LOGICAL :: ok_a,ok_b
    REAL(8) :: K,r0
    CHARACTER(7) :: type1,type2

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    K= inp % K
    r0= inp % r0
    type1 = inp % type_a
    type2 = inp % type_b
    count=0
    p=0
    DO n=1,nprot
       m=protl(count+1)
       ok_a=.FALSE.
       ok_b=.FALSE.
       ip_a=0
       ip_b=0
       DO ii=1,m
          i=protl(count+1+ii)
          IF(type1 == beta(i)) THEN
             ip_a=i
             ok_a=.TRUE.
             EXIT
          END IF
       END DO
       DO ii=1,m
          i=protl(count+1+ii)
          IF(type2 == beta(i)) THEN
             ip_b=i
             ok_b=.TRUE.
             EXIT
          END IF
       END DO
       IF(ok_a .AND. ok_b) THEN
          p=p+1
       ELSE IF(.NOT. ok_a .AND. .NOT. ok_b) THEN
          CONTINUE
       ELSE
          STOP
       END IF
       count=count+1+m
    END DO

    ALLOCATE(pop(p))
    p=0
    count=0
    DO n=1,nprot
       m=protl(count+1)
       DO ii=1,m
          i=protl(count+1+ii)
          IF(type1 == beta(i)) THEN
             p=p+1
             pop(p) % ind_a = i
             EXIT
          END IF
       END DO
       DO ii=1,m
          i=protl(count+1+ii)
          IF(type2 == beta(i)) THEN
             pop(p) % ind_b = i
             EXIT
          END IF
       END DO
       pop(p) % K = 1000.0D0*K*4.184/(unite*avogad)
       pop(p) % r0 = r0
       count=count+1+m
    END DO
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

  END SUBROUTINE Stretch_Init
  SUBROUTINE force(xp0,yp0,zp0,fpx,fpy,fpz,U_extra)
    IMPLICIT none
    
    REAL(8) :: xp0(*),yp0(*),zp0(*),fpx(*),fpy(*),fpz(*),U_extra
    INTEGER :: i,la,lb
    REAL(8)  :: xr1,xr2,yr1,yr2,zr1,zr2,x21,y21,z21,rs21,uux1,uux2&
         &,uuy1,uuy2,uuz1,uuz2,ubond(2)
    REAL(8) ::  qforce,r0,K
    INTEGER :: n_1,n_2,n_3
    INTEGER, SAVE :: q=0
    

    n_1=SIZE(pop1)
    n_2=SIZE(pop2)
    n_3=SIZE(pop3)
    U_extra=0.0D0
    DO i=1,n_1
       K=pop1(i) % K
       r0=pop1(i) % r0
       la=pop1(i) % ind_a
       lb=pop1(i) % ind_b
       xr1=xp0(la)
       yr1=yp0(la)
       zr1=zp0(la)
       xr2=xp0(lb)
       yr2=yp0(lb)
       zr2=zp0(lb)
       
       x21=xr2-xr1
       y21=yr2-yr1
       z21=zr2-zr1
       rs21=DSQRT(x21**2+y21**2+z21**2)
       
       qforce=-2.0D0*K*(r0-rs21)
       uux1=x21/rs21
       uuy1=y21/rs21
       uuz1=z21/rs21
       uux2=-uux1
       uuy2=-uuy1
       uuz2=-uuz1
       fpx(la)=fpx(la)+qforce*uux1
       fpy(la)=fpy(la)+qforce*uuy1
       fpz(la)=fpz(la)+qforce*uuz1
       
       fpx(lb)=fpx(lb)+qforce*uux2
       fpy(lb)=fpy(lb)+qforce*uuy2
       fpz(lb)=fpz(lb)+qforce*uuz2
       U_extra=U_extra+K*(rs21-r0)**2
    END DO
    IF(MOD(q,100) == 0) WRITE(*,*) rs21
    DO i=1,n_2
       K=pop2(i) % K
       r0=pop2(i) % r0
       la=pop2(i) % ind_a
       lb=pop2(i) % ind_b
       xr1=xp0(la)
       yr1=yp0(la)
       zr1=zp0(la)
       xr2=xp0(lb)
       yr2=yp0(lb)
       zr2=zp0(lb)
       
       x21=xr2-xr1
       y21=yr2-yr1
       z21=zr2-zr1
       rs21=DSQRT(x21**2+y21**2+z21**2)
       
       qforce=-2.0D0*K*(r0-rs21)
       uux1=x21/rs21
       uuy1=y21/rs21
       uuz1=z21/rs21
       uux2=-uux1
       uuy2=-uuy1
       uuz2=-uuz1
       fpx(la)=fpx(la)+qforce*uux1
       fpy(la)=fpy(la)+qforce*uuy1
       fpz(la)=fpz(la)+qforce*uuz1
       
       fpx(lb)=fpx(lb)+qforce*uux2
       fpy(lb)=fpy(lb)+qforce*uuy2
       fpz(lb)=fpz(lb)+qforce*uuz2
       U_extra=U_extra+K*(rs21-r0)**2
    END DO
    IF(MOD(q,100) == 0) WRITE(*,*) rs21
    DO i=1,n_3
       K=pop3(i) % K
       r0=pop3(i) % r0
       la=pop3(i) % ind_a
       lb=pop3(i) % ind_b
       xr1=xp0(la)
       yr1=yp0(la)
       zr1=zp0(la)
       xr2=xp0(lb)
       yr2=yp0(lb)
       zr2=zp0(lb)
       
       x21=xr2-xr1
       y21=yr2-yr1
       z21=zr2-zr1
       rs21=DSQRT(x21**2+y21**2+z21**2)
       
       qforce=-2.0D0*K*(r0-rs21)
       uux1=x21/rs21
       uuy1=y21/rs21
       uuz1=z21/rs21
       uux2=-uux1
       uuy2=-uuy1
       uuz2=-uuz1
       fpx(la)=fpx(la)+qforce*uux1
       fpy(la)=fpy(la)+qforce*uuy1
       fpz(la)=fpz(la)+qforce*uuz1
       
       fpx(lb)=fpx(lb)+qforce*uux2
       fpy(lb)=fpy(lb)+qforce*uuy2
       fpz(lb)=fpz(lb)+qforce*uuz2
       U_extra=U_extra+K*(rs21-r0)**2
    END DO
    IF(MOD(q,100) == 0) WRITE(*,*) rs21
    q=q+1
  END SUBROUTINE force
END MODULE Module_Extra_Forces
