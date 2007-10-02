PROGRAM CompDipDiff

!!$***********************************************************************
!!$   Time-stamp: <2007-09-18 14:47:40 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Sep 14 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program Utilities ----*


!!$======================== DECLARATIONS ================================*

  IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

!!$----------------------- VARIABLES IN COMMON --------------------------*

!!$------------------------- LOCAL VARIABLES ----------------------------*
  
  INTEGER, SAVE :: kread=5,kbin=10
  CHARACTER(80) :: filename
  INTEGER :: m,n,ii,nn,mm,i
  INTEGER, ALLOCATABLE, DIMENSION(:) :: index1,index_sv
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: Mask

  TYPE Dyna
     REAL(8) :: x,y,z
     REAL(8) :: dip_x,dip_y,dip_z
  END TYPE Dyna

  TYPE(Dyna), DIMENSION(:,:), ALLOCATABLE, SAVE :: Coord
  REAL(8) :: dipn_x,dipn_y,dipn_z,dipo_x,dipo_y,dipo_z,corra
  INTEGER :: pp,oo,ia,nn_end
  REAL(8), DIMENSION(:), ALLOCATABLE :: corr_1,corr_2,diff
  INTEGER, DIMENSION(:), ALLOCATABLE :: ncorr
  REAL(8) :: fstep,fstep_o,dfstep,norma,x0,y0,z0,x1,y1,z1,dist

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

  
  CALL Open_Input(kbin)

  READ(kbin) m,n

!!$  n=30000
  ALLOCATE(Coord(n,m))
  ALLOCATE(Mask(n,m))
  ALLOCATE(index_sv(m+1))
  ALLOCATE(corr_1(0:n-1),corr_2(0:n-1),diff(0:n-1),ncorr(0:n-1))
  corr_1=0.0D0
  corr_2=0.0D0
  diff=0.0D0
  ncorr=0

  index_sv(1)=m
  READ(kbin) (index_sv(1+ii),ii=1,m)


  fstep=0.0D0
  nn=1
  DO WHILE(nn <= n)
     fstep_o=fstep
     READ(kbin) fstep
     READ(kbin) (Coord(nn,i) % x, Coord(nn,i) % y, Coord(nn,i) % z,&
          & Coord(nn,i) % dip_x, Coord(nn,i) % dip_y, Coord(nn,i) % dip_z,&
          & i=1,m)

     DO i=1,m
        dipn_x=Coord(nn,i) % dip_x
        dipn_y=Coord(nn,i) % dip_y
        dipn_z=Coord(nn,i) % dip_z
        norma=1.0D0/SQRT(dipn_x**2+dipn_y**2+dipn_z**2)
        Coord(nn,i) % dip_x=dipn_x*norma
        Coord(nn,i) % dip_y=dipn_y*norma
        Coord(nn,i) % dip_z=dipn_z*norma
     END DO

     READ(kbin) (mask(nn,mm),mm=1,m)
     IF(MOD(nn,1000) == 0) WRITE(*,*) nn
     nn=nn+1
  END DO
  dfstep=fstep-fstep_o
  DO i=1,m
     ia=index_sv(1+i)
     WRITE(*,'(''----> Atom  '',i4,'' to follow'')') ia
     DO nn=1,n
        dipn_x=Coord(nn,i) % dip_x
        dipn_y=Coord(nn,i) % dip_y
        dipn_z=Coord(nn,i) % dip_z
        x0=Coord(nn,i) % x
        y0=Coord(nn,i) % y
        z0=Coord(nn,i) % z
        IF(.NOT. Mask(nn,i)) THEN
           nn_end=nn+2000
           IF(nn_end > n) nn_end=n
           DO oo=nn,nn_end
              IF(Mask(oo,i)) EXIT
              dipo_x=Coord(oo,i) % dip_x
              dipo_y=Coord(oo,i) % dip_y
              dipo_z=Coord(oo,i) % dip_z
              corra=dipn_x*dipo_x+dipn_y*dipo_y+dipn_z*dipo_z
              pp=oo-nn
              corr_1(pp)=corr_1(pp)+corra
              corr_2(pp)=corr_2(pp)+1.5D0*corra**2-0.5D0

              x1=Coord(oo,i) % x
              y1=Coord(oo,i) % y
              z1=Coord(oo,i) % z

              dist=(x0-x1)**2+(y0-y1)**2+(z0-z1)**2
              diff(pp)=diff(pp)+dist

              ncorr(pp)=ncorr(pp)+1
           END DO
        END IF
     END DO
  END DO
  
  DO i=0,n-1
     norma=1.0D0/DBLE(ncorr(i))
     IF(ncorr(i) /= 0) WRITE(*,'(f12.3,3f25.18)') DBLE(i)*dfstep&
          &,corr_1(i)*norma,corr_2(i)*norma,diff(i)*norma
  END DO

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

CONTAINS
  SUBROUTINE Open_Input(kread)
    IMPLICIT none
    
!!$----------------------------- ARGUMENTS ------------------------------*

    INTEGER :: kread

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: IARGC,nargs,nchars
    CHARACTER(80) :: buffer,fileinput,errmsg,command
    LOGICAL, SAVE :: ok=.FALSE.

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*


    nargs=IARGC()
    IF(nargs >= 1) THEN
       ok=.TRUE.
       CALL GETARG(1, buffer)
       fileinput=TRIM(buffer)
       OPEN(unit=kread,file=fileinput,form='UNFORMATTED')
    END IF
    
    IF(.NOT. ok) THEN
       CALL GETARG(0, buffer)
       command = TRIM(buffer)
       nchars=LEN_TRIM(buffer)
       errmsg='Usage: '// command(1:nchars) //' input [ > output ]'
       WRITE(*,'(a)') errmsg
       STOP
    END IF

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*


  END SUBROUTINE Open_Input
END PROGRAM CompDipDiff
