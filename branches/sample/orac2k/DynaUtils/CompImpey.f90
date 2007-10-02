PROGRAM CompImpey

!!$***********************************************************************
!!$   Time-stamp: <2007-09-19 15:58:47 marchi>                           *
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
  
  INTEGER, SAVE :: kread=5,kbin=10,kprint=6,kout=11
  CHARACTER(80) :: filename
  INTEGER :: m,n,ii,nn,mm,i,nn_old,nn_new
  INTEGER, ALLOCATABLE, DIMENSION(:) :: index1,index_sv
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: Mask

  TYPE Dyna
     REAL(8) :: x,y,z
     REAL(8) :: dip_x,dip_y,dip_z
  END TYPE Dyna

  TYPE(Dyna), DIMENSION(:,:), ALLOCATABLE, SAVE :: Coord
  REAL(8) :: dipn_x,dipn_y,dipn_z,dipo_x,dipo_y,dipo_z,corra
  INTEGER :: pp,oo,ia,counter
  REAL(8), DIMENSION(:), ALLOCATABLE, SAVE :: corr_1
  REAL(8) :: fstep,fstep_o,dfstep,norma,normb,x0,y0,z0,x1,y1,z1,dist
  LOGICAL :: ok,stdout
  CHARACTER(80), SAVE :: file_out = 'IMPEY.dat'

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

  
  CALL Open_Input
  INQUIRE(unit=kread,NAMED=stdout)

  OPEN(unit=kout,FILE=file_out,form='FORMATTED')

  READ(kbin) m,n

  ALLOCATE(Coord(n,m))
  ALLOCATE(Mask(n,m))
  ALLOCATE(index_sv(m+1))
  ALLOCATE(corr_1(0:n-1))
  corr_1=0.0D0

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
     READ(kbin) (mask(nn,mm),mm=1,m)
     IF(MOD(nn,1000) == 0) WRITE(*,*) nn
     nn=nn+1
  END DO
  dfstep=fstep-fstep_o

  DO i=1,m
     ia=index_sv(1+i)
     WRITE(*,'(''----> Molecule  '',i4,'' to follow'')') ia
     nn=1
     nn_old=1
     DO WHILE(nn < n)
        DO WHILE(Mask(nn,i) .AND. nn < n)
           nn=nn+1
        END DO
        nn_new=nn-1
        IF(nn == n) nn_new=n

        DO oo=nn_old,nn_new
           corr_1(nn_new-oo)=corr_1(nn_new-oo)+DBLE(oo-nn_old+1)/DBLE(n-nn_new+oo)
        END DO
        nn=nn+1
        nn_old=nn
     END DO
  END DO
  counter=0
  DO i=1,m
     ia=index_sv(1+i)
     ok=.TRUE.
     DO nn=1,n
        IF(.NOT. Mask(nn,i)) ok=.FALSE.
     END DO
     IF(ok) THEN
        counter=counter+1
        WRITE(*,'(''Attached Molecules No.'',i6,'' Seq. No. '',i6)') counter,ia
     END IF
  END DO
  WRITE(*,'(''----> N(t) '')') 
  DO i=0,n-1
     IF(stdout) THEN
        IF(corr_1(i) /= 0.0D0) WRITE(*,'(f12.3,3f25.18)') DBLE(i)*dfstep&
             &,corr_1(i)
     ELSE
        IF(corr_1(i) /= 0.0D0) WRITE(kout,'(f12.3,3f25.18)') DBLE(i)&
             &*dfstep,corr_1(i)
     END IF
  END DO

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

CONTAINS
  SUBROUTINE Open_Input
    IMPLICIT none
    
!!$----------------------------- ARGUMENTS ------------------------------*

    INTEGER :: IARGC,nargs,nchars
    CHARACTER(80) :: buffer,fileinput,errmsg,command
    LOGICAL, SAVE :: ok=.FALSE.

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*


    nargs=IARGC()
    IF(nargs > 0) THEN
       ok=.TRUE.
       CALL GETARG(1, buffer)
       fileinput=TRIM(buffer)
       OPEN(unit=kbin,file=fileinput,form='UNFORMATTED')
    END IF
    
    IF(.NOT. ok) THEN
       CALL GETARG(0, buffer)
       command = TRIM(buffer)
       nchars=LEN_TRIM(buffer)
       errmsg='Usage: '// command(1:nchars) //' DynamicFile [Input] [ > output ]'
       WRITE(*,'(a)') errmsg
       STOP
    END IF

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*


  END SUBROUTINE Open_Input
END PROGRAM CompImpey
