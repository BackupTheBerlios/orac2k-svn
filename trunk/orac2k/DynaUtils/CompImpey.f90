PROGRAM CompImpey

!!$***********************************************************************
!!$   Time-stamp: <2008-03-20 16:13:01 marchi>                           *
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

  REAL(4) :: x,y,z
  REAL(4) :: dip_x,dip_y,dip_z
  TYPE Traj
     INTEGER :: n0
     INTEGER, ALLOCATABLE :: Ind(:)
  END type Traj
  TYPE(Traj), ALLOCATABLE, SAVE :: Tr(:)

  REAL(8) :: dipn_x,dipn_y,dipn_z,dipo_x,dipo_y,dipo_z,corra
  INTEGER :: pp,oo,ia,counter,n0,counter_ia
  REAL(8), DIMENSION(:), ALLOCATABLE, SAVE :: corr_1
  REAL(8) :: fstep,fstep_o,dfstep,norma,normb,x0,y0,z0,x1,y1,z1,dist,ratio
  LOGICAL :: ok,stdout
  CHARACTER(80), SAVE :: file_out = 'IMPEY.dat'

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

  
  CALL Open_Input
  INQUIRE(unit=kread,NAMED=stdout)

  OPEN(unit=kout,FILE=file_out,form='FORMATTED')

  READ(kbin) m,n

  ALLOCATE(Tr(n))
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

     READ(kbin) fstep,n0
     Tr(nn) % n0=n0
     ALLOCATE(tr(nn) % Ind(n0))
     READ(kbin) (tr(nn) % Ind(i),i=1,n0)     
     READ(kbin) (x, y, z, dip_x, dip_y, dip_z,i=1,n0)
     nn=nn+1
  END DO
  dfstep=fstep-fstep_o

!!$
!!$ --- Make Mask
!!$

  mask=.FALSE.
  DO nn=1,n
     DO ii=1,Tr(nn) % n0
        i=tr(nn) % Ind(ii)
        mask(nn,i)=.TRUE.
     END DO
  END DO

  DO i=1,m
     ia=index_sv(1+i)
     IF(MOD(i,1000) == 0) WRITE(*,'(''----> Molecule  '',i4,'' to follow'')') ia
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
  counter=0
  DO i=1,m
     ia=index_sv(1+i)
     counter_ia=0
     DO nn=1,n
        IF(Mask(nn,i)) THEN
           counter_ia=counter_ia+1
        END IF
     END DO
     ratio=DBLE(counter_ia)/DBLE(n)
     IF( ratio > 0.5D0 .AND. ratio < 0.6D0) THEN
        counter=counter+1
        WRITE(*,'(''Att. more than 80 %, Molecules No.'',i6,'' Seq. No. '',i6)') counter,ia
     END IF
  END DO


  WRITE(*,'(''----> N(t) '')') 
  WRITE(*,*) n
  DO i=0,n-1
     IF(stdout) THEN
        IF(corr_1(i) /= 0.0D0) WRITE(*,'(f12.3,3f25.18)') DBLE(i)*dfstep/1000.0D0&
             &,corr_1(i)
     ELSE
        IF(corr_1(i) /= 0.0D0) WRITE(kout,'(f12.3,3f25.18)') DBLE(i)&
             &*dfstep/1000.0D0,corr_1(i)
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
