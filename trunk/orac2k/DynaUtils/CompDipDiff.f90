PROGRAM CompDipDiff

!!$***********************************************************************
!!$   Time-stamp: <2008-03-21 11:09:19 marchi>                           *
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
  
  INTEGER, SAVE :: kread=-1,kbin=-1,Corr_Max=-1,n_Max=-1,kunit=-1
  CHARACTER(80) :: filename='DipDiff.dat'
  INTEGER :: m,n,ii,nn,mm,i
  INTEGER, ALLOCATABLE, DIMENSION(:) :: index1,index_sv
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: Mask

  TYPE Dyna
     REAL(4) :: x,y,z
     REAL(4) :: dip_x,dip_y,dip_z
  END TYPE Dyna
  TYPE Traj
     INTEGER :: n0
     INTEGER, DIMENSION(:), ALLOCATABLE :: Ind
     TYPE(Dyna), DIMENSION(:), ALLOCATABLE :: Coord
  END type Traj

  TYPE(Traj), DIMENSION(:), ALLOCATABLE, SAVE :: tr

  REAL(8) :: dipn_x,dipn_y,dipn_z,dipo_x,dipo_y,dipo_z,corra
  INTEGER :: pp,oo,ia,nn_end
  REAL(8), DIMENSION(:), ALLOCATABLE :: corr_1,corr_2,diff,data2,dip1,dip2
  INTEGER, DIMENSION(:), ALLOCATABLE :: ncorr,ncorr2
  REAL(8) :: fstep,fstep_o,dfstep,norma,x0,y0,z0,x1,y1,z1,dist,aux1
  INTEGER :: n0,n1,n2,j,jj
  LOGICAL :: ok_a,ok_b
  REAL(4), ALLOCATABLE :: xyz(:,:),dip(:,:)
  INTEGER :: kmax,k

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

  
  CALL Open_Input

  IF(kread /= -1) THEN
     READ(kread,*) Corr_Max,N_Max
     WRITE(*,*) 'Max Correlation: ',Corr_Max,'Max Configuration No.: ',n_Max
  END IF

  READ(kbin) m,n

  IF(N_Max > 0) n=n_Max

  ALLOCATE(xyz(3,n),dip(3,n))
  ALLOCATE(Mask(n,m))

  ALLOCATE(Tr(n))
  ALLOCATE(index_sv(m+1))
  ALLOCATE(corr_1(0:n-1),corr_2(0:n-1),diff(0:n-1),ncorr(0:n-1))
  ALLOCATE(data2(0:n-1),ncorr2(0:n-1),dip1(0:n-1),dip2(0:n-1))
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
     READ(kbin) fstep,n0
     Tr(nn) % n0=n0
     ALLOCATE(tr(nn) % Coord (n0))
     ALLOCATE(tr(nn) % Ind(n0))

     READ(kbin) (tr(nn) % Ind(i),i=1,n0)
     READ(kbin) (tr(nn) % Coord(i) % x, tr(nn) % Coord(i) % y, tr(nn)&
          & % Coord(i) % z,tr(nn) % Coord(i) % dip_x, tr(nn) %&
          & Coord(i) % dip_y, tr(nn) % Coord(i) % dip_z,i=1,n0)

     DO i=1,n0
        dipn_x=tr(nn) % Coord(i) % dip_x
        dipn_y=tr(nn) % Coord(i) % dip_y
        dipn_z=tr(nn) % Coord(i) % dip_z
        norma=1.0D0/SQRT(dipn_x**2+dipn_y**2+dipn_z**2)
        tr(nn) % Coord(i) % dip_x=dipn_x*norma
        tr(nn) % Coord(i) % dip_y=dipn_y*norma
        tr(nn) % Coord(i) % dip_z=dipn_z*norma
     END DO
     IF(MOD(nn,1000) == 0) WRITE(*,*) nn
     nn=nn+1
  END DO

  mask=.FALSE.
  DO nn=1,n
     DO i=1,Tr(nn) % n0
        ii=tr(nn) % Ind(i)
        mask(nn,ii)=.TRUE.
     END DO
  END DO

  dfstep=fstep-fstep_o
  ncorr=0
  diff=0.0D0
  DO ii=1,m
     ncorr2=0
     data2=0.0D0
     dip1=0.0D0
     dip2=0.0D0
     ia=index_sv(1+ii)
     DO nn=1,n
        xyz(:,nn)=0.0D0
        dip(:,nn)=0.0D0
        IF(.NOT. mask(nn,ii)) CYCLE
        ok_b=.FALSE.
        n1=tr(nn) % n0
        DO j=1,n1
           jj=Tr(nn) % Ind (j)
           IF(jj == ii) THEN
              xyz(1,nn)=tr(nn) % Coord(j) % x
              xyz(2,nn)=tr(nn) % Coord(j) % y
              xyz(3,nn)=tr(nn) % Coord(j) % z
              dip(1,nn)=tr(nn) % Coord(j) % dip_x
              dip(2,nn)=tr(nn) % Coord(j) % dip_y
              dip(3,nn)=tr(nn) % Coord(j) % dip_z
              ok_b=.TRUE.
              EXIT
           END IF
        END DO
        IF(.NOT. ok_b) THEN
           WRITE(*,*) 'It shouldn''t have stopped here!!'
           STOP
        END IF
     END DO
     IF(MOD(ii,500) == 0) WRITE(*,*) ii,ia
     DO nn=1,n
        IF(mask(nn,ii)) THEN
           dipn_x=dip(1,nn)
           dipn_y=dip(2,nn)
           dipn_z=dip(3,nn)
           x0=xyz(1,nn)
           y0=xyz(2,nn)
           z0=xyz(3,nn)
           nn_end=n
           IF(Corr_Max > 0) THEN
              nn_end=nn+Corr_Max
           END IF
           IF(nn_end > n) nn_end=n
           DO oo=nn,nn_end
              pp=oo-nn
              IF(.NOT. mask(oo,ii)) EXIT
              dipo_x=dip(1,oo)
              dipo_y=dip(2,oo)
              dipo_z=dip(3,oo)
              
              corra=dipn_x*dipo_x+dipn_y*dipo_y+dipn_z*dipo_z
              dip1(pp)=dip1(pp)+corra
              dip2(pp)=dip2(pp)+1.5D0*corra**2-0.5D0              

              x1=xyz(1,oo)
              y1=xyz(2,oo)
              z1=xyz(3,oo)
              dist=(x0-x1)**2+(y0-y1)**2+(z0-z1)**2
              data2(pp)=data2(pp)+dist
              ncorr2(pp)=ncorr2(pp)+1
           END DO
        END IF
     END DO
     Kmax=-1
     DO k=1,n-1
        IF(data2(k) .EQ. 0.0D0) THEN
           Kmax=k-1
           GOTO 20000
        END IF
     END DO
20000 CONTINUE
     IF(Kmax .EQ. -1) Kmax=n-1
     IF(Kmax .NE. 0) THEN
        DO k=0,Kmax
           aux1=DFLOAT(ncorr2(k))
           data2(k)=data2(k)/aux1
           diff(k)=diff(k)+data2(k)

           dip1(k)=dip1(k)/aux1
           dip2(k)=dip2(k)/aux1
           corr_1(k)=corr_1(k)+dip1(k)
           corr_2(k)=corr_2(k)+dip2(k)

           ncorr(k)=ncorr(k)+1
        END DO
        DO k=Kmax+1,n-1
           diff(k)=diff(k)+data2(Kmax)
           corr_1(k)=corr_1(k)+dip1(Kmax)
           corr_2(k)=corr_2(k)+dip2(Kmax)
           ncorr(k)=ncorr(k)+1
        END DO
     END IF
  END DO
  WRITE(kunit,'(''#      Time         Dip Rank 1                Dip Ra&
       &nk 2             <|r(t)-r(0)|^2>'')')
  IF(Corr_Max < 0) Corr_Max=n
  DO i=0,Corr_Max-1
     norma=1.0D0/DBLE(ncorr(i))
     IF(ncorr(i) /= 0) WRITE(kunit,'(f12.3,3f25.18)') DBLE(i)*dfstep&
          &/1000.0D0,corr_1(i)*norma,corr_2(i)*norma,diff(i)*norma
  END DO

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

CONTAINS
  SUBROUTINE Open_Input
    IMPLICIT none
    
!!$----------------------------- ARGUMENTS ------------------------------*

!!$------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER :: IARGC,nargs,nchars,n,n_in,n_out,nargs_e
    CHARACTER(80) :: buffer,fileinput,command,labels
    CHARACTER(120) :: errmsg
    LOGICAL, SAVE :: ok=.FALSE.

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*


    nargs=IARGC()
    n=0
    DO WHILE(n < nargs) 
       n=n+1
       CALL GETARG(n, buffer)
       IF(TRIM(buffer) == '>') THEN
          n=n-1
          EXIT
       END IF
    END DO
    nargs=n
    n=0
    n_in=0
    n_out=0
    DO WHILE(n < nargs) 
       n=n+1
       CALL GETARG(n, buffer)
       labels=TRIM(buffer)
       IF(labels(1:2) == '-i') THEN
          n_in=n+1
       END IF
       IF(labels(1:2) == '-o') THEN
          n_out=n+1
       END IF
    END DO
    IF(n_out == 0) THEN
       kunit=12
       OPEN(unit=kunit,file=filename,form='FORMATTED')
    END IF

    SELECT CASE(nargs)
    CASE(1)
       CALL GETARG(1, buffer)
       fileinput=TRIM(buffer)
       kbin=10
       OPEN(unit=kbin,file=fileinput,form='UNFORMATTED')
    CASE(3)
       CALL GETARG(3, buffer)
       fileinput=TRIM(buffer)
       kbin=10
       OPEN(unit=kbin,file=fileinput,form='UNFORMATTED')
       CALL GETARG(1, buffer)
       labels=TRIM(buffer)
       n=0
       DO WHILE(n < 2) 
          n=n+1
          CALL GETARG(n, buffer)
          labels=TRIM(buffer)
          IF(labels(1:2) == '-i') THEN
             n=n+1
             CALL GETARG(n, buffer)
             kread=5
             fileinput=TRIM(buffer)
             OPEN(unit=kread,file=filename,form='FORMATTED')
          ELSE IF(labels(1:2) == '-o') THEN
             n=n+1
             CALL GETARG(n, buffer)
             kunit=12
             fileinput=TRIM(buffer)
             OPEN(unit=kunit,file=fileinput,form='FORMATTED')
          END IF
       END DO
    CASE(5)
       CALL GETARG(5, buffer)
       fileinput=TRIM(buffer)
       kbin=10
       OPEN(unit=kbin,file=fileinput,form='UNFORMATTED')
       
       CALL GETARG(1, buffer)
       labels=TRIM(buffer)
       n=0
       DO WHILE(n < 4) 
          n=n+1
          CALL GETARG(n, buffer)
          labels=TRIM(buffer)
          IF(labels(1:2) == '-i') THEN
             n=n+1
             CALL GETARG(n, buffer)
             kread=5
             fileinput=TRIM(buffer)
             OPEN(unit=kread,file=fileinput,form='FORMATTED')
          ELSE IF(labels(1:2) == '-o') THEN
             n=n+1
             CALL GETARG(n, buffer)
             kunit=12
             fileinput=TRIM(buffer)
             OPEN(unit=kunit,file=fileinput,form='FORMATTED')
          END IF
       END DO
    CASE DEFAULT
       CALL GETARG(0, buffer)
       command = TRIM(buffer)
       nchars=LEN_TRIM(buffer)
       errmsg='Usage: '// command(1:nchars) //' [-i input&
            &1 -o DipDiff] dynamic_file [ > output ]'
       WRITE(*,'(a)') TRIM(errmsg)
       STOP
    END SELECT

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*


  END SUBROUTINE Open_Input
END PROGRAM CompDipDiff
