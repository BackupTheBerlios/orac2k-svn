SUBROUTINE do_pfft3d(isign,z,na1,na2,na3,nda1,nda2,nda3,iproca,nproca)
!***********************************************************************
!                                                                      *
!     Interface for Goedecker's Parallel FFT's                         *
!                                                                      *
!     Subroutine gpfft computes forward and backward complex           *
!     discrete Fourier transform. The transform is not normalized.     *
!     to obtain a normalized transform the output must be divided      *
!     by n1*n2*n3. Otherwise a call of gpfft with isign = 1 followed   *
!     by a call with isign = -1 will multiply the sequence by n1*n2*n3 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     ARGUMENTS:                                                       *
!                                                                      *
!   isign  : +1  forward transform of a complex periodic sequence      *
!            -1  backward transform of a complex periodic sequence     *
!   n1     : The length of the transform along the x direction         *
!   n2     : The length of the transform along the y direction         *
!   n3     : The length of the transform along the z direction         *
!   nd1    : The first physical dimensions of the complex array, z,    *
!            to be transformed                                         *
!   nproc  : Number of processors doing the transform                  *
!   z      : A complex array of length z(nd1,nd2,nd3) which contains   *
!            the sequence to be transformed                            *
!                                                                      *
!---- Description of the Distributed Data -----------------------------*
!                                                                      *
!   Consider a 3D matrix A of size n1-by-n2-by-n3.  n1, n2, and n3     *
!   are the sizes of the matrix A along the X, Y, and Z dimensions,    *
!   respectively. The nprocs processors are only partitioned along     *
!   the direction Z. The input matrix A is distributed along the Z     *
!   dimensions and each processors handles n1*n2*npz elements of A     *
!   (npz=n3/nprocs).                                                   *
!                                                                      *
!---- Tests performed so far ------------------------------------------*
!                                                                      *
!   The routine has been tested on a T3E machine. It also works        *
!   on serial machines on single processor. It needs be compiled       *
!   with a precompiler instruction -DPARALLEL and -DNOPARALLEL         *
!   on parallel and serial systems, respectively.                      *
!                                                                      *
!***********************************************************************

  IMPLICIT none
  INTEGER, INTENT(in) :: isign
  INTEGER, INTENT(in), OPTIONAL :: na1,na2,na3,iproca,nproca
  REAL(8), DIMENSION (*), INTENT(inout) :: z
  INTEGER, INTENT(out), OPTIONAL :: nda1,nda2,nda3
  
  REAL(8), DIMENSION (:), ALLOCATABLE :: zt
  REAL(8), DIMENSION (:), ALLOCATABLE, SAVE :: zw

  REAL :: fn1,fn2,fn3
  INTEGER, SAVE :: n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nsize,ncache
  INTEGER, SAVE :: iproc,nproc
  REAL(8) :: time1_avg,time2_avg,time_avg
  INTEGER, SAVE :: ntime_avg=0

  IF(isign .EQ. 0) THEN
     n1=na1
     n2=na2
     n3=na3
     iproc=iproca
     nproc=nproca
     CALL set_gpfft
     RETURN
  END IF
  

  ALLOCATE(zt(nsize))
  CALL gpfft3d(n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nsize,iproc,nproc &
       ,z,zt,isign,zw,ncache)
  DEALLOCATE(zt)
  RETURN
  CONTAINS
    SUBROUTINE set_gpfft
      ncache=1024*4
      nd1=n1
      nd2=n2
      nd3=n3
      md1=n1
      md2=n2
      md3=n3
      nsize=2*MAX(nd1,md1)*MAX(nd2,md2)*MAX(nd3,md3)/nproc
      nda1=MAX(nd1,md1)
      nda2=MAX(nd2,md2)
      nda3=MAX(nd3,md3)/nproc
      ALLOCATE(zw(ncache))
      IF(iproc .EQ. 0) WRITE(*,*) 'Using Goedecker''s fft Code'
      RETURN
    END SUBROUTINE set_gpfft
END Subroutine do_pfft3d
