Module PmeRfft3d
#ifdef PARALLEL
  USE mpi
#endif
  Use Errors, Only: Add_Errors=>Add, Print_Errors, error_args, errmsg_f
  Real(8), Dimension (:,:,:), Allocatable, Save, Private :: zzr
  Real(8), Dimension (:), Allocatable, Private :: work
  Double Complex, Dimension (:,:,:), Allocatable, Private :: zzc
  
  Real(8), Private :: fn1,fn2,fn3
  Integer, Save, Private :: n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nsize,ncache=-1 &
       ,n_n(3),n3_local,n3_start,n3_end,n2_start,n2_local,total_work
  Integer(8), Save:: plan_forward,plan_backward
!!$  Include 'fftw_f77.i'
!!$!     This file contains Parameter statements for various constants
!!$!     that can be passed to FFTW routines.  You should include
!!$!     this file in any FORTRAN program that calls the fftw_f77
!!$!     routines (either directly or with an #include statement
!!$!     if you use the C preprocessor).

      integer FFTW_FORWARD,FFTW_BACKWARD
      parameter (FFTW_FORWARD=-1,FFTW_BACKWARD=1)

      integer FFTW_REAL_TO_COMPLEX,FFTW_COMPLEX_TO_REAL
      parameter (FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1)

      integer FFTW_ESTIMATE,FFTW_MEASURE
      parameter (FFTW_ESTIMATE=0,FFTW_MEASURE=1)

      integer FFTW_OUT_OF_PLACE,FFTW_IN_PLACE,FFTW_USE_WISDOM
      parameter (FFTW_OUT_OF_PLACE=0)
      parameter (FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16)

      integer FFTW_THREADSAFE
      parameter (FFTW_THREADSAFE=128)

!     Constants for the MPI wrappers:
      integer FFTW_TRANSPOSED_ORDER, FFTW_NORMAL_ORDER
      integer FFTW_SCRAMBLED_INPUT, FFTW_SCRAMBLED_OUTPUT
      parameter(FFTW_TRANSPOSED_ORDER=1, FFTW_NORMAL_ORDER=0)
      parameter(FFTW_SCRAMBLED_INPUT=8192)
      parameter(FFTW_SCRAMBLED_OUTPUT=16384)
Contains
  Function do_rfft3d(isign,z,na1,na2,na3,na3_start,na3_local&
       &,na2_start,na2_local,nda1,nda2,nda3,My_Comm) Result(out)

!!$***********************************************************************
!!$                                                                      *
!!$     Interface for Fftw version 2.1.5                                 *
!!$                                                                      *
!!$     Use Fftw routines to compute real-to-complex and complex-to-     *
!!$     real Dfts. The transform is not normalized!                      *
!!$     to obtain a normalized transform the output must be divided      *
!!$     by n1*n2*n3. Otherwise a call of gpfft with isign = 1 followed   *
!!$     by a call with isign = -1 will multiply the sequence by n1*n2*n3 *
!!$                                                                      *
!!$----------------------------------------------------------------------*
!!$                                                                      *
!!$     ARGUMENTS:                                                       *
!!$                                                                      *
!!$   isign  : +1  forward transform of a complex periodic sequence      *
!!$            -1  backward transform of a complex periodic sequence     *
!!$   z      : A complex array of length z(nd1,nd2,nd3) which contains   *
!!$            the sequence to be transformed                            *
!!$   na1    : The length of the transform along the x direction         *
!!$   na2    : The length of the transform along the y direction         *
!!$   na3    : The length of the transform along the z direction         *
!!$   nad1   : On return the x physical dimension of the complex array   *
!!$   nad2   : On return the y physical dimension of the complex array   *
!!$   nad3   : On return the z physical dimension of the complex array   *
!!$   iproca : Processor identifier                                      *
!!$   nproca : Number of processors doing the transform                  *
!!$                                                                      *
!!$---- Description of the Distributed Data -----------------------------*
!!$                                                                      *
!!$   Consider a 3D matrix A of size n1-by-n2-by-n3.  n1, n2, and n3     *
!!$   are the sizes of the matrix A along the X, Y, and Z dimensions,    *
!!$   respectively. The nprocs processors are only partitioned along     *
!!$   the direction Z. The input matrix A is distributed along the Z     *
!!$   dimensions and each processors handles n1*n2*npz elements of A     *
!!$   (npz=n3/nprocs).                                                   *
!!$                                                                      *
!!$***********************************************************************
    
    Implicit none
    Integer, Intent(in), Optional  :: isign,My_Comm
    Integer, Intent(in), Optional :: na1,na2,na3
    Real(8), Dimension (:,:,:), Intent(inout) :: z
    Integer, Intent(out), Optional :: nda1,nda2,nda3,na3_local&
         &,na3_start,na2_start,na2_local
    Logical :: out

    Real(8), Save :: times=0.0,elapse,cpu1,cpu2,dummy
    
    Integer, Save :: calls=0
    Integer :: i,j,k,ierr
    Integer :: My_Node
    
    out=.True.
    If(isign .Eq. 0) Then
       n1=na1
       n2=na2
       n3=na3
       n_n(1)=n1
       n_n(2)=n2
       n_n(3)=n3
       Call set_fftw
       na3_start=n3_start
       na3_local=n3_local
       na2_start=n2_start
       na2_local=n2_local
       Return
    End If
    call Mpi_comm_rank(Mpi_comm_world, My_node, ierr )
    If(ncache == -1) Then
       errmsg_f='do_rfft3d not initialized. Abort now!'
       Call Add_Errors(-1,errmsg_f)
       out=.False.
    End If

    calls=calls+1
#if defined PARALLEL
    If(isign == 1) Then
       Call rfftwnd_f77_mpi(plan_forward,1,z,work,0,Fftw_transposed_order)
    Else If(isign == -1) Then
       Call rfftwnd_f77_mpi(plan_backward,1,z,work,0,Fftw_transposed_order)
    End If
#else
    If(isign == 1) Then
       Allocate(zzr(n1,n2,n3))
       Allocate(zzc(n1/2+1,n2,n3))
       zzr(1:n1,1:n2,1:n3)=z(1:n1,1:n2,1:n3)
       Call rfftwnd_f77_one_real_to_complex(plan_forward,zzr,zzc)
       Do k=1,n3
          Do j=1,n2
             Do i=1,n1/2+1
                z((i-1)*2+1,j,k)=Dreal(zzc(i,j,k))
                z((i-1)*2+2,j,k)=Dimag(zzc(i,j,k))
             End Do
          End Do
       End Do
       Deallocate(zzr,zzc)
    Else If(isign == -1) Then
       Allocate(zzr(n1,n2,n3))
       Allocate(zzc(n1/2+1,n2,n3))
       Do k=1,n3
          Do j=1,n2
             Do i=1,n1/2+1
                zzc(i,j,k)=Dcmplx(z((i-1)*2+1,j,k),z((i-1)*2+2,j,k))
             End Do
          End Do
       End Do
       Call rfftwnd_f77_one_complex_to_real(plan_backward,zzc,zzr)
       z(1:n1,1:n2,1:n3)=zzr(1:n1,1:n2,1:n3)
       Deallocate(zzr,zzc)
    End If
#endif

  Contains
    Subroutine set_fftw
#if defined PARALLEL
      nda1=(n1/2+1)*2
      nda2=n2
      ncache=1
      Call rfftw3d_f77_mpi_create_plan(plan_backward,My_Comm&
           &,n1,n2,n3,Fftw_complex_to_real,Fftw_measure)
      Call rfftw3d_f77_mpi_create_plan(plan_forward ,My_Comm&
           &,n1,n2,n3,Fftw_real_to_complex,Fftw_measure)
      Call rfftwnd_f77_mpi_local_sizes(plan_forward,n3_local,n3_start&
           &,n2_local,n2_start,total_work)
      Allocate(work(total_work))
      n3_start=n3_start+1
      n2_start=n2_start+1
      n3_end=n3_start+n3_local-1
      nda3=n3_local
      nd1=nda1
      nd2=nda2
      nd3=n3_local
#else
      nda1=(n1/2+1)*2
      nda2=n2
      nda3=n3
      ncache=1
      Call rfftwnd_f77_create_plan(plan_forward,3,n_n,Fftw_real_to_complex &
           ,Fftw_estimate)
      Call rfftwnd_f77_create_plan(plan_backward,3,n_n,Fftw_complex_to_real &
           ,Fftw_estimate)
      n3_start=1
      n3_end=n3
      n3_local=n3
      n2_start=1
      n2_local=n2
#endif
      Return
    End Subroutine set_fftw
  End Function do_rfft3d
  Subroutine apply_symmetry(z)
    Implicit none
    Double Complex, Dimension (:,:,:) :: z

    Integer :: n,i,j,k
    Integer :: nf1,ic,jc,kc
    Logical :: k_l,j_l,i_l

    nf1 = n1/2+1

    Do k=1,n3
       kc=Mod(-k+n3+1,n3)+1
       Do j=1,n2
          jc=Mod(-j+n2+1,n2)+1
          Do i=1,nf1
             ic=Mod(-i+n1+1,n1)+1
             z(ic,jc,kc)=Conjg(z(i,j,k))
          End Do
       End Do
    End Do
  End Subroutine apply_symmetry
End Module PmeRfft3d
