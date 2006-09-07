      SUBROUTINE P_buffer_attach

      IMPLICIT none

      INTEGER MPIBufferSize,MPIBufferSize_p,ierr
      PARAMETER (MPIBufferSize_p=1024*1024)
      REAL*8 MPIBuffer(MPIBufferSize_p)
      COMMON /mpibuffer/ MPIBuffer,MPIBufferSize

*----------------------- EXECUTABLE STATEMENTS ------------------------*

*---Set MPIBuffer to 8Mb
      MPIBufferSize=8*MPIBufferSize_p

      CALL MPI_BUFFER_ATTACH(MPIBuffer,MPIBufferSize,ierr)
      
      RETURN
      END
