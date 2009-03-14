      SUBROUTINE P_buffer_detach

      IMPLICIT none

      INTEGER MPIBufferSize,MPIBufferSize_p,ierr
      PARAMETER (MPIBufferSize_p=1024)
      REAL*8 MPIBuffer(MPIBufferSize_p)
      COMMON /mpibuffer/ MPIBufferSize,MPIBuffer


*----------------------- EXECUTABLE STATEMENTS ------------------------*

*---Set MPIBuffer to 8Mb

      CALL MPI_BUFFER_DETACH(MPIBuffer,MPIBufferSize,ierr)

      RETURN
      END
