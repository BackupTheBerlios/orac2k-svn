      INTEGER MPIBufferSize,MPIBufferSize_p
      PARAMETER (MPIBufferSize_p=1024)
      REAL*8 MPIBuffer(MPIBufferSize_p)
      COMMON /mpibuffer/ MPIBufferSize,MPIBuffer
