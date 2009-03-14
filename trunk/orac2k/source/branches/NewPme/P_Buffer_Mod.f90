MODULE P_Buffer_Mod
  REAL(8), DIMENSION(:), ALLOCATABLE, SAVE :: MPIBuffer
  INTEGER, SAVE :: MPIBuffer_Size=1024*1024*32
CONTAINS
  SUBROUTINE P_buffer_attach
    IMPLICIT none
    INTEGER :: natom
    
    INTEGER :: ierr

!!$*----------------------- EXECUTABLE STATEMENTS ------------------------*
!!$
!!$*---Set MPIBuffer to 8Mb

    
    ALLOCATE(MPIBuffer(MPIBuffer_Size))

    CALL MPI_BUFFER_ATTACH(MPIBuffer,MPIBuffer_Size,ierr)
    IF(ierr /= 0) THEN
       WRITE(*,*) 'Stops in P_Buffer_attach!'
       STOP
    END IF
  END SUBROUTINE P_buffer_attach
  SUBROUTINE P_buffer_detach

    IMPLICIT none

    INTEGER :: ierr

!!$*----------------------- EXECUTABLE STATEMENTS ------------------------*
!!$
!!$*---Set MPIBuffer to 8Mb

    CALL MPI_BUFFER_DETACH(MPIBuffer,MPIBuffer_Size,ierr)
    DEALLOCATE(MPIBuffer)
  END SUBROUTINE P_buffer_detach
END MODULE P_Buffer_Mod
