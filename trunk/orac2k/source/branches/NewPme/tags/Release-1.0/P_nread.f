      INTEGER FUNCTION P_nread(buf,isize,ipartner,itype,null)

************************************************************************
*   Time-stamp: <04/11/11 14:21:15 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Jul  9 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8 buf(*)
      INTEGER isize,ipartner,itype,null

*----------------------- VARIABLES IN COMMON --------------------------*

#ifdef MPI      
      include 'mpif.h'
#if defined _CRAY_ | defined T3E
      INTEGER*8 status(mpi_status_size),size,comm,isizea,ipartnera
     &     ,itypea,ierrora
#else
      integer status(mpi_status_size),size,comm,isizea,ipartnera
     &     ,itypea,ierrora
#endif

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER ierror

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      ipartnera=ipartner
      itypea=itype
      isizea=isize

      IF(ipartnera .EQ. -1) ipartnera = MPI_ANY_SOURCE
      IF(itypea .EQ. -1) itypea = MPI_ANY_TAG
      size=MPI_BYTE
      comm=MPI_COMM_WORLD
      CALL MPI_RECV(buf,isizea,size,ipartnera,itypea,comm,STATUS,ierrora
     &     )
      IF(ierrora .NE. 0) WRITE(6,'(''MPI_RECV error:'',i9)') ierrora
      P_nread = ierrora
#endif

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
