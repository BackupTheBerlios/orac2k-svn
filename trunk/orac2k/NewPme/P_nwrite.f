      INTEGER FUNCTION P_nwrite(buf,isize,ipartner,itype,null)

************************************************************************
*   Time-stamp: <2009-03-14 18:32:26 marchi>                             *
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
      INTEGER isize,ipartner,null,itype


*----------------------- VARIABLES IN COMMON --------------------------*

#ifdef __MPI      
      include 'mpif.h'

*------------------------- LOCAL VARIABLES ----------------------------*

#if defined _CRAY_ | defined T3E
      INTEGER*8 status(mpi_status_size),size,comm,isizea,ipartnera
     &     ,itypea,ierrora
#else
      integer status(mpi_status_size),size,comm,isizea,ipartnera
     &     ,itypea,ierrora
#endif

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      ipartnera=ipartner
      itypea=itype
      isizea=isize

      size=MPI_BYTE
      comm=MPI_COMM_WORLD

      CALL MPI_BSEND(buf,isizea,size,ipartnera,itypea,comm,ierrora)
      P_nwrite = ierrora      
#endif
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
