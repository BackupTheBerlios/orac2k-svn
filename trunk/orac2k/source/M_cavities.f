      SUBROUTINE M_cavities(nprocs,mqq)

************************************************************************
*   Time-stamp: <1999-10-12 13:46:20 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed May 19 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nprocs,mqq

*----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'parst.h'
      INCLUDE 'voronoi.h'
      INCLUDE 'analysis.h'

*------------------------- LOCAL VARIABLES ----------------------------*

#if defined SHMEM
#define _memory_ M_memory_s
#else
#define _memory_ M_memory
#endif
      INTEGER naux,M_get_length,mma

*----------------------- EXECUTABLE STATEMENTS ------------------------*

#if defined DYNAMIC_MEM
      IF(nprocs .EQ. 1) THEN
         mma=m1
      ELSE
         mma=m1/(nprocs-1)
      END IF

      naux=M_get_length(maxcav_bin*maxcav_nres,8)
      CALL _memory_(ip_cavity_h,naux)
      
      naux=M_get_length(maxcav_atom*mma,8)
      CALL _memory_(ip_cavity_r,naux)
      mqq=maxcav_atom*mma
      
      naux=M_get_length(2*mma,4)
      CALL _memory_(ip_cavity_n,naux)
#endif

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
