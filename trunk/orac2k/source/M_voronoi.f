      SUBROUTINE M_voronoi(nprocs)

************************************************************************
*   Time-stamp: <01/07/09 08:32:27 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu May 13 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nprocs

*----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'parst.h'
      INCLUDE 'voronoi.h'

*------------------------- LOCAL VARIABLES ----------------------------*

#if defined T3E
#define _memory_ M_memory_s
#else
#define _memory_ M_memory
#endif
      INTEGER len,M_get_length,mma

*----------------------- EXECUTABLE STATEMENTS ------------------------*

#if defined DYNAMIC_MEM  
      IF(nprocs .EQ. 1) THEN
         mma=m1
      ELSE
         mma=m1/(nprocs-1)
      END IF
      len=M_get_length(m1*pig_nnl,4)
      CALL _memory_(ip_ig_nnl,len)

      len=M_get_length(mma*pnnlpp_vor,4)
      CALL _memory_(ip_nnlpp_vor,len)

      len=M_get_length(mma*pnnlpp_vor,8)
      CALL _memory_(ip_area_vor,len)

      len=M_get_length(m1,8)
      CALL _memory_(ip_volume_vor,len)
#endif


*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
