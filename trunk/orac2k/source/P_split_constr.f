      SUBROUTINE P_split_constr(node,nprocs,ncube,nstart,nend,work2_h
     &     ,iret,errmsg)

************************************************************************
*   Time-stamp: <2009-03-09 12:45:08 marchi>                             *
*                                                                      *
*   Divide intramolecular interaction according to atom decomposition  *
*   Each processor might do a little more work than necessary.         *
*   All intramolecular interactions including any atom in the range    *
*   are included.                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Jul 11 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nstart,nend,node,nprocs,ncube,work2_h(*),iret
      CHARACTER*80 errmsg

*----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'parst.h'
      INCLUDE 'cpropar.h'

      INTEGER i,i_iret

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      iret=0
      errmsg=' '

*=======================================================================
*---- Do bond constraints ----------------------------------------------
*=======================================================================

      CALL P_get_intra(lcnstr,2,lconstr,lconstr_x,nstart,nend,iret)

#if defined PARALLEL
      CALL P_get_iret(iret,node,nprocs,ncube,nbyte)
#endif
      IF(iret .EQ. 1) THEN
         errmsg=
     &  ' Shared constrained bond larger than expected. Change n_h'
     & / /' in MTSMD '
         RETURN
      END IF

*---- Eliminate redundances --------------------------------------------

#if defined PARALLEL
      CALL P_omit_intra(node,nprocs,ncube,nbyte,lconstr_x,work2_h,iret)
#endif
      
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
