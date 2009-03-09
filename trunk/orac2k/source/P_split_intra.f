      SUBROUTINE P_split_intra(node,nprocs,ncube,nstart,nend,work2_h
     &     ,iret,errmsg)

************************************************************************
*   Time-stamp: <2009-03-09 12:45:09 marchi>                             *
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

      INTEGER i,i_iret,mma

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      iret=0
      errmsg=' '
*=======================================================================
*---- Do bond stretching -----------------------------------------------
*=======================================================================

      CALL P_get_intra(lstrtch,2,lstretch,lbnd_x,nstart,nend,iret)
#if defined PARALLEL
      CALL P_get_iret(iret,node,nprocs,ncube,nbyte)
#endif
      IF(iret .EQ. 1) THEN
         errmsg=
     &  ' Shared stretch bond larger than expected. Change n_h'
     & / /' in MTSMD '
         RETURN
      END IF

#if defined PARALLEL
*---- Eliminate redundances --------------------------------------------

      CALL P_omit_intra(node,nprocs,ncube,nbyte,lbnd_x,work2_h,iret)
#endif      
      
*=======================================================================
*---- Do angle bending -------------------------------------------------
*=======================================================================

      CALL P_get_intra(lbndg,3,lbend,lbndg_x,nstart,nend,iret)
#if defined PARALLEL
      CALL P_get_iret(iret,node,nprocs,ncube,nbyte)
#endif      
      IF(iret .EQ. 1) THEN
         errmsg=
     &  ' Shared bending bond larger than expected. Change n_h'
     & / /' in MTSMD '
         RETURN
      END IF

#if defined PARALLEL
*---- Eliminate redundances --------------------------------------------
      
      CALL P_omit_intra(node,nprocs,ncube,nbyte,lbndg_x,work2_h
     &     ,iret)
#endif      

*=======================================================================
*---- Do proper torsion ------------------------------------------------
*=======================================================================

      CALL P_get_intra(ltor,4,ltors,ltor_x,nstart,nend,iret)
#if defined PARALLEL
      CALL P_get_iret(iret,node,nprocs,ncube,nbyte)
#endif      
      IF(iret .EQ. 1) THEN
         errmsg=
     &  ' Shared p-torsion bond larger than expected. Change n_h'
     & / /' in MTSMD '
         RETURN
      END IF

#if defined PARALLEL
*---- Eliminate redundances --------------------------------------------
      
      CALL P_omit_intra(node,nprocs,ncube,nbyte,ltor_x,work2_h
     &     ,iret)
#endif      

*=======================================================================
*---- Do improper torsion ----------------------------------------------
*=======================================================================

      CALL P_get_intra(litr,4,litor,litr_x,nstart,nend,iret)
#if defined PARALLEL
      CALL P_get_iret(iret,node,nprocs,ncube,nbyte)
#endif      
      IF(iret .EQ. 1) THEN
         errmsg=
     &  ' Shared i-torsion bond larger than expected. Change n_h'
     & / /' in MTSMD '
         RETURN
      END IF

#if defined PARALLEL
*---- Eliminate redundances --------------------------------------------
      
      CALL P_omit_intra(node,nprocs,ncube,nbyte,litr_x,work2_h
     &     ,iret)
#endif      

*=======================================================================
*---- Do 1-4 interactions ----------------------------------------------
*=======================================================================

      CALL P_get_intra(int14,2,int14p,int14_x,nstart,nend,iret)
#if defined PARALLEL
      CALL P_get_iret(iret,node,nprocs,ncube,nbyte)
#endif      
      IF(iret .EQ. 1) THEN
         errmsg=
     &  ' Shared 1--4 bond larger than expected. Change n_h'
     & / /' in MTSMD '
         RETURN
      END IF

*---- Eliminate redundances --------------------------------------------
      
#if defined PARALLEL
      CALL P_omit_intra(node,nprocs,ncube,nbyte,int14_x,work2_h
     &     ,iret)
#endif      

*=======================================================================
*---- Do 1-3 interactions ----------------------------------------------
*=======================================================================

      CALL P_get_intra(int13,2,int13p,int13_x,nstart,nend,iret)
#if defined PARALLEL
      CALL P_get_iret(iret,node,nprocs,ncube,nbyte)
#endif      
      IF(iret .EQ. 1) THEN
         errmsg=
     &  ' Shared 1--3 bond larger than expected. Change n_h'
     & / /' in MTSMD '
         RETURN
      END IF

#if defined PARALLEL
*---- Eliminate redundances --------------------------------------------
      CALL P_omit_intra(node,nprocs,ncube,nbyte,int13_x,work2_h
     &     ,iret)
#endif      

*=======================================================================
*---- Do intragroup interactions ---------------------------------------
*=======================================================================

      CALL P_get_intra(ingrp,2,ingrpp,ingrp_x,nstart,nend,iret)
#if defined PARALLEL
      CALL P_get_iret(iret,node,nprocs,ncube,nbyte)
#endif      
      IF(iret .EQ. 1) THEN
         errmsg=
     &  ' Shared 1--3 bond larger than expected. Change n_h'
     & / /' in MTSMD '
         RETURN
      END IF

#if defined PARALLEL
*---- Eliminate redundances --------------------------------------------
      
      CALL P_omit_intra(node,nprocs,ncube,nbyte,ingrp_x,work2_h
     &     ,iret)
#endif      

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
