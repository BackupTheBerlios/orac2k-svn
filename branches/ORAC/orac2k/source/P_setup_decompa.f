      SUBROUTINE P_setup_decompa(nstart_h,nend_h,nlocal_h,nstart_ah
     &     ,nend_ah,nlocal_ah,nstart_uh,nend_uh,nlocal_uh,ss_index,nbun
     &     ,ngrp,grppt,mres,node,nprocs,ncube,nbyte)

************************************************************************
*   Time-stamp: <2005-01-29 11:42:14 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Fri Mar  5 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nstart_h,nend_h,nlocal_h,nstart_ah,nend_ah,nlocal_ah
     &     ,nstart_uh,nend_uh,nlocal_uh,mres(2,*),nbun
     &     ,node,nprocs,ncube,nbyte,ngrp,ss_index(*),grppt(2,*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,iret,sum,n_s(128),n_e(128),count,n_int
      CHARACTER*80 tag,errmsg

*----------------------- EXECUTABLE STATEMENTS ------------------------*

*=======================================================================
*---- Setup unit  decomposition ----------------------------------------
*=======================================================================
      
#if defined PARALLEL
      iret=0
      sum=0
      DO i=1,nbun
         sum=sum+mres(2,i)-mres(1,i)+1
      END DO
      n_int=sum/nprocs
      i=0
      count=0
      DO WHILE (i .LT. nbun)
         sum=0
         DO WHILE (sum .LT. n_int .AND. i .LT. nbun)
            i=i+1
            sum=sum+mres(2,i)-mres(1,i)+1
         END DO
         count=count+1
         n_e(count)=i
      END DO
      IF(node .EQ. 0) THEN
         nstart_uh=1
         nend_uh=n_e(1)
      ELSE
         nstart_uh=n_e(node)+1
         nend_uh=n_e(node+1)
      END IF
      IF(node+1 .EQ. nprocs) nend_uh=nbun
      nlocal_uh=nend_uh-nstart_uh+1

*=======================================================================
*---- Setup group decomposition ----------------------------------------
*=======================================================================

      nstart_h=mres(1,nstart_uh)
      nend_h=mres(2,nend_uh)
      nlocal_h=nend_h-nstart_h+1

*=======================================================================
*---- Setup atom decomposition -----------------------------------------
*=======================================================================

      nstart_ah=grppt(1,nstart_h)
      nend_ah=grppt(2,nend_h)
      nlocal_ah=nend_ah-nstart_ah+1
      
*=======================================================================
*---- Check decomposition ----------------------------------------------
*=======================================================================
      
      tag=' Unit  '
      CALL P_check_decomp(nstart_uh,nend_uh,node,nbyte,nprocs,tag
     &     ,iret,errmsg)
      CALL P_get_errmsg(iret,errmsg,80,node,nprocs,ncube,nbyte)
      tag=' Group '
      CALL P_check_decomp(nstart_h,nend_h,node,nbyte,nprocs,tag
     &     ,iret,errmsg)
      CALL P_get_errmsg(iret,errmsg,80,node,nprocs,ncube,nbyte)
      tag=' Atom '
      CALL P_check_decomp(nstart_ah,nend_ah,node,nbyte,nprocs,tag
     &     ,iret,errmsg)
      CALL P_get_errmsg(iret,errmsg,80,node,nprocs,ncube,nbyte)

#endif

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
