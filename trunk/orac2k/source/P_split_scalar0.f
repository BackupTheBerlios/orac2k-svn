      SUBROUTINE P_split_scalar0(node,nprocs,nstart_grp,nend_grp
     &     ,nlocal_grp,nstart,nend,nlocal,ngrp,grppt)

************************************************************************
*   Time-stamp: <2007-11-20 15:47:06 marchi>                             *
*                                                                      *
*                                                                      *
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

      INTEGER node,nprocs,nstart_grp,nend_grp,nstart,nend,nlocal_grp
     &     ,nlocal,ngrp,grppt(2,*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,count,sum,n_end(1024),n_int

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      n_int=ngrp/nprocs

      i=0
      count=0
      DO WHILE (i .LT. ngrp) 
         sum=0
         DO WHILE (sum .LT. n_int .AND. i .LT. ngrp)
            i=i+1
            sum=sum+1
         END DO
         count=count+1
         n_end(count)=i
      END DO

      IF(node .EQ. 0) THEN
         nstart_grp=1
         nend_grp=n_end(1)
      ELSE
         nstart_grp=n_end(node)+1
         nend_grp=n_end(node+1)
      END IF
      IF(node+1 .EQ. nprocs) nend_grp=ngrp
      nlocal_grp=nend_grp-nstart_grp+1

      nstart=grppt(1,nstart_grp)
      nend=grppt(2,nend_grp)
      nlocal=nend-nstart+1

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
