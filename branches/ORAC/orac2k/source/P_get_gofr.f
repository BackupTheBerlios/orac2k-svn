      SUBROUTINE P_get_gofr(krdf,krdfa,work,maxint,offset,node,nprocs
     &     ,ncube,nbyte)

************************************************************************
*   Time-stamp: <04/11/09 12:20:25 marchi>                             *
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

      INTEGER krdf(*),work(*),krdfa(*),offset,node,ncube,nprocs,nbyte
     &     ,maxint

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER ntot,nstart,nend,nlocal,off

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      off=(offset+3)*maxint
      CALL icopy(off,krdf,1,krdfa,1)
      ntot=(offset+3)*maxint/nprocs
      nstart=node*ntot+1
      nend=(node+1)*ntot
      nlocal=nend-nstart+1

      CALL P_fold_i(krdfa,work,nstart,nend,nlocal,node,nprocs,ncube
     &     ,nbyte,nbyte)

      CALL P_expand_i(krdfa,nstart,nend,nlocal,node,nprocs)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
