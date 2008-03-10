      SUBROUTINE P_atoms_split_inter(node,nprocs,ncube,nbyte,ntap
     &     ,nstart_2,nend_2,nlocal_2,cnstp,cnst_protp,cnst_protl,iret
     &     ,errmsg)

************************************************************************
*   Time-stamp: <2007-11-27 17:43:31 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Jan 28 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER node,nprocs,ncube,nbyte,ntap,nstart_2,nend_2,nlocal_2
     &     ,cnst_protp,iret,cnstp(2,*),cnst_protl(*)
      CHARACTER*80 errmsg

*----------------------- VARIABLES IN COMMON --------------------------*

      INTEGER n_end(1024)
      COMMON /rag1/ n_end

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,m,count,count1,sum,tot_int,n_int,at_max,k,la,lb

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      n_int=cnst_protp/nprocs

      i=0
      count=0
      count1=0
      DO WHILE (i .LT. cnst_protp) 
         sum=0
         DO WHILE (sum .LT. n_int .AND. i .LT. cnst_protp)
            i=i+1
            m=cnst_protl(count+1)
            sum=sum+1
            count=count+m+1
         END DO
         count=count-m-1
         count1=count1+1
         at_max=-1
         DO j=1,m
            k=cnst_protl(count+1+j)
            la=cnstp(1,k)
            lb=cnstp(2,k)
            IF(la .GT. at_max) at_max=la
            IF(lb .GT. at_max) at_max=lb
         END DO
         count=count+m+1
         n_end(count1)=at_max
      END DO

      IF(node .EQ. 0) THEN
         nstart_2=1
         nend_2=n_end(1)
      ELSE
         nstart_2=n_end(node)+1
         nend_2=n_end(node+1)
      END IF
      IF(node+1 .EQ. nprocs) nend_2=ntap
      nlocal_2=nend_2-nstart_2+1

*=======================================================================
*---- Split bond constraints for non bonded loop -----------------------
*=======================================================================

      CALL P_get_prot(cnst_protp,cnstp,cnst_protl,nstart_2,nend_2
     &     ,iret)
#if defined PARALLEL
      CALL P_get_iret(iret,node,nprocs,ncube,nbyte)
#endif
      IF(iret .EQ. 1) THEN
         errmsg=
     &  ' Shared stretch bond larger than expected. Change n_h'
     & / /' in MTSMD '
         RETURN
      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
