      SUBROUTINE P_split_scalar(node,nprocs,ncube,nbyte,nstart_grp
     &     ,nend_grp,nlocal_grp,nstart,nend,nlocal,ngrp,grppt,nnlpp,cpu
     &     ,ncpu)

************************************************************************
*   Time-stamp: <2004-12-10 19:07:54 marchi>                             *
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

      INTEGER node,nprocs,ncube,nbyte,nstart_grp,nend_grp,nstart,nend
     &     ,nlocal_grp,nlocal,ngrp,ncpu,nnlpp(*),grppt(2,*)
      REAL*8  cpu(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,n,m,j,count,sum,n_s(128),n_e(128),nei,diff,nei0,iret
      REAL*8  factor(128),suma,n_int
      CHARACTER*80 errmsg,tag
      LOGICAL near0
      COMMON /rag1/ factor,n_s,n_e

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      iret=0
      suma=0.0D0
      DO i=1,nprocs
         suma=suma+cpu(i)
      END DO
      IF(near0(suma)) THEN
         DO i=1,nprocs
            factor(i)=0.0D0
         END DO
      ELSE
         suma=suma/DBLE(nprocs)
         DO i=1,nprocs
            factor(i)=(suma-cpu(i))/suma
         END DO
      END IF

      CALL P_expand_i(nnlpp,nstart_grp,nend_grp,nlocal_grp,node
     &     ,nprocs)

      IF(near0(suma)) THEN
         sum=0
         DO i=1,ngrp
            sum=sum+nnlpp(i)
         END DO
         n_int=sum/nprocs
         i=0
         count=0
         DO WHILE (i .LT. ngrp) 
            sum=0
            DO WHILE (sum .LT. n_int .AND. i .LT. ngrp)
               i=i+1
               sum=sum+nnlpp(i)
            END DO
            count=count+1
            n_e(count)=i
         END DO
         IF(node .EQ. 0) THEN
            nstart_grp=1
            nend_grp=n_e(1)
         ELSE
            nstart_grp=n_e(node)+1
            nend_grp=n_e(node+1)
         END IF
         IF(node+1 .EQ. nprocs) nend_grp=ngrp
         nlocal_grp=nend_grp-nstart_grp+1
      ELSE
         DO i=1,nprocs
            n_s(i)=0
            n_e(i)=0
         END DO
         n_s(node+1)=nstart_grp
         n_e(node+1)=nend_grp
         CALL P_merge_veci(n_s,nprocs)
         CALL P_merge_veci(n_e,nprocs)
         
         DO i=1,nprocs
            nei=0
            DO n=n_s(i),n_e(i)
               nei=nei+nnlpp(n)
            END DO
            IF(factor(i) .LT. -0.9D0) factor(i)=-0.5D0
            diff=nei*factor(i)
            nei0=nei
            n=n_e(i)
            IF(nei0+diff .LT. nei) THEN
               DO WHILE(nei .GT. nei0+diff)
                  nei=nei-nnlpp(n)
                  n=n-1
               END DO
            ELSE
               DO WHILE(nei .LT. nei0+diff .AND. n .LT. ngrp)
                  nei=nei+nnlpp(n)
                  n=n+1
               END DO
            END IF
            m=n-n_e(i)
            n_e(i)=n
            DO j=i+1,nprocs
               n_s(j)=n_s(j)+m
               n_e(j)=n_e(j)+m
            END DO
         END DO


         nstart_grp=n_s(node+1)
         nend_grp=n_e(node+1)
         IF(node+1 .EQ. nprocs) nend_grp=ngrp
         nlocal_grp=nend_grp-nstart_grp+1
      END IF

      nstart=grppt(1,nstart_grp)
      nend=grppt(2,nend_grp)
      nlocal=nend-nstart+1
      tag='nbnd group'
      CALL P_check_decomp(nstart_grp,nend_grp,node,nbyte,nprocs,tag,iret
     &     ,errmsg)
      CALL P_get_errmsg(iret,errmsg,80,node,nprocs,ncube,nbyte)

      tag='nbnd atom'
      CALL P_check_decomp(nstart,nend,node,nbyte,nprocs,tag,iret,errmsg)
      CALL P_get_errmsg(iret,errmsg,80,node,nprocs,ncube,nbyte)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
