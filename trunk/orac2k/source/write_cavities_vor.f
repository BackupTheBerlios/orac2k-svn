      SUBROUTINE write_cavities_vor(kcavities,cavities_file,prsymb,mend
     &     ,jmax_cav,bin_size_cav,rmax_size_cav,cavity_h,mq,node,nprocs
     &     ,ncube)

************************************************************************
*   Time-stamp: <04/11/11 10:56:42 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Mon Apr 13 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER kcavities,jmax_cav,mq,node,nprocs,ncube
      INTEGER mend(*)
      REAL*8  cavity_h(mq,*),bin_size_cav,rmax_size_cav
      CHARACTER*80 cavities_file
      CHARACTER*8  prsymb(*)
      
*------------------------- LOCAL VARIABLES ----------------------------*

      INCLUDE 'parst.h'
#if defined SHMEM
#define _memory_ M_memory_s
#define _free_   M_free_s
#else
#define _memory_ M_memory
#define _free_   M_free
#endif
#if defined PARALLEL
#if defined DYNAMIC_MEM
      INTEGER len,M_get_length
      REAL*8 worka(*),workb(mq,*)
      POINTER (ip_worka,worka),(ip_workb,workb)
#else
      REAL*8 worka(nores*nores),workb(nores*nores)
      COMMON /rag1/ worka,workb
#endif
#endif
      INTEGER i,m,nmax,nstart,nend,nlocal,n_int,nntot,GET_IND
      REAL*8  dx,y
      LOGICAL old
      DATA old/.FALSE./
      GET_IND(y)=INT(y)+(SIGN(1.0D0,y-INT(y))-1.0D0)/2+1

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      IF(old) OPEN(unit=kcavities,file=cavities_file,form
     &     ='FORMATTED',status='OLD')
      
      IF(node .EQ. 0) WRITE(kcavities,100)

      dx=bin_size_cav
      nmax=GET_IND(rmax_size_cav/dx)

#if defined PARALLEL
#if defined DYNAMIC_MEM
      len=M_get_length(mq*jmax_cav,8)
      CALL _memory_(ip_worka,len)
      CALL _memory_(ip_workb,len)
#endif
      DO i=1,jmax_cav
         workb(1,i)=cavity_h(1,i)
         DO m=1,mq-1
            workb(1+m,i)=cavity_h(1+m,i)
         END DO
      END DO
      nntot=jmax_cav*mq
      n_int=jmax_cav/nprocs

      nstart=node*n_int*mq+1
      nend=(node+1)*n_int*mq
      IF(node+1 .EQ. nprocs) nend=mq*jmax_cav
      nlocal=nend-nstart+1

      IF(nprocs .NE. 1) THEN
         CALL P_fold_r8(nntot,workb,nstart,nend,nlocal,node,nprocs)
         CALL P_expand_r8(workb,nstart,nend,nlocal,node,nprocs)
      END IF


      IF(node .EQ. 0) THEN
         DO i=1,jmax_cav
            WRITE(kcavities,'(a8,f20.1)') prsymb(mend(i)),workb(1,i)
            DO m=1,nmax
               WRITE(kcavities,'(f12.3,2x,f20.1)') (m-1)*dx,workb(1+m,i)
            END DO
         END DO
      END IF
      CALL P_barrier

#if defined DYNAMIC_MEM
      CALL _free_(ip_worka)
      CALL _free_(ip_workb)
#endif
#else
      DO i=1,jmax_cav
         WRITE(kcavities,'(a8,f20.1)') prsymb(mend(i)),cavity_h(1,i)
         DO m=1,nmax
            WRITE(kcavities,'(f12.3,2x,f20.1)') (m-1)*dx,cavity_h(1+m,i)
         END DO
      END DO
#endif

      old=.TRUE.

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

100   FORMAT('********************************************************',
     &       '***********************'/
     &       '*                                                       ',
     &       '                      *'/
     &       '*                           CAVITIES OUTPUT FILE        ',
     &       '                      *'/
     &       '*                                                       ',
     &       '                      *'/
     &       '*      The cavity distribution function is written for e',
     &       'ach residue.          *'/
     &       '*      The function are written according to the occurre',
     &       'nce of the            *'/
     &       '*      residue in the main sequence                     ',
     &       '                      *'/
     &       '*                                                       ',
     &       '                      *'/
     &       '*      The first line of each distribution function cont',
     &       'ains                  *'/
     &       '*      1st Field :  Residue type                        ',
     &       '                      *'/
     &       '*      2nd Field :  Total number of point averages (grid',
     &       ' point + No. conf.)   *'/
     &       '*                                                       ',
     &       '                      *'/
     &       '*      The histogram lines of each distribution function',
     &       ' contains             *'/
     &       '*      1st Field :  Cavity Radius (AA)                  ',
     &       '                      *'/
     &       '*      2nd Field :  No. of Samples occurring at that cav',
     &       'ity radius            *'/
     &       '*                                                       ',
     &       '                      *'/
     &       '********************************************************',
     &       '***********************'/)
      RETURN
      END
