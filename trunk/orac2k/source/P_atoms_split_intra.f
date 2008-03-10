      SUBROUTINE P_atoms_split_intra(node,nprocs,ncube,nbyte,nstart_1
     &     ,nend_1,nlocal_1,nstart_ex0,nend_ex0,nlocal_ex0,nstart_ex
     &     ,nend_ex,nlocal_ex,nprot,protl,cnstp,cnst_protp,cnst_protl
     &     ,cnst_protp_1,cnst_protl_1,cnst_protp_ex,cnst_protl_ex,ntap
     &     ,lstrtch,lstretch,lbndg,lbend,atomp,iret,errmsg)

************************************************************************
*   Time-stamp: <2008-03-10 14:44:51 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Jan 19 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER node,nprocs,ncube,nbyte,nstart_1,nend_1,nlocal_1,nprot
     &     ,nstart_ex0,nend_ex0,nlocal_ex0,nstart_ex,nend_ex,nlocal_ex
     &     ,cnst_protp,cnst_protp_1,cnst_protp_ex,ntap,iret,lstretch
     &     ,lbend,cnstp(2,*),cnst_protl(*),cnst_protl_1(*)
     &     ,cnst_protl_ex(*),lstrtch(2,*),lbndg(3,*),atomp(*),protl(*)
      CHARACTER*80 errmsg

*----------------------- VARIABLES IN COMMON --------------------------*

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,m,at_max,at_min,nato_c,count,count1,sum_n,sum_o,la,lb
     &     ,lc,k,ac_min,ac_max,c_min,c_max,nstart_c,nend_c,c_proc
     &     ,get_pointer_protl,aux
      CHARACTER*80 tag
      LOGICAL, DIMENSION(:), ALLOCATABLE :: mask
      INTEGER, DIMENSION(:), ALLOCATABLE :: p_max,p_min,nstart_d,nend_d
      LOGICAL ok

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      ALLOCATE(mask(nprot))
      ALLOCATE(p_max(nprot),p_min(nprot))
      ALLOCATE(nstart_d(0:nprocs),nend_d(0:nprocs))
      p_max=-1
      p_min=1E+8

      count=0
      DO i=1,nprot
         m=protl(count+1)
         DO k=1,m
            j=protl(count+1+k)
            IF(p_max(i) < j) p_max(i)=j
            IF(p_min(i) > j) p_min(i)=j
         END DO
         count=count+1+m
      END DO

      mask=.FALSE.
      DO i=1,nprot
         IF(i /= 1 .AND. i /= nprot) THEN
            IF(p_max(i-1) < p_min(i) .AND. p_min(i+1) > p_max(i)) THEN
               mask(i)=.TRUE.
            END IF
         ELSE
            mask(i)=.TRUE.
         END IF
      END DO

*=======================================================================
*---- Split atoms to compute intramolecular interactions ---------------
*=======================================================================

      at_min=1e8
      at_max=-1
      DO i=1,lstretch
         la=lstrtch(1,i)
         lb=lstrtch(2,i)
         IF(la .LT. at_min) at_min=la
         IF(lb .LT. at_min) at_min=lb
         IF(la .GT. at_max) at_max=la
         IF(lb .GT. at_max) at_max=lb
      END DO
      DO i=1,lbend
         la=lbndg(1,i)
         lb=lbndg(2,i)
         lc=lbndg(3,i)
         IF(la .LT. at_min) at_min=la
         IF(lb .LT. at_min) at_min=lb
         IF(lc .LT. at_min) at_min=lc
         IF(la .GT. at_max) at_max=la
         IF(lb .GT. at_max) at_max=lb
         IF(lc .GT. at_max) at_max=lc
      END DO
<<<<<<< P_atoms_split_intra.f

      IF(at_max .NE. 0 .AND. at_min .NE. 1e8) THEN
=======
      IF(at_max .NE. -1 .AND. at_min .NE. 1e8) THEN
>>>>>>> 1.4
         nato_c=(atomp(at_max)-atomp(at_min)+1)/nprocs
         DO i=0,nprocs-1
            nstart_d(i)=nato_c*i+1
            nend_d(i)=nato_c*(i+1)
         END DO
         nend_d(nprocs-1)=atomp(at_max)

         DO i=1,nprocs-1
            ok=.FALSE.
            count=0
            DO WHILE(.NOT.ok)
               IF(mask(nend_d(i-1)) .AND. mask(nstart_d(i))) THEN
                  ok=.TRUE.
               ELSE IF(count < 40) THEN
                  nend_d(i-1)=nend_d(i-1)+1
                  nstart_d(i)=nstart_d(i)+1
                  count=count+1
               ELSE
                  ok=.TRUE.
               END IF
            END DO
         END DO
         nstart_c=nstart_d(node)
         nend_c=nend_d(node)
      ELSE
         nato_c=0
         nstart_c=1
         nend_c=0
      END IF

      count=0
      DO i=1,nstart_c-1
         m=protl(count+1)
         count=count+1+m
      END DO
      ac_min=1e8
      ac_max=-1
      DO i=nstart_c,nend_c
         m=protl(count+1)
         DO k=1,m
            j=protl(count+1+k)
            IF(j .LT. ac_min) ac_min=j
            IF(j .GT. ac_max) ac_max=j
         END DO
         count=count+m+1
      END DO
      
      IF(ac_min .NE. 1e8 .AND. ac_max .NE. -1) THEN
         nstart_1=ac_min
         nend_1=ac_max
      ELSE
         nstart_1=1
         nend_1=0
      END IF
      nlocal_1=nend_1-nstart_1+1
         
#if defined PARALLEL
      tag='bnd atom'
      CALL P_check_decomp(nstart_1,nend_1,node,nbyte,nprocs,tag,iret
     &     ,errmsg)
      CALL P_get_errmsg(iret,errmsg,80,node,nprocs,ncube,nbyte)
#endif

      DEALLOCATE(mask)
      DEALLOCATE(p_max,p_min)
      DEALLOCATE(nstart_d,nend_d)

*=======================================================================
*---- Set nstart and nend for atoms non involved in intramolecular -----
*---- interactions -----------------------------------------------------
*=======================================================================

      IF(at_max .EQ. -1) THEN
         nstart_ex=1
         nend_ex=ntap
         nlocal_ex=ntap
      ELSE
         IF(ntap .NE. at_max) THEN
            nato_c=(atomp(ntap)-atomp(at_max+1)+1)/nprocs
         
            nstart_c=nato_c*node+atomp(at_max)+1
            nend_c=nato_c*(node+1)+atomp(at_max)
            IF(node+1 .EQ. nprocs) nend_c=atomp(ntap)
         
            count=0
            DO i=1,nstart_c-1
               m=protl(count+1)
               count=count+1+m
            END DO
            ac_min=1e8
            ac_max=-1
            DO i=nstart_c,nend_c
               m=protl(count+1)
               DO k=1,m
                  j=protl(count+1+k)
                  IF(j .LT. ac_min) ac_min=j
                  IF(j .GT. ac_max) ac_max=j
               END DO
               count=count+m+1
            END DO
         
            nstart_ex=ac_min
            nend_ex=ac_max
            nlocal_ex=nend_ex-nstart_ex+1
         
         ELSE
            nstart_ex=ntap+1
            nend_ex=0
            nlocal_ex=0
         END IF
      END IF
      nstart_ex0=nstart_ex
      nend_ex0=nend_ex
      nlocal_ex0=nlocal_ex

#if defined PARALLEL
      IF(ntap .NE. at_max .AND. at_max .NE. -1) THEN
         tag='ext atom'
         CALL P_check_decomp(nstart_ex,nend_ex,node,nbyte,nprocs,tag
     &        ,iret,errmsg)
         CALL P_get_errmsg(iret,errmsg,80,node,nprocs,ncube,nbyte)
      END IF
#endif

*=======================================================================
*---- Set cnst_protl_1 and cnst_protl_ex equal to cnst_protl -----------
*=======================================================================

      cnst_protp_1=cnst_protp
      cnst_protp_ex=cnst_protp
      count=0
      DO i=1,cnst_protp
         m=cnst_protl(1+count)
         count=count+m+1
      END DO

      DO i=1,count
         cnst_protl_1(i)=cnst_protl(i)
         cnst_protl_ex(i)=cnst_protl(i)
      END DO

*=======================================================================
*---- Split bond constraints for bonded loop ---------------------------
*=======================================================================

      CALL P_get_prot(cnst_protp_1,cnstp,cnst_protl_1,nstart_1,nend_1
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

*=======================================================================
*---- Split bond constraints for atoms nont involved in intramolecular -
*---- interactions -----------------------------------------------------
*=======================================================================

      CALL P_get_prot(cnst_protp_ex,cnstp,cnst_protl_ex,nstart_ex
     &     ,nend_ex,iret)
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
