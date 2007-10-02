      SUBROUTINE P_omit_intra(node,nprocs,ncube,nbyte,intra_x,work,iret)

************************************************************************
*   Time-stamp: <04/11/11 17:43:46 marchi>                             *
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

      INTEGER node,nprocs,ncube,iret,nbyte,intra_x(*),work(*)

*----------------------- VARIABLES IN COMMON --------------------------*

      include 'mpif.h'
      integer status(mpi_status_size)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,null,idest,isrc,ierr,n,m,P_nread,P_nwrite,off,k,j,l,na
     &     ,nb
      LOGICAL ok

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      CALL P_Buffer_attach

      n=(intra_x(1)+1)*nbyte
      na=intra_x(1)
      DO i=node+1,nprocs-1
         idest=i
         ierr=P_nwrite(n,nbyte,idest,0,null)
         ierr=P_nwrite(intra_x,n,idest,1,null)
      END DO
      DO i=0,node-1
         isrc=i
         ierr=P_nread(m,nbyte,isrc,0,null)
         ierr=P_nread(work,m,isrc,1,null)
         nb=work(1)

*---- Set to -intra_x all elements equal to those in work

         DO j=1,na
            DO k=1,nb
               IF(work(1+k) .EQ. intra_x(1+j)) THEN
                  intra_x(1+j)=-intra_x(1+j)
               END IF
            END DO
         END DO
      END DO
      CALL P_Buffer_detach

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
