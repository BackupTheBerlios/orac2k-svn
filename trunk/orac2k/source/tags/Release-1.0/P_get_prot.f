      SUBROUTINE P_get_prot(nprot,cnst,cnst_protl,nstart,nend,iret)

************************************************************************
*   Time-stamp: <98/07/14 15:40:23 marchi>                             *
*                                                                      *
*   Get for each processor the right number og groups                  *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Mon Jul 13 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER iret,nprot,nstart,nend,cnst(2,*),cnst_protl(*)
      
*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,k,ka,la,lb,cnstp,count1,count2,count3
      LOGICAL ok

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      count1=0
      count2=0
      count3=0
      DO i=1,nprot
         cnstp=cnst_protl(1+count1)
         ok=.FALSE.
         DO ka=1,cnstp
            k=cnst_protl(1+count1+ka)
            la=cnst(1,k)
            lb=cnst(2,k)
            IF(la .GT. nstart .AND. la .LT. nend .OR. lb .GT. nstart
     &           .AND. lb .LT. nend) ok=.TRUE.
         END DO
         IF(ok) THEN
            count3=count3+1
            cnst_protl(1+count2)=cnstp
            DO ka=1,cnstp
               cnst_protl(1+count2+ka)=cnst_protl(1+count1+ka)
            END DO
            count2=count2+cnstp+1
         END IF
         count1=count1+cnstp+1
      END DO
      nprot=count3

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
