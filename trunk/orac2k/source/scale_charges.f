      SUBROUTINE scale_charges(kprint,nprot_charges,prot_charges,chrge
     &     ,protl,nprot,res,scharge,UnCharge)

************************************************************************
*   Time-stamp: <01/02/24 16:28:17 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Jan 30 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER kprint,nprot,nprot_charges,protl(*),prot_charges(*),res(*)
      REAL*8  chrge(*)
      LOGICAL scharge,UnCharge

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,j1,j2,ncount,noff,natomm,Nstart,Nend,m_old,m_new,nato
      LOGICAL Near0
      REAL*8  sum,sum1

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      natomm=0
      ncount=0
      nato=0
      sum1=0.0D0
      DO i=1,nprot
         noff=protl(ncount+1)
         sum=0.0D0
         DO j=ncount+2,ncount+1+noff
            j1=protl(j)
            sum=sum+chrge(j1)
         END DO
         IF(scharge) THEN
            DO j=1,nprot_charges
               IF(i .EQ. prot_charges(j)) THEN
                  sum1=sum1+sum
                  natomm=natomm+noff
               END IF
            END DO
         END IF
         IF(DABS(sum) .GT. 1.0D-6) THEN
            WRITE(kprint,'(5x,''Found Charge'',f10.5,
     &           '' on Solute Molecule '',i6)') sum,i
         END IF
         ncount=ncount+noff+1
         nato=nato+noff
      END DO
      IF(UnCharge) THEN
         m_old=0
         Nstart=0
         DO i=1,nato
            m_new=res(i)
            IF(m_old .NE. m_new) THEN
               Nend=i-1
               IF(Nend .NE. 0 .AND. Nend-Nstart+1 .GT. 3) THEN
                  sum=0.0D0
                  DO j=Nstart,Nend
                     sum=sum+chrge(j)
                  END DO
                  IF(.NOT. Near0(sum)) sum=sum/DBLE(Nend-Nstart+1)
                  DO j=Nstart,Nend
                     chrge(j)=chrge(j)-sum
                  END DO
               END IF
               Nstart=i
               m_old=m_new
            END IF
         END DO
         WRITE(kprint,'(5x,a)')
     &        'Residue UnCharging Done               ---->'
         natomm=0
         ncount=0
         sum1=0.0D0
         DO i=1,nprot
            noff=protl(ncount+1)
            sum=0.0D0
            DO j=ncount+2,ncount+1+noff
               j1=protl(j)
               sum=sum+chrge(j1)
            END DO
            IF(DABS(sum) .GT. 1.0D-6) THEN
               WRITE(kprint,'(5x,''Found New Charge'',f10.5,
     &              '' on Solute Molecule '',i6)') sum,i
            END IF
            ncount=ncount+noff+1
            nato=nato+noff
         END DO
      END IF
      IF(scharge) THEN
         sum1=sum1/DBLE(natomm)
         ncount=0
         DO i=1,nprot
            noff=protl(ncount+1)
            DO j=1,nprot_charges
               IF(i .EQ. prot_charges(j)) THEN
                  DO j1=ncount+2,ncount+1+noff
                     j2=protl(j1)
                     chrge(j2)=chrge(j2)-sum1
                  END DO
               END IF
            END DO
            ncount=ncount+noff+1
         END DO

        
         WRITE(kprint,'(5x,a)')
     &        'Charge Scaling Done               ---->'
         ncount=0
         sum1=0.0D0
         DO i=1,nprot
            noff=protl(ncount+1)
            sum=0.0D0
            DO j=ncount+2,ncount+1+noff
               j1=protl(j)
               sum=sum+chrge(j1)
            END DO
            DO j=1,nprot_charges
               IF(i .EQ. prot_charges(j)) THEN
                  sum1=sum1+sum
                  natomm=natomm+noff
               END IF
            END DO
            IF(DABS(sum) .GT. 1.0D-6) THEN
               WRITE(kprint,'(5x,''Found New Charge'',f10.5,
     &              '' on Solute Molecule '',i6)') sum,i
            END IF
            ncount=ncount+noff+1
         END DO

      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
