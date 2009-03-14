      SUBROUTINE CalcXrmsSecStruct(anxca,anxbc,anxhe,anxal,SecStructure
     &     ,SecStructTotal,SecPointer,grppt,mres,protl,wca,whe,wbc,xp0
     &     ,yp0,zp0,xp1,yp1,zp1,nato,errca,errhe,errbc,erral,drpca,drpbc
     &     ,drphe,drpal)

************************************************************************
*   Time-stamp: <01/02/24 18:33:59 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Jun 29 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none
      INCLUDE  'parst.h'

*----------------------------- ARGUMENTS ------------------------------*

      LOGICAL anxca,anxbc,anxhe,anxal
      REAL*8  wca(*),whe(*),wbc(*),xp0(*),yp0(*),zp0(*),xp1(*),yp1(*)
     &     ,zp1(*),errca(*),errhe(*),errbc(*),erral(*),drpca(*)
     &     ,drpbc(*),drphe(*),drpal(*)
      INTEGER SecStructTotal,SecPointer(*),protl(*),nato,mres(2,*)
     &     ,grppt(2,*)
      LOGICAL SecStructure

*----------------------- VARIABLES IN COMMON --------------------------*

      REAL*8  xyz(3,m1),xyz0(3,m1),xyzfit(3,m1),wca2(m1),whe2(m1)
     &     ,wbc2(m1),work(m1),wal2(m1),qt(4)
      COMMON /rag1/ xyz,xyz0,xyzfit,work,wca2,wbc2,whe2,wal2

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER n,m,l,i,ii,ncount,nc2,m_Sec,n_Sec,NSec,nn

*----------------------- EXECUTABLE STATEMENTS ------------------------*

!=======================================================================
!----- Calculate XRMS from template ------------------------------------
!=======================================================================

      IF(SecStructure) THEN
         NSec=0
         DO n_Sec=1,SecStructTotal
            m_Sec=SecPointer(NSec+1)
            l=0
            DO nn=1,m_Sec
               n=SecPointer(NSec+1+nn)
               DO ii=mres(1,n),mres(2,n)
                  DO i=grppt(1,ii),grppt(2,ii)
                     l=l+1
                     xyz(1,l)=xp0(i)
                     xyz(2,l)=yp0(i)
                     xyz(3,l)=zp0(i)
                     xyz0(1,l)=xp1(i)
                     xyz0(2,l)=yp1(i)
                     xyz0(3,l)=zp1(i)
                     wca2(l)=wca(i)
                     wbc2(l)=wbc(i)
                     whe2(l)=whe(i)
                     wal2(l)=1.0D0
                  END DO
               END DO
            END DO
            m=l
            IF(anxca) THEN
               CALL normal(wca2,m)
               CALL xfit(xyz,xyz0,xyzfit,qt,wca2,work,m,errca(n_Sec))
               CALL caldr(xyz,xyzfit,wca2,work,m)
               l=0
               DO nn=1,m_Sec
                  n=SecPointer(NSec+1+nn)
                  DO ii=mres(1,n),mres(2,n)
                     DO i=grppt(1,ii),grppt(2,ii)
                        l=l+1
                        drpca(i)=drpca(i)+work(l)
                     END DO
                  END DO
               END DO
            END IF
         
            IF(anxbc) THEN
               CALL normal(wbc2,m)
               CALL xfit(xyz,xyz0,xyzfit,qt,wbc2,work,m,errbc(n_Sec))
               CALL caldr(xyz,xyzfit,wbc2,work,m)
               l=0
               DO nn=1,m_Sec
                  n=SecPointer(NSec+1+nn)
                  DO ii=mres(1,n),mres(2,n)
                     DO i=grppt(1,ii),grppt(2,ii)
                        l=l+1
                        drpbc(i)=drpbc(i)+work(l)
                     END DO
                  END DO
               END DO
            END IF

            IF(anxhe) THEN
               CALL normal(whe2,m)
               CALL xfit(xyz,xyz0,xyzfit,qt,whe2,work,m,errhe(n_Sec))
               CALL caldr(xyz,xyzfit,whe2,work,m)
               l=0
               DO nn=1,m_Sec
                  n=SecPointer(NSec+1+nn)
                  DO ii=mres(1,n),mres(2,n)
                     DO i=grppt(1,ii),grppt(2,ii)
                        l=l+1
                        drphe(i)=drphe(i)+work(l)
                     END DO
                  END DO
               END DO
            END IF

            IF(anxal) THEN
               CALL normal(wal2,m)
               CALL xfit(xyz,xyz0,xyzfit,qt,wal2,work,m,erral(n_Sec))
               CALL caldr(xyz,xyzfit,wal2,work,m)
               l=0
               DO nn=1,m_Sec
                  n=SecPointer(NSec+1+nn)
                  DO ii=mres(1,n),mres(2,n)
                     DO i=grppt(1,ii),grppt(2,ii)
                        l=l+1
                        drpal(i)=drpal(i)+work(l)
                     END DO
                  END DO
               END DO
            END IF
            NSec=NSec+m_Sec+1
         END DO
      ELSE
         n=0
         nc2=0
         ncount=1
2000     CONTINUE
         nc2=nc2+1
         m=protl(ncount)
         DO l=1,m
            i=protl(ncount+l)
            xyz(1,l)=xp0(i)
            xyz(2,l)=yp0(i)
            xyz(3,l)=zp0(i)
            xyz0(1,l)=xp1(i)
            xyz0(2,l)=yp1(i)
            xyz0(3,l)=zp1(i)
            wca2(l)=wca(i)
            wbc2(l)=wbc(i)
            whe2(l)=whe(i)
            wal2(l)=1.0D0
         END DO
         IF(anxca) THEN
            CALL normal(wca2,m)
            CALL xfit(xyz,xyz0,xyzfit,qt,wca2,work,m,errca(nc2))
            CALL caldr(xyz,xyzfit,wca2,work,m)
            DO l=1,m
               i=protl(ncount+l)
               drpca(i)=drpca(i)+work(l)
            END DO
         END IF

         IF(anxbc) THEN
            CALL normal(wbc2,m)
            CALL xfit(xyz,xyz0,xyzfit,qt,wbc2,work,m,errbc(nc2))
            CALL caldr(xyz,xyzfit,wbc2,work,m)
            DO l=1,m
               i=protl(ncount+l)
               drpbc(i)=drpbc(i)+work(l)
            END DO
         END IF

         IF(anxhe) THEN
            CALL normal(whe2,m)
            CALL xfit(xyz,xyz0,xyzfit,qt,whe2,work,m,errhe(nc2))
            CALL caldr(xyz,xyzfit,whe2,work,m)
            DO l=1,m
               i=protl(ncount+l)
               drphe(i)=drphe(i)+work(l)
            END DO
         END IF
         
         IF(anxal) THEN
            CALL normal(wal2,m)
            CALL xfit(xyz,xyz0,xyzfit,qt,wal2,work,m,erral(nc2))
            CALL caldr(xyz,xyzfit,wal2,work,m)
            DO l=1,m
               i=protl(ncount+l)
               drpal(i)=drpal(i)+work(l)
            END DO
         END IF

         ncount=ncount+m+1
         n=n+m
         IF(n .LT. nato) GOTO 2000
      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
