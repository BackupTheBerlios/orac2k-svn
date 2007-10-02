      SUBROUTINE calc_xslt(anxca,anxbc,anxhe,anxal,anprot,annpro,anpoint
     &     ,protl,wca,whe,wbc,xp0,yp0,zp0,xp1,yp1,zp1,nato
     &     ,errca,errhe,errbc,erral)

************************************************************************
*   Time-stamp: <99/03/02 15:14:13 marchi>                             *
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
     &     ,zp1(*),errca(*),errhe(*),errbc(*),erral(*)
      INTEGER annpro,anpoint(2,*),protl(*),nato
      LOGICAL anprot

*----------------------- VARIABLES IN COMMON --------------------------*

      REAL*8  xyz(3,m1),xyz0(3,m1),xyzfit(3,m1),wca2(m1),whe2(m1)
     &     ,wbc2(m1),work(m1),wal2(m1),qt(4)
      COMMON /rag1/ xyz,xyz0,xyzfit,work,wca2,wbc2,whe2,wal2

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER n,m,l,i,ncount,nc2,nrm,aux2
      REAL*8 aux
      DATA aux2/0/
      
*----------------------- EXECUTABLE STATEMENTS ------------------------*
       IF (aux2 .EQ. 0) THEN
           aux2 = 1
           call zero0(work,m1)
           call zero0(wca2,m1)
           call zero0(wbc2,m1)
           call zero0(whe2,m1)
           call zero0(wal2,m1)
       END IF
!=======================================================================
!----- Calculate XRMS from template: fit all solute atoms --------------
!=======================================================================

      DO i=1,nato
         xyz(1,i)=xp0(i)
         xyz(2,i)=yp0(i)
         xyz(3,i)=zp0(i)
         xyz0(1,i)=xp1(i)
         xyz0(2,i)=yp1(i)
         xyz0(3,i)=zp1(i)
         wca2(i)=wca(i)
         wbc2(i)=wbc(i)
         whe2(i)=whe(i)
         wal2(i)=1.0D0
      END DO

      IF(anxca) THEN
         CALL normal(wca2,nato)
         CALL xfit(xyz,xyz0,xyzfit,qt,wca2,work,nato,aux)
         CALL caldr(xyz,xyzfit,wca2,work,nato)

         IF(anprot) THEN
            DO n=1,annpro
               errca(n) =0.0d0
               nrm = 0
               DO i=anpoint(1,n),anpoint(2,n)
                  errca(n)=errca(n)+work(i)
                  if(work(i) .GT. 1.0d-6) nrm = nrm + 1
               END DO
               if (nrm .ne. 0) then
                  errca(n)= errca(n)/ DFLOAT(nrm)
                  errca(n)= dsqrt(errca(n))
               end if
            END DO
         ELSE
            n=0
            nc2=0
            ncount=1
2000        CONTINUE
            nc2=nc2+1
            errca(nc2) =0.0d0
            nrm = 0
            m=protl(ncount)
            DO l=1,m
               i=protl(ncount+l)
               errca(nc2)=errca(nc2)+work(i)
               if(work(i) .GT. 1.0d-6) nrm = nrm + 1
            END DO
            if (nrm .ne. 0) then
               errca(nc2)= errca(nc2)/ DFLOAT(nrm)
               errca(nc2)= dsqrt(errca(nc2))
            end if
            ncount=ncount+m+1
            n=n+m
            IF(n .LT. nato) GOTO 2000
         END IF

      END IF

      IF(anxbc) THEN
         CALL normal(wbc2,nato)
         CALL xfit(xyz,xyz0,xyzfit,qt,wbc2,work,nato,aux)
         CALL caldr(xyz,xyzfit,wbc2,work,nato)

         IF(anprot) THEN
            DO n=1,annpro
               errbc(n) =0.0d0
               nrm = 0
               DO i=anpoint(1,n),anpoint(2,n)
                  errbc(n)=errbc(n)+work(i)
                  if(work(i) .GT. 1.0d-6) nrm = nrm + 1
               END DO
               if (nrm .ne. 0) then
                  errbc(n)= errbc(n)/ DFLOAT(nrm)
                  errbc(n)= dsqrt(errbc(n))
               end if
            END DO
         ELSE
            n=0
            nc2=0
            ncount=1
2001        CONTINUE
            nc2=nc2+1
            errbc(nc2) =0.0d0
            nrm = 0
            m=protl(ncount)
            DO l=1,m
               i=protl(ncount+l)
               errbc(nc2)=errbc(nc2)+work(i)
               if(work(i) .GT. 1.0d-6) nrm = nrm + 1
            END DO
            if (nrm .ne. 0) then
               errbc(nc2)= errbc(nc2)/ DFLOAT(nrm)
               errbc(nc2)= dsqrt(errbc(nc2))
            end if
            ncount=ncount+m+1
            n=n+m
            IF(n .LT. nato) GOTO 2001
         END IF

      END IF

      IF(anxhe) THEN
         CALL normal(whe2,nato)
         CALL xfit(xyz,xyz0,xyzfit,qt,whe2,work,nato,aux)
         CALL caldr(xyz,xyzfit,whe2,work,nato)

         IF(anprot) THEN
            DO n=1,annpro
               errhe(n) =0.0d0
               nrm = 0
               DO i=anpoint(1,n),anpoint(2,n)
                  errhe(n)=errhe(n)+work(i)
                  if(work(i) .GT. 1.0d-6) nrm = nrm + 1
               END DO
               if (nrm .ne. 0) then
                  errhe(n)= errhe(n)/ DFLOAT(nrm)
                  errhe(n)= dsqrt(errhe(n))
               end if
            END DO
         ELSE
            n=0
            nc2=0
            ncount=1
2002        CONTINUE
            nc2=nc2+1
            errhe(nc2) =0.0d0
            nrm = 0
            m=protl(ncount)
            DO l=1,m
               i=protl(ncount+l)
               errhe(nc2)=errhe(nc2)+work(i)
               if(work(i) .GT. 1.0d-6) nrm = nrm + 1
            END DO
            if (nrm .ne. 0) then
               errhe(nc2)= errhe(nc2)/ DFLOAT(nrm)
               errhe(nc2)= dsqrt(errhe(nc2))
            end if
            ncount=ncount+m+1
            n=n+m
            IF(n .LT. nato) GOTO 2002
         END IF

      END IF

      IF(anxal) THEN
         CALL normal(wal2,nato)
         CALL xfit(xyz,xyz0,xyzfit,qt,wal2,work,nato,aux)
         CALL caldr(xyz,xyzfit,wal2,work,nato)

         IF(anprot) THEN
            DO n=1,annpro
               erral(n) =0.0d0
               nrm = 0
               DO i=anpoint(1,n),anpoint(2,n)
                  erral(n)=erral(n)+work(i)
                  if(work(i) .GT. 1.0d-6) nrm = nrm + 1
               END DO
               if (nrm .ne. 0) then
                  erral(n)= erral(n)/ DFLOAT(nrm)
                  erral(n)= dsqrt(erral(n))
               end if
            END DO
         ELSE
            n=0
            nc2=0
            ncount=1
2003        CONTINUE
            nc2=nc2+1
            erral(nc2) =0.0d0
            nrm = 0
            m=protl(ncount)
            DO l=1,m
               i=protl(ncount+l)
               erral(nc2)=erral(nc2)+work(i)
               if(work(i) .GT. 1.0d-6) nrm = nrm + 1
            END DO
            if (nrm .ne. 0) then
               erral(nc2)= erral(nc2)/ DFLOAT(nrm)
               erral(nc2)= dsqrt(erral(nc2))
            end if
            ncount=ncount+m+1
            n=n+m
            IF(n .LT. nato) GOTO 2003
         END IF

      END IF


*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
