      SUBROUTINE write_hbonds_vor(hbonds_tot,hbonds_res,fstep,khbonds
     &     ,a_mask,d_mask,ss_index,nbun,nres,o1,nato)

************************************************************************
*   Time-stamp: <97/07/10 18:56:36 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed Jul  9 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER o1,nato,nbun,khbonds
      INTEGER ss_index(*),nres(o1,*)
      REAL*8  fstep
      LOGICAL hbonds_tot,hbonds_res,a_mask(*),d_mask(*)

*----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'parst.h'
      INCLUDE 'voronoi.h'
      INTEGER count_res(nores),count(3),index(nores)
      COMMON /rag1/ count_res,index
      
*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,typei,typeij,type_resi,type_resj,m,j,j1,map,n

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO i=1,nbun
         count_res(i)=0
      END DO
      DO i=1,3
         count(i)=0
      END DO

      DO i=1,nato
         typei=ss_index(i)
         type_resi=nres(i,1)
         IF(a_mask(i)) THEN
            m=nnlpp_vor(1,i)
            DO j=1,m
               j1=nnlpp_vor(1+j,i)
               IF(d_mask(j1)) THEN
                  type_resj=nres(j1,1)
                  typeij=typei+ss_index(j1)-1
                  count(typeij)=count(typeij)+1
                  count_res(type_resi)=count_res(type_resi)+1
                  count_res(type_resj)=count_res(type_resj)+1
               END IF
            END DO
         ELSE IF(d_mask(i)) THEN
            m=nnlpp_vor(1,i)
            DO j=1,m
               j1=nnlpp_vor(1+j,i)
               IF(a_mask(j1)) THEN
                  type_resj=nres(j1,1)
                  typeij=typei+ss_index(j1)-1
                  count(typeij)=count(typeij)+1
                  count_res(type_resi)=count_res(type_resi)+1
                  count_res(type_resj)=count_res(type_resj)+1
               END IF
            END DO
         END IF
      END DO
      DO i=1,3
         count(i)=count(i)/2
      END DO
      IF(hbonds_tot)
     &     WRITE(khbonds,'(a1,f11.2,a16,i7,a9,i7,a9,i7)') 'T',fstep,
     &     ' Hbonds Vor. SLT-SLT ',count(1),' SLV-SLV ',count(3)
     &     ,' SLT-SLV ',count(2)
         
      IF(hbonds_res) THEN
         map=0
         DO i=1,nbun
            IF(count_res(i) .GT. 0) THEN
               map=map+1
               index(map)=i
            END IF
         END DO
         DO j=1,map,4
            n=3
            IF(map-j .LT. 3) n=map-j
            WRITE(khbonds,'(a1,f11.2,1x,a8,4(1x,i6,1x,i6))') 'T',fstep,
     &           ' Hbonds Vor. ',(count_res(index(i)),index(i),i=j,j+n)
         END DO
      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
