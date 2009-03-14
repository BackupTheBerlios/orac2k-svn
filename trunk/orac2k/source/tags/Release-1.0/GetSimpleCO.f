      SUBROUTINE GetSimpleCO(co,oc)
      IMPLICIT none
      REAL*8 co(3,3),oc(3,3)
      INTEGER i,j

      DO i=1,3
         DO j=1,3
            co(i,j)=0.0D0
            oc(i,j)=0.0D0
         END DO
      END DO

      DO i=1,3
         co(i,i)=1000.0D0
         oc(i,i)=1.0D0/1000.0D0
      END DO

      RETURN
      END
