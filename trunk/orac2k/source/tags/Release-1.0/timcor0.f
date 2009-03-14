      SUBROUTINE timcor0(natob,nintv)
      IMPLICIT none
      INTEGER natob,nintv

      INCLUDE 'vacf.h'
      INTEGER n,m,i,nsevere
      CHARACTER*80 errmsg


      nsevere=0
      IF(nintv .GT. nc) THEN
          errmsg='Physical dimensions of veloci'
     x//'ty array exceeded. Redimension _NUM_CORR_ in config.h.'
          CALL xerror(errmsg,80,1,30)
          nsevere=nsevere+1
      END IF 
      IF(natob .GT. nb) THEN
          errmsg='Physical dimensions of veloci'
     x//'ty array exceeded. Redimension _NAT_SOLU_ in config.h'
          CALL xerror(errmsg,80,1,30)
          nsevere=nsevere+1
      END IF
      IF(nsevere .GT. 0) THEN
          CALL xerror(errmsg,80,1,2)
       END IF

       DO n=1,nintv
          vafp(n-1,1)=0.0D0 
          weip(n-1,1)=0.0D0
          vafp(n-1,2)=0.0D0 
          weip(n-1,2)=0.0D0
          DO m=1,natob
             DO i=1,3
                vpl(i,m,n)=0.0D0
             END DO
          END DO
       END DO

      RETURN
      END
