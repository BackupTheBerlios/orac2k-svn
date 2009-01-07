MODULE FactorizeNo
  IMPLICIT NONE
  PRIVATE
  PUBLIC FactorizeNo_
  INTEGER, ALLOCATABLE, TARGET, SAVE :: Factors(:)
CONTAINS
  FUNCTION FactorizeNo_(No, N_Maxa) RESULT(out)
    IMPLICIT None
    INTEGER, OPTIONAL :: N_Maxa
    INTEGER :: No
    INTEGER, POINTER :: out(:)

    INTEGER, PARAMETER :: pq_Max=2048
    INTEGER :: p(pq_Max),q(pq_Max)
    INTEGER :: i,j,n,ii
    INTEGER :: m0,count0
    INTEGER, POINTER :: m(:)=>NULL()
    INTEGER :: N_max

    n_Max=pq_Max
    IF(PRESENT(N_Maxa)) N_Max=N_Maxa

    IF(ALLOCATED(Factors)) DEALLOCATE(Factors)
    out=>NULL()

    !First generate all prime numbers less than N_Max - these
    !beginning prime numbers are used for factorization:
    
    II=0
    p=0
    q=1
    DO I=2,40
       IF (q(I)  == 1) THEN
          j=1
          DO WHILE (I*J <= N_Max)
             q(i*j)=0
             j=j+1
          ENDDO
          ii=ii+1
          p(ii)=i
       ENDIF
    ENDDO

!!$  !Next generate and list all prime number factors for each
!!$  !consecutive natural number N:
    
    count0=1
    n=No
    DO I=1,II
       DO WHILE (MOD(N,P(I)).EQ.0)
          m0=p(i)
          count0=count0+1
          n=INT(DBLE(n)/DBLE(p(I)))
       ENDDO
       IF((N /= 1) .AND. (i == ii)) THEN
          count0=count0+1
       END IF
    END DO

    ALLOCATE(Factors(count0))
    m0=1
    count0=1
    Factors(count0)=1
    n=No
    DO I=1,II
       DO WHILE (MOD(N,P(I)).EQ.0)
          m0=p(i)
          count0=count0+1
          Factors(count0)=m0
          n=INT(DBLE(n)/DBLE(p(I)))
       ENDDO
       IF((N /= 1) .AND. (i == ii)) THEN
          count0=count0+1
          Factors(count0)=m0
       END IF
    END DO
    out=>Factors
  END FUNCTION FactorizeNo_
  SUBROUTINE FactorizeNo__Destroy
    DEALLOCATE(Factors)
  END SUBROUTINE FactorizeNo__Destroy
END MODULE FactorizeNo
