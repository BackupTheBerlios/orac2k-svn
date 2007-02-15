MODULE FactorizeNo
  USE Node
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
    
    IF(.NOT. Node_()) STOP
    m0=1
    CALL Node__Push(m0)
    n=No
    DO I=1,II
       DO WHILE (MOD(N,P(I)).EQ.0)
          m0=p(i)
          CALL Node__Push(m0)
          n=INT(DBLE(n)/DBLE(p(I)))
       ENDDO
       IF((N /= 1) .AND. (i == ii)) THEN
          CALL Node__Push(n)
       END IF
    END DO
    count0=Node__Size()
    ALLOCATE(Factors(count0))
    count0=0
    DO WHILE(Node__Pop(m))
       count0=count0+1
       Factors(count0)=m(1)
    END DO
    out=>Factors
    CALL Node__Delete()
  END FUNCTION FactorizeNo_
  SUBROUTINE FactorizeNo__Destroy
    DEALLOCATE(Factors)
  END SUBROUTINE FactorizeNo__Destroy
END MODULE FactorizeNo
