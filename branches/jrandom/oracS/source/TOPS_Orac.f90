SUBROUTINE Transform_ORAC(Residue)
  IMPLICIT NONE 
  TYPE(Resid) :: Residue
  INTEGER :: count,n,m
  CHARACTER(len=max_char) :: line
  CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE :: lines

  count=0
  n=0
  DO WHILE(n < SIZE(Residue%line))
     n=n+1
     IF(LEN_TRIM(Residue%line(n)) == 0) THEN
        CYCLE
     END IF
     CALL parse(Residue%line(n), strngs)
     IF(MY_FAM('atoms',strngs(1))) THEN
        n=n+1
        m=0
        DO WHILE(n < SIZE(Residue%line))
           CALL parse(Residue%line(n), strngs)
           IF(MY_FAM('end',strngs(1))) EXIT
           IF(m == 0 .AND. (.NOT. MY_FAM('grou',strngs(1)))) THEN
              count=count+1
           END IF
           count=count+1
           n=n+1
           m=m+1
        END DO
     END IF
     IF(MY_FAM('bonds',strngs(1))) THEN
        n=n+1
        DO WHILE(n < SIZE(Residue%line))
           CALL parse(Residue%line(n), strngs)
           IF(MY_FAM('end',strngs(1))) EXIT
           count=count+1
           n=n+1
        END DO
     END IF
     IF(MY_FAM('imphd',strngs(1))) THEN
        n=n+1
        DO WHILE(n < SIZE(Residue%line))
           CALL parse(Residue%line(n), strngs)
           IF(MY_FAM('end',strngs(1))) EXIT
           count=count+1
           n=n+1
        END DO
     END IF
     IF(MY_FAM('don',strngs(1)) .OR. MY_FAM('acc',strngs(1))) THEN
        count=count+1
     END IF
     IF(MY_FAM('term',strngs(1))) THEN
        IF(TRIM(strngs(2)) /= '*' .AND. TRIM(strngs(3)) /= '*') THEN
           count=count+1
        END IF
     END IF
  END DO
  ALLOCATE(lines(count))
  count=0
  n=0
  DO WHILE(n < SIZE(Residue%line))
     n=n+1
     IF(LEN_TRIM(Residue%line(n)) == 0) THEN
        CYCLE
     END IF
     CALL parse(Residue%line(n), strngs)
     IF(MY_FAM('atoms',strngs(1))) THEN
        n=n+1
        m=0
        DO WHILE(n < SIZE(Residue%line))
           CALL parse(Residue%line(n), strngs)
           IF(MY_FAM('end',strngs(1))) EXIT
           IF(m == 0 .AND. (.NOT. MY_FAM('grou',strngs(1)))) THEN
              count=count+1
              lines(count)='group'
           END IF
           count=count+1
           IF(MY_FAM('grou',strngs(1))) THEN
              lines(count)=TRIM(Residue%line(n))
           ELSE
              line=TRIM(Residue%line(n))
              lines(count)='atom '//TRIM(line)
           END IF
           n=n+1
           m=m+1
        END DO
     END IF
     IF(MY_FAM('bonds',strngs(1))) THEN
        n=n+1
        DO WHILE(n < SIZE(Residue%line))
           CALL parse(Residue%line(n), strngs)
           IF(MY_FAM('end',strngs(1))) EXIT
           count=count+1
           line=TRIM(Residue%line(n))
           lines(count)='bond '//TRIM(line)
           n=n+1
        END DO
     END IF
     IF(MY_FAM('imphd',strngs(1))) THEN
        n=n+1
        DO WHILE(n < SIZE(Residue%line))
           CALL parse(Residue%line(n), strngs)
           IF(MY_FAM('end',strngs(1))) EXIT
           count=count+1
           line=TRIM(Residue%line(n))
           lines(count)='impr '//TRIM(line)
           n=n+1
        END DO
     END IF
     IF(MY_FAM('don',strngs(1))) THEN
        count=count+1
        lines(count)='dono '//TRIM(strngs(3))//' '//TRIM(strngs(2))
     END IF
     IF(MY_FAM('acc',strngs(1))) THEN
        count=count+1
        lines(count)='acce '//TRIM(strngs(3))//' '//TRIM(strngs(2))
     END IF
     IF(MY_FAM('term',strngs(1))) THEN
        IF(TRIM(strngs(2)) /= '*' .AND. TRIM(strngs(3)) /= '*') THEN
           count=count+1
           lines(count)=TRIM(Residue%line(n))
        END IF
     END IF

  END DO

  DEALLOCATE(Residue%line)
  ALLOCATE(Residue%line(count))
  Residue%line=lines
  DEALLOCATE(lines)
  CALL Transform_CHARMM(Residue)
END SUBROUTINE Transform_ORAC
