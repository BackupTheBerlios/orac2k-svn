!!$/---------------------------------------------------------------------\
!!$   Copyright  © 2006-2007 Massimo Marchi <Massimo.Marchi at cea.fr>   |
!!$                                                                      |
!!$    This software is a computer program named oracDD whose            |
!!$    purpose is to simulate and model complex molecular systems.       |
!!$    The code is written in fortran 95 compliant with Technical        |
!!$    Report TR 15581, and uses MPI-1 routines for parallel             |
!!$    coding.                                                           |
!!$                                                                      |
!!$    This software is governed by the CeCILL license under             |
!!$    French law and abiding by the rules of distribution of            |
!!$    free software.  You can  use, modify and/ or redistribute         |
!!$    the software under the terms of the CeCILL icense as              |
!!$    circulated by CEA, CNRS and INRIA at the following URL            |
!!$    "http://www.cecill.info".                                         |
!!$                                                                      |
!!$    As a counterpart to the access to the source code and rights      |
!!$    to copy, modify and redistribute granted by the license,          |
!!$    users are provided only with a limited warranty and the           |
!!$    software's author, the holder of the economic rights, and         |
!!$    the successive licensors have only limited liability.             |
!!$                                                                      |
!!$    The fact that you are presently reading this means that you       |
!!$    have had knowledge of the CeCILL license and that you accept      |
!!$    its terms.                                                        |
!!$                                                                      |
!!$    You should have received a copy of the CeCILL license along       |
!!$    with this program; if not, you can collect copies on the URL's    |
!!$    "http://www.cecill.info/licences/Licence_CeCILL_V2-en.html"       |
!!$    "http://www.cecill.info/licences/Licence_CeCILL_V2-fr.html"       |
!!$                                                                      |
!!$----------------------------------------------------------------------/
    SUBROUTINE Orac_
      INTEGER :: n,m,nword,c
      CHARACTER(len=max_char) :: line,aux
      CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE :: lines
      LOGICAL :: ok

      IF(.NOT. Node_()) STOP
      ok=.FALSE.
      DO n=1,SIZE(Ri % line)
         IF(LEN_TRIM(Ri % line(n)) == 0) CYCLE
         nword=Myparse_(Ri % line(n))
         IF(MY_FAM('atoms',strngs(1))) THEN
            nword=Myparse_(Ri % line(n+1))
            IF(.NOT. MY_FAM('grou',strngs(1))) THEN
               line='group'
               CALL Node__Push(line)
            END IF
            ok=.TRUE.
         ELSE IF(ifend(strngs(1)) .AND. ok) THEN
            EXIT
         ELSE IF(MY_FAM('group',strngs(1)) .AND. ok) THEN
            line='group'
            CALL Node__Push(line)
         ELSE IF(ok) THEN
            aux='ATOM'
            line=ST_Concat(1,strngs)
            line=TRIM(aux)//' '//TRIM(line)
            CALL Node__Push(line)
         END IF
      END DO

      ok=.FALSE.
      DO n=1,SIZE(Ri % line)
         IF(LEN_TRIM(Ri % line(n)) == 0) CYCLE
         nword=Myparse_(Ri % line(n))
         IF(MY_Fxm('bonds',strngs(1))) THEN
            aux='BOND'
            ok=.TRUE.
         ELSE IF(ifend(strngs(1)) .AND. ok) THEN
            EXIT
         ELSE IF(ok) THEN
            line=ST_Concat(1,strngs)
            line=TRIM(aux)//' '//TRIM(line)
            CALL Node__Push(line)
         END IF
      END DO

      ok=.FALSE.
      DO n=1,SIZE(Ri % line)
         IF(LEN_TRIM(Ri % line(n)) == 0) CYCLE
         nword=Myparse_(Ri % line(n))
         IF(MY_Fxm('imphd',strngs(1))) THEN
            aux='IMPR'
            ok=.TRUE.
         ELSE IF(ifend(strngs(1)) .AND. ok) THEN
            EXIT
         ELSE IF(ok) THEN
            line=ST_Concat(1,strngs)
            line=TRIM(aux)//' '//TRIM(line)
            CALL Node__Push(line)
         END IF
      END DO
      DO n=1,SIZE(Ri % line)
         IF(LEN_TRIM(Ri % line(n)) == 0) CYCLE
         nword=Myparse_(Ri % line(n))
         IF(MY_Fxm('acc',strngs(1))) THEN
            aux='ACCE'
            line=TRIM(aux)//' '//TRIM(strngs(3))//' '//TRIM(strngs(2))
            CALL Node__Push(line)
         ELSE IF(MY_Fxm('don',strngs(1))) THEN
            aux='DONO'
            line=TRIM(aux)//' '//TRIM(strngs(3))//' '//TRIM(strngs(2))
            CALL Node__Push(line)
         ELSE IF(MY_Fxm('termatom',strngs(1))) THEN
            aux='TERM'
            IF(TRIM(strngs(2)) /= '*' .AND. TRIM(strngs(3)) /= '*') THEN
               line=TRIM(aux)//' '//TRIM(strngs(2))//' '//TRIM(strngs(3))
               CALL Node__Push(line)
            END IF
         END IF
      END DO
      n=Node__Size()
      DEALLOCATE(Residue(i_L) % line)
      ALLOCATE(Residue(i_L) % line(n))
      n=0
      DO WHILE(Node__Pop(line))
         n=n+1
         CALL TRANLC(line)
         Residue(i_L) % line(n) = TRIM(line)
      END DO
    END SUBROUTINE Orac_
    FUNCTION ifend(line) RESULT(out)
      CHARACTER(len=max_char) :: line
      LOGICAL :: out
      
      INTEGER :: n,c
      
      out=.FALSE.
      c=0
      DO n=1,SIZE(Oracs)
         IF(MY_Fxm(Oracs(n),strngs(1))) c=c+1
      END DO
      IF(c /= 0) out=.TRUE.
    END FUNCTION ifend
