!!$/---------------------------------------------------------------------\
!!$                                                                      |
!!$  Copyright (C) 2006-2007 Massimo Marchi <Massimo.Marchi@cea.fr>      |
!!$                                                                      |
!!$      This program is free software;  you  can  redistribute  it      |
!!$      and/or modify it under the terms of the GNU General Public      |
!!$      License version 2 as published  by  the  Free  Software         |
!!$      Foundation;                                                     |
!!$                                                                      |
!!$      This program is distributed in the hope that  it  will  be      |
!!$      useful, but WITHOUT ANY WARRANTY; without even the implied      |
!!$      warranty of MERCHANTABILITY or FITNESS  FOR  A  PARTICULAR      |
!!$      PURPOSE.   See  the  GNU  General  Public License for more      |
!!$      details.                                                        |
!!$                                                                      |
!!$      You should have received a copy of the GNU General  Public      |
!!$      License along with this program; if not, write to the Free      |
!!$      Software Foundation, Inc., 59  Temple  Place,  Suite  330,      |
!!$      Boston, MA  02111-1307  USA                                     |
!!$                                                                      |
!!$\---------------------------------------------------------------------/
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
