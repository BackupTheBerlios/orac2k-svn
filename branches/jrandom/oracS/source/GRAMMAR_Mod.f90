MODULE GRAMMAR_Mod
  USE ERROR_Mod, ONLY: Add_error=>Add,Print_Errors
  USE TYPES
  USE Class_Tree
  USE STRINGS_Mod, ONLY: MY_Fam
  USE STRPAK
  IMPLICIT NONE 
  PRIVATE 
  PUBLIC Create,Grammar_Collect, Read_Grammar, Grammar_String
  CHARACTER(len=*), PARAMETER :: lst1='{}',lst3=']['
  LOGICAL, SAVE :: Grammar_Collect=.TRUE.
  CHARACTER(len=max_data), SAVE :: Grammar_String
  CHARACTER(len=max_pars), DIMENSION(:), ALLOCATABLE :: strngs
  CHARACTER(len=max_char_long) :: errmsg_f
CONTAINS
  RECURSIVE SUBROUTINE Create(string,parent)
    IMPLICIT NONE 
    CHARACTER(len=*) :: string,parent
    INTEGER, DIMENSION(:), ALLOCATABLE :: loc
    CHARACTER(len=1) :: a
    INTEGER :: nxt,iflag,m,c,n,q,bgn,stp,n_str
    CHARACTER(len=max_data) :: strngs,str1,str2,str3,str4

    m=LEN_TRIM(string)
    n=0
    c=0
    q=0
    DO WHILE(n < m) 
       n=n+1
       a=string(n:n)
       IF(a == lst1(1:1) .OR. a == lst1(2:2) ) THEN
          IF(a == lst1(1:1)) c=c+1
          IF(c == 1) q=q+1
          IF(a == lst1(2:2)) c=c-1
       END IF
       CYCLE
    END DO
    IF(q == 0) RETURN
    ALLOCATE(loc(q))
    n=0
    c=0
    q=0
    DO WHILE(n < m) 
       n=n+1
       a=string(n:n)
       IF(a == lst1(1:1) .OR. a == lst1(2:2) ) THEN
          IF(a == lst1(1:1)) c=c+1
          IF(c == 1) THEN
             q=q+1; loc(q)=n
          END IF
          IF(a == lst1(2:2)) c=c-1
       END IF
       CYCLE
    END DO
    IF(MOD(SIZE(loc),2) /= 0) THEN
       errmsg_f='Syntax mistake inreading grammar an open { must be followed'&
            &' by }. Instead it was found: '//TRIM(string)
       CALL Add_Error(-1,errmsg_f)
       CALL Print_Errors()
    END IF
    DO n=1,SIZE(loc),2
       bgn=loc(n)+1
       stp=loc(n+1)-1
       strngs=string(bgn:stp)
       str1=strngs
       IF(MY_Fam('{',strngs)) THEN
          str2=str1
          CALL STR_TRIM('{',str1)
          CALL STR_XLET(str2); CALL STR_TRIM('}',str2)
          CALL STR_XLET(str2)
          str1=TRIM(str1)//' '//str2
       END IF
       n_str=1
       IF(grammar_collect) THEN
          CALL get_links(n_str)
          IF(n_str > 1) THEN
             CALL add_node(str3,Parent,(/str4/))
          ELSE
             CALL add_node(str3,Parent)
          END IF
          CALL Create(strngs,TRIM(str3))
       ELSE
          CALL add_node(str1,Parent)
          CALL Create(strngs,TRIM(str1))
       END IF
    END DO
    DEALLOCATE(loc)
  CONTAINS
    SUBROUTINE get_links(n_str)
      INTEGER :: s,n_str
       nxt=1
       c=0
       str1=TRIM(ADJUSTL(str1))
       DO 
          str2=' '
          CALL Token(1,lst3,str1,nxt,str2,iflag)
          IF(iflag /= 0) EXIT
          c=c+1
       END DO
       n_str=c
       nxt=1
       str3=' '
       str2=' '
       CALL Token(1,lst3,str1,nxt,str2,iflag)
       str3=TRIM(str3)//' '//TRIM(str2)
       SELECT CASE(n_str)
       CASE(3)
          str2=' '
          CALL Token(1,lst3,str1,nxt,str2,iflag)
          str4=str2
          str4=TRIM(ADJUSTL(str4))
          str2=' '
          CALL Token(1,lst3,str1,nxt,str2,iflag)
          str3=TRIM(str3)//' '//TRIM(str2)
          str3=TRIM(ADJUSTL(str3))
       CASE(2)
          CALL Token(1,lst3,str1,nxt,str2,iflag)
          str4=str2
          str4=TRIM(ADJUSTL(str4))
       END SELECT

       nxt=1
       CALL Token(1,' ',str3,nxt,str2,iflag)
       str3=TRIM(ADJUSTL(str2))       
     END SUBROUTINE get_links
  END SUBROUTINE Create
  FUNCTION Read_Grammar() RESULT(out)
    LOGICAL :: out

    CHARACTER(len=max_char) :: red
    INTEGER :: io,iopt
    LOGICAL :: exist
    
    CALL CHANNEL(io)
    INQUIRE(file='GRAMMAR.dat',EXIST=exist)
    IF(.NOT. exist) THEN
       errmsg_f='GRAMMAR.dat file does not exist. No documentation&
            & will be avaiable'
       CALL Add_Error(1,errmsg_f)
       out=.FALSE.
       RETURN
    END IF
    OPEN(unit=io,file='GRAMMAR.dat', form='FORMATTED', status='OLD')
    Grammar_String=' '
    DO
       READ(io,'(a)',iostat=iopt) red
       IF(iopt /= 0) EXIT
       Grammar_String=TRIM(Grammar_String)//' '//red
    END DO
    out=.TRUE.
  END FUNCTION Read_Grammar
END MODULE GRAMMAR_Mod
