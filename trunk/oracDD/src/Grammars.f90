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
MODULE Grammars
  USE Errors, ONLY: Add_errors=>Add,Print_Errors
  USE Types
  USE Tree
  USE Strings, ONLY: MY_Fam
  USE STRPAK
  IMPLICIT NONE 
  PRIVATE 

!!$
!!$--- Shared objects
!!$
  
  PUBLIC Grammars_,Grammars__Inputs,Grammars__Docs
  TYPE(trees), POINTER, SAVE :: Grammars__Inputs,Grammars__Docs

!!$
!!$--- Private objects
!!$
  
  CHARACTER(len=max_data), SAVE :: Grammars__String
  CHARACTER(len=*), PARAMETER :: lst1='{}',lst3=']['
  LOGICAL :: Grammar_Collect
  CHARACTER(len=max_pars), DIMENSION(:), ALLOCATABLE :: strngs
  CHARACTER(len=max_char_long) :: errmsg_f
CONTAINS
  SUBROUTINE Grammars_(Inputs__string)
    CHARACTER(len=max_data) :: Inputs__String

    Grammar_Collect=.FALSE.
    CALL Tree__Start(Grammars__Inputs)
    CALL Create(Inputs__String,'')

    IF(Read_Grammar()) THEN
       Grammar_Collect=.TRUE.
       CALL Tree__Start(Grammars__Docs)
       CALL Create(Grammars__String,'')
    END IF

  END SUBROUTINE Grammars_

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
       CALL Add_Errors(-1,errmsg_f)
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
             CALL Tree__add_node(str3,Parent,(/str4/))
          ELSE
             CALL Tree__add_node(str3,Parent)
          END IF
          CALL Create(strngs,TRIM(str3))
       ELSE
          CALL Tree__add_node(str1,Parent)
          CALL Create(strngs,TRIM(str1))
       END IF
    END DO
    DEALLOCATE(loc)
  CONTAINS
    SUBROUTINE get_links(n_str)
      INTEGER :: n_str
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
            & will be available'
       CALL Add_Errors(1,errmsg_f)
       out=.FALSE.
       RETURN
    END IF
    OPEN(unit=io,file='GRAMMAR.dat', form='FORMATTED', status='OLD')
    Grammars__String=' '
    DO
       READ(io,'(a)',iostat=iopt) red
       IF(iopt /= 0) EXIT
       Grammars__String=TRIM(Grammars__String)//' '//red
    END DO
    out=.TRUE.
  END FUNCTION Read_Grammar
END MODULE Grammars
