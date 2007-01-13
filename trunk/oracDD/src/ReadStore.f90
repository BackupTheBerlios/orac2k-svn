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
MODULE ReadStore
  USE Node
  USE Constants
  USE STRPAK
  USE Strings, ONLY: My_Fam
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors,errmsg_f
  IMPLICIT NONE 
  PRIVATE
  PUBLIC ReadStore_, RS__String, ReadStore__Delete
  CHARACTER(len=max_char), ALLOCATABLE, SAVE :: RS__string(:)
CONTAINS
  FUNCTION ReadStore_(line) RESULT(out)
    LOGICAL :: out
    CHARACTER(len=*) :: line
    LOGICAL :: ok
    CHARACTER(len=max_char) :: linea
    CHARACTER(len=max_char), DIMENSION(:), POINTER :: line1=>NULL()
    INTEGER :: n,io, count_out,count_a,iopt

    IF(ALLOCATED(RS__string)) DEALLOCATE(RS__String)
    out=.TRUE.
    CALL Channel(io)
    INQUIRE(file=line,EXIST=ok)
    IF(ok) THEN
       OPEN(unit=io,file=line,form='FORMATTED',status='OLD')
    ELSE
       out=.FALSE.
       errmsg_f='Trying to read file '''//TRIM(line)//''',  which does not exist'
       CALL Add_Errors(-1,errmsg_f)
       RETURN
    END IF

    IF(.NOT. Node_()) STOP
    DO 
       READ(unit=io,fmt='(a80)',IOSTAT=iopt) linea
       IF(iopt /= 0) EXIT
       ok=.FALSE.
       DO n=1,SIZE(Used)
          IF(My_Fam(TRIM(Used(n)),linea)) ok=.TRUE.
       END DO
       IF(.NOT. ok) CYCLE
       CALL Node__Push(linea)
    END DO
    count_out=Node__Size()
    ALLOCATE(RS__string(count_out))
    
    count_a=0
    DO WHILE(Node__Pop(line1))
       count_A=count_A+1
       RS__string(count_a)=ADJUSTL(TRIM(line1(1)))
    END DO
    CALL Node__Delete()
    CLOSE(io)
  END FUNCTION ReadStore_
  SUBROUTINE ReadStore__Delete
    IF(ALLOCATED(RS__string)) DEALLOCATE(RS__String)
  END SUBROUTINE ReadStore__Delete
END MODULE ReadStore
