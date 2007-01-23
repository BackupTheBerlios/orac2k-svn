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
