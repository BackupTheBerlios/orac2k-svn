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
SUBROUTINE Read_it(knlist,kprint,nsevere,nword,strngs,iret,errmsg&
     &,read_err)
  
!!$======================== DECLARATIONS ================================*

  IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*

  INTEGER :: knlist,kprint,nword,nsevere,iret,read_err
  CHARACTER(80) :: strngs(:),errmsg
  
!!$------------------------- LOCAL VARIABLES ----------------------------*

  CHARACTER(80) :: line
  REAL(8) :: dummy
  LOGICAL :: exist
  
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

!!$
!!$==== Command DENSITY  ================================================
!!$

  line(79:80)='  '
  read_err=0
  FRAGM__=.TRUE.
  DO
     READ(knlist,'(a78)',END=600) line(1:78)
     CALL wrenc(kprint,line)
     IF(line(1:1) .EQ. '#') CYCLE
     CALL parser(line,strngs,nword)
!!$
!!$==== Subcommand PDB_FRAGM ===========================================
!!$
     SELECT CASE(strngs(1))
     CASE DEFAULT
        errmsg=err_unr(2)//strngs(2)//err_end(1:14)//err_end(16:20)
        CALL xerror(errmsg,80,1,30)
        nsevere = nsevere + 1
        
     CASE('print' )
        CALL Read_String(strngs(2),dummy)
        npdb=IDINT(dummy)
        IF(nword == 2) THEN
           filename='PDB_FRAGM.pdb'
           INQUIRE(FILE=filename,EXIST=exist)
           IF(exist) THEN
              CALL openf(kpdb,filename,'FORMATTED','OLD',0)
           ELSE
              CALL openf(kpdb,filename,'FORMATTED','NEW',0)
           END IF
        ELSE IF(nword == 4) THEN
           IF(strngs(3) .EQ. 'OPEN') THEN
              CALL uscrpl(strngs(4),80)
              filename=strngs(4)
              INQUIRE(FILE=filename,EXIST=exist)
              IF(exist) THEN
                 CALL openf(kpdb,strngs(4),'FORMATTED','OLD',0)
              ELSE
                 CALL openf(kpdb,strngs(4),'FORMATTED','NEW',0)
              END IF
           ELSE
              errmsg='OPEN keyword not found'
              CALL xerror(errmsg,80,1,30)
              nsevere = nsevere + 1
           END IF
        ELSE
           errmsg=err_args(1)//'2'
           CALL xerror(errmsg,80,1,30)
           nsevere=nsevere+1
        END IF
!!$
!!$--------------------------------------------------------------------
!!$
       CASE(' ')
          CYCLE

       CASE('END')
          EXIT
       END SELECT
    END DO
    RETURN

600 read_err=1
    RETURN

20  CONTINUE
    iret=1
    errmsg='internal reading error: wrong format?? TAB character??'
    CALL xerror(errmsg,80,1,2)
    
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

  END SUBROUTINE Read_it
