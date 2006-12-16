  SUBROUTINE Store_Topology(linea)
    IMPLICIT NONE 
    CHARACTER(len=max_char) :: linea
    
    CHARACTER(len=max_char) :: Type, Residue
    LOGICAL :: end_of_file=.FALSE.,exist
    INTEGER :: io
    
    CALL CHANNEL(io)
    ktpg_read=io
    ftpg_read=linea
    INQUIRE(FILE=ftpg_read,EXIST=exist)
    IF(.NOT. exist) THEN
       errmsg_f=TRIM(error_other % g(1))//' : '//TRIM(ftpg_read)//' does not exist'
       CALL Add_Errors(-1,errmsg_f)
       RETURN
    END IF
    OPEN(unit=ktpg_read,file=ftpg_read,form='FORMATTED')
    CALL Create_Keywords('TOPOLOGY')

    CALL ST_Init
    DO WHILE(.NOT. end_of_file)
       CALL ST_TopPar(io,end_of_file)
    END DO
    ALLOCATE(Topology(counter))
    end_of_file=.FALSE.

    CALL ST_Init
    DO WHILE(.NOT. end_of_file)
       CALL ST_TopPar(io,Residue,Type,end_of_file)
       CALL ST_GetTop(Topology,Residue,Type)
    END DO
  END SUBROUTINE Store_Topology
  SUBROUTINE Store_Parameters(linea)
    IMPLICIT NONE 
    CHARACTER(len=max_char) :: linea
    
    CHARACTER(len=max_char) :: Type, Residue
    LOGICAL :: end_of_file=.FALSE.,exist
    INTEGER :: io
    
    CALL CHANNEL(io)
    kpar_read=io
    fpar_read=linea
    INQUIRE(FILE=fpar_read,EXIST=exist)
    IF(.NOT. exist) THEN
       errmsg_f=TRIM(error_other % g(1))//' : '//TRIM(ftpg_read)//' does not exist'
       CALL Add_Errors(-1,errmsg_f)
       RETURN
    END IF
    OPEN(unit=kpar_read,file=fpar_read,form='FORMATTED')
    CALL Create_Keywords('PARAMETERS')
    CALL ST_Init
    DO WHILE(.NOT. end_of_file)
       CALL ST_TopPar(io,end_of_file)
    END DO
    ALLOCATE(Paras(counter))
    end_of_file=.FALSE.
    CALL ST_Init
    DO WHILE(.NOT. end_of_file)
       CALL ST_TopPar(io,Residue,end_of_file)
       CALL ST_GetPar(Paras,Residue)
    END DO
  END SUBROUTINE Store_Parameters
