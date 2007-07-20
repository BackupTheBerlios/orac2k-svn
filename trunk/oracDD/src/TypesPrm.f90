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
MODULE TypesPrm

!!$***********************************************************************
!!$   Time-stamp: <2007-01-11 18:50:30 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Wed Jan 10 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC ----*

  USE Constants
  USE Resid
  USE MyParse
  USE Strings
  USE Errors, ONLY: Add_Errors=>Add, Print_Errors, errmsg_f
  USE Parameters
  USE Print_Defs
  IMPLICIT none
  PRIVATE
  PUBLIC TypesPrm_, TypesPrm__Types, TypesPrm__Number, TypesPrm__Write&
       &, TypesPrm__Read
  CHARACTER(len=max_atm), ALLOCATABLE, TARGET, SAVE :: TypesPrm__Types(:)
CONTAINS
  FUNCTION TypesPrm_() RESULT(out)
    CHARACTER(len=max_atm), POINTER :: out(:)
    INTEGER :: n,ntypes,i_bonds,count_out,count_A,iflags,nlines,nword
    CHARACTER(len=max_char), DIMENSION(:), ALLOCATABLE :: rewrite
    CHARACTER(len=max_char) :: lab0

    out=>NULL()

    IF(.NOT. Types__Valid()) THEN
       errmsg_f='Topology and Parameter files must be read &
            &before force field parameters can be defined'
       CALL Add_Errors(-1,errmsg_f); out=>NULL()
    END IF

    i_Bonds=-1
    DO n=1,SIZE(Topology)
       IF(My_Fxm('HEADER', Topology (n) % Residue )) THEN; i_bonds=n; EXIT; END IF
    END DO

    IF(i_Bonds == -1) THEN
       errmsg_f='Can''t find atom type labels in the input file'
       CALL Add_Errors(-1,errmsg_f)
       out=>NULL()
    END IF
    ntypes=SIZE(Topology(i_bonds) % line)
    count_a=0
    DO n=1,ntypes
       IF(My_FxM('mass',Topology (i_Bonds) % line(n))) THEN
          count_a=count_a+1
       END IF
    END DO
    ALLOCATE(TypesPrm__Types(count_a))
    count_a=0
    DO n=1,ntypes
       IF(My_FxM('mass',Topology (i_Bonds) % line(n))) THEN
          count_a=count_a+1
          nword=MyParse_(Topology (i_Bonds) % line(n))
          IF(nword < 3) THEN
             errmsg_f='Wrong format on the atom type labels file.&
                  & At least four fields are expected. It was found: '//&
                  &TRIM(Topology (i_Bonds) % line(n))
             CALL Add_Errors(-1,errmsg_f)
             out=>NULL()
             RETURN
          END IF
          TypesPrm__Types(count_a) = TRIM(strngs(3))
       END IF
    END DO
    out=>TypesPrm__Types
    WRITE(kprint,*) 'Total atomic types No. =====>',SIZE(out)
  CONTAINS
    FUNCTION Types__Valid() RESULT(out)
      LOGICAL :: out
      
      out=.FALSE.
      IF((ALLOCATED(Topology)) .AND. (ALLOCATED(Paras))) out=.TRUE.
      
    END FUNCTION Types__Valid
  END FUNCTION TypesPrm_
  FUNCTION TypesPrm__Number(label,noerror) RESULT(out)
    CHARACTER(len=*) :: label
    LOGICAL, OPTIONAL :: noerror
    INTEGER  :: out
    CHARACTER(len=max_char) :: lab1
    INTEGER :: i2

    IF(TRIM(label) == 'x' .OR. TRIM(label) == '*') THEN
       out = -1
       RETURN
    END IF
    DO i2=1,SIZE(TypesPrm__Types)
       IF(TRIM(label) == TRIM(TypesPrm__Types(i2))) THEN
          out = i2
          RETURN
       END IF
    END DO
    out = -999
    IF(PRESENT(noerror)) RETURN
    errmsg_f='Atom type '//TRIM(label)//' not on parameter list'
    CALL Add_Errors(-1,errmsg_f)
  END FUNCTION TypesPrm__Number
  SUBROUTINE TypesPrm__Write
    WRITE(kbinary) SIZE(TypesPrm__Types)
    WRITE(kbinary) TypesPrm__Types
  END SUBROUTINE TypesPrm__Write
  FUNCTION TypesPrm__Read() RESULT(out)
    CHARACTER(len=max_atm), POINTER :: out(:)
    INTEGER :: m

    READ(kbinary) m
    ALLOCATE(TypesPrm__Types(m))
    READ(kbinary) TypesPrm__Types
    out=>TypesPrm__Types
  END FUNCTION TypesPrm__Read
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE TypesPrm
