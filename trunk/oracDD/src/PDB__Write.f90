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

SUBROUTINE PDB__Write(unit,PDB__Coords)
!!$***********************************************************************
!!$   Time-stamp: <2007-01-14 17:36:13 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Sun Jan 14 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program oracDD ----*
  
  TYPE(AtomPDB) :: PDB__Coords(:)
  INTEGER :: unit
  CHARACTER(len=3) :: Res
  CHARACTER(len=4) :: AtmName
  REAL(8) :: x,y,z,occ,tmp
  INTEGER :: Serial,ResSeq
  INTEGER :: n
  
  occ=1.0D0; tmp=0.0D0

  WRITE(unit,'(a6,3f9.3,3f7.2,1x,a4)')  'CRYST1',a,b,c,alpha,beta,gamma,'P  1'
  DO n=1,SIZE(PDB__Coords)
     x=PDB__Coords(n) % x
     y=PDB__Coords(n) % y
     z=PDB__Coords(n) % z
     Serial=PDB__Coords(n) % Serial
     AtmName=TRIM(ADJUSTL(Tpg % atm (Serial) % a % beta))
     AtmName=ADJUSTR(AtmName)
     Res=ADJUSTL(Tpg % atm (Serial) % a % Res)
     ResSeq=Tpg % atm (Serial) % a % Res_No
     CALL TRANUC(Res)
     CALL TRANUC(AtmName)
     IF(Serial < 100000) THEN
        IF(ResSeq < 10000) THEN
           WRITE(unit,'(A6,I5,1x,A4,1x,A3,2X,I4,4X,3F8.3,2F6.2)') &
                &'ATOM  ',Serial,AtmName,Res,ResSeq,x,y,z,occ,tmp
        ELSE
           WRITE(unit,'(A6,I5,1x,A4,1x,A3,1X,I5,4X,3F8.3,2F6.2)') &
                &'ATOM  ',Serial,AtmName,Res,ResSeq,x,y,z,occ,tmp
        END IF
     ELSE
        IF(ResSeq < 10000) THEN
           WRITE(unit,'(A5,I6,1x,A4,1x,A3,2X,I4,4X,3F8.3,2F6.2)') &
                &'ATOM ',Serial,AtmName,Res,ResSeq,x,y,z,occ,tmp
        ELSE
           WRITE(unit,'(A5,I6,1X,A4,1x,A3,1X,I5,4X,3F8.3,2F6.2)') &
                &'ATOM ',Serial,AtmName,Res,ResSeq,x,y,z,occ,tmp
        END IF
     END IF
  END DO
END SUBROUTINE PDB__Write
