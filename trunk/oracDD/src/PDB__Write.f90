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
SUBROUTINE PDB__Write(PDB__Coords)

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
  CHARACTER(len=3) :: Res
  CHARACTER(len=4) :: AtmName
  REAL(8) :: x,y,z,occ,tmp
  INTEGER :: Serial,ResSeq
  INTEGER :: n
  
  occ=1.0D0; tmp=0.0D0!!$/---------------------------------------------------------------------\
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


  DO n=1,SIZE(PDB__Coords)
     x=PDB__Coords(n) % x
     y=PDB__Coords(n) % y
     z=PDB__Coords(n) % z
     Serial=PDB__Coords(n) % Serial
     AtmName=TRIM(Tpg % atm (Serial) % a % beta)
     AtmName=ADJUSTL(AtmName)
     Res=ADJUSTL(Tpg % atm (Serial) % a % Res)
     ResSeq=Tpg % atm (Serial) % a % Res_No
     CALL TRANUC(Res)
     CALL TRANUC(AtmName)
     WRITE(99,'(A6,I5,2X,A4,A3,2X,I4,4X,3F8.3,2F6.2)') &
          &'ATOM  ',Serial,AtmName,Res,ResSeq,x,y,z,occ,tmp
  END DO
END SUBROUTINE PDB__Write
