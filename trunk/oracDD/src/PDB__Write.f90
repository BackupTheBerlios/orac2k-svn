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
