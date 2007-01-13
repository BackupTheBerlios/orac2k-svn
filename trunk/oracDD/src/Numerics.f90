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
MODULE Numerics

!!$***********************************************************************
!!$   Time-stamp: <2007-01-04 16:53:01 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Wed Dec 20 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program ORAC  ----*
  IMPLICIT none
  PRIVATE
  PUBLIC :: MatInv, Determinant
  REAL(8), SAVE :: Determinant
CONTAINS
  SUBROUTINE MatInv(A,B)
    REAL(8), DIMENSION(:,:) :: A
    REAL(8), DIMENSION(:,:) :: B

    INTEGER, DIMENSION(:), ALLOCATABLE :: l
    INTEGER, DIMENSION(:), ALLOCATABLE :: m
    INTEGER :: o,p
    REAL(8) :: d
    
    o=SIZE(A,1)
    p=SIZE(A,2)
    ALLOCATE(l(o),m(o))
    b=a
    IF(o < 20) THEN
       CALL DMINV(b,o,d,l,m)
    END IF
    IF(d == 0.0D0) b=0.0D0
    Determinant=d
    DEALLOCATE(l,m)
  CONTAINS
    INCLUDE 'dminv.f'
  END SUBROUTINE MatInv
END MODULE Numerics
