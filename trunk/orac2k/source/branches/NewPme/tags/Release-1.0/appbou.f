      SUBROUTINE appbou(xp0,yp0,zp0,xpg,ypg,zpg,pmass,nstart,nend,grppt)

************************************************************************
*                                                                      *
*     Recalculate positions of the groups                              *
*                                                                      *
*     XP0     :  Protein coordinates in simulation box units.   (I/O)  *
*     YP0        >> REAL*8 XP0(*) <<                                   *
*     ZP0                                                              *
*     XPG     :  Coordinates of the groups in simulation box    (I/O)  *
*     YPG        units.                                                *
*     ZPG        >> REAL*8 XPG(*) <<                                   *
*     PMASS   :  List of atomic masses.                           (I)  *
*                >> REAL*8 PMASS(*) <<                                 *
*     NGRP    :  Number of groups in the protein                  (I)  *
*     GRPPT   :  Pointer to the first an last atom of a group     (I)  *
*                >> INTEGER GRPPT(2,*) <<                              *
*     RGGRUP  :  If SECTION is .TRUE. this is the number of       (I)  *
*                groups which are contained in the primary             *
*                region                                                *
*     RGGMAP  :  If SECTION is .TRUE. this is the list of atom    (I)  *
*                numbers contained in the primary region               *
*                >> INTEGER RGGMAP(*) <<                               *
*     SCTION  :  Flag which indicates if a primary and secondary  (I)  *
*                simulation regions have been defined                  *
*                >> LOGICAL SCTION <<                                  *
*                                                                      *
*---- Last update  25/1092 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi CECAM, Orsay France 1992               *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER nstart,nend,grppt(2,*)
      REAL*8  pmass(*),xp0(*),yp0(*),zp0(*),xpg(*),ypg(*),zpg(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8  sumx,sumy,sumz
      INTEGER n,m

*==================== EXECUTABLE STATEMENTS ============================

      DO n=nstart,nend
         sumx=0.0D0
         sumy=0.0D0
         sumz=0.0D0
         DO m=grppt(1,n),grppt(2,n)
            sumx=sumx+pmass(m)*xp0(m)
            sumy=sumy+pmass(m)*yp0(m)
            sumz=sumz+pmass(m)*zp0(m)
         END DO
         xpg(n)=sumx
         ypg(n)=sumy
         zpg(n)=sumz
      END DO

*==================== END OF EXECUTABLE STATEMENTS =====================

      RETURN
      END
