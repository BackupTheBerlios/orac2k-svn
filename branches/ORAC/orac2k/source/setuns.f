      SUBROUTINE setuns(time)

************************************************************************
*                                                                      *
*     Set up the unit system used in the program. Values are           *
*     assigned to some fundamental physical constants. All             *
*     these variables are left untouched throughout the program.       *
*                                                                      *
*--------------------- ARGUMENT ---------------------------------------*
*                                                                      *
*     TIME    :  Timestep size.                                        *
*                                                                      *
*----- Last update  04/08/89/ -----------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*                                                                      *
*     EXTERNALS NONE                                                   *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      REAL*8   time

*-------------------- VARIABLES IN COMMONS -----------------------------

      INCLUDE 'unit.h'
      REAL*8 hartree

*==================== EXECUTABLE STATEMENTS ============================

      pi=4.0D0*DATAN(1.0D0)
      twopi=2.0D0*pi
      avogad=6.0225d23
      boltz=1.38054d-23
      gascon=8.3143d0
      planck=6.6256d-34
      elechg=1.602d-19
      epso=8.854d-12
      boxl=2.0D0
      unitm=1.0D0/(avogad*1000.0D0)
      unitl=1.0D-10
      unitt=1.D-15
      unite=unitm*(unitl/unitt)**2
      unitc=4.0D0*pi*epso*unitl*unite/(elechg*elechg)
      unitp=(unite/unitl**3)/1.0D6
      efact=unite*avogad
      ttstep=time
      hartree = 4.35981 * 1.0d-18
      unitepot=unite/DSQRT(unitc)/hartree
      lbohr = 0.52917706d0
      unitefield=unitepot*lbohr
*==================== END OF EXECUTABLE STATEMENTS =====================

      RETURN
      END
