      SUBROUTINE gauss(tmp,mass,nato,nmol,vx,vy,vz)

************************************************************************
*                                                                      *
*                                                                      *
*     GAUSS will generate a set of gaussianly distributed              *
*     random velocities with average zero and standard deviation       *
*     DSQRT(tmp).                                                      *
*                                                                      *
*                                                                      *
*     TMP     :  Set temperature.                                 (I)  *
*     OC      :  Transform the coordinates to simulation          (I)  *
*                box frame.                                            *
*                >> real*8 OC(3,3) <<                                  *
*     MASS    :  List of atomic masses.                           (I)  *
*                >> real*8 MASS(NATO*NMOL) <<                          *
*     NATO    :  Number of atoms in the molecule with non         (I)  *
*                zero mass.                                            *
*     NCNST   :  Number of constraints per molecule.              (I)  *
*     NMOL    :  Number of molecules.                             (I)  *
*     VX      :  List of atomic velocities.                       (O)  *
*     VY         >> real*8 VX(NATO*NMOL), ... <<                       *
*     VZ                                                               *
*                                                                      *
*---- Last update 04/26/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNAL RANF.                                                   *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER nato,nmol
      REAL*8  mass(*),vx(*),vy(*),vz(*),tmp

*-------------------- LOCAL VARIABLES ----------------------------------

      EXTERNAL ranf
      INTEGER i,nf,nts
      REAL*8  u1,u2,u3,u4,sig,dummy,t1,t2,t3,ranf,tvel

*-------------------- VARIABLES IN COMMONS -----------------------------

      INCLUDE 'unit.h'
      LOGICAL near0
      DATA dummy/0.0D+0/

*==================== EXECUTABLE STATEMENTS ============================

      nts=nato*nmol
      DO 10 i=1,nts
          tvel=0.0D0
          IF(.NOT. near0(mass(i))) tvel=boltz*tmp/(unite*mass(i))
          sig=DSQRT(tvel)
          u1=ranf(dummy)
          u2=ranf(dummy)
          u3=ranf(dummy)
          u4=ranf(dummy)
          IF(u1 .EQ. 0.0D0) THEN
              u1=1.0D-15
          END IF
          IF(u2 .EQ. 0.0D0) THEN
              u2=1.0D-15
          END IF
          IF(u3 .EQ. 0.0D0) THEN
              u3=1.0D-15
          END IF
          IF(u4 .EQ. 0.0D0) THEN
              u4=1.0D-15
          END IF
          t1=DSQRT(-2.0d0*DLOG(u1))*DCOS(2.0d0*pi*u2)
          t2=DSQRT(-2.0d0*DLOG(u1))*DSIN(2.0d0*pi*u2)
          t3=DSQRT(-2.0d0*DLOG(u3))*DCOS(2.0d0*pi*u4)
          vx(i)=t1*sig
          vy(i)=t2*sig
          vz(i)=t3*sig
10    CONTINUE

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
