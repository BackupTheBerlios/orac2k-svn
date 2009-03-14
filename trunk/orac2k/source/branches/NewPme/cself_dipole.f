      SUBROUTINE cself_dipole(nato,alphal_c,alphal_d,charge,dipole
     &     ,volume,self,self_dip)

************************************************************************
*                                                                      *
*                                                                      *
*     Compute the Ewald self term.                                     *
*                                                                      *
*---- Last update 06/12/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNALS  NONE                                                  *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER nato
      REAL*8  charge(*),alphal_c,alphal_d,self,dipole(3,*),volume
     &     ,self_dip

*------------------ VARIABLES IN COMMON --------------------------------

      INCLUDE 'unit.h'

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,n
      REAL*8  fsrtal,tot,cor_dip,TotCharge,VolCor

*==================== EXECUTABLE STATEMENTS ============================


*=======================================================================
*---- Compute the self term --------------------------------------------
*=======================================================================
      
      fsrtal=0.0D0
      cor_dip=0.0D0
      TotCharge=0.0D0
      DO n=1,nato
         fsrtal=fsrtal+charge(n)**2
         cor_dip=cor_dip+dipole(1,n)**2+dipole(2,n)**2
     &        +dipole(3,n)**2
         TotCharge=TotCharge+charge(n)
      END DO
      tot=-(alphal_c/DSQRT(pi))*fsrtal
      self_dip=-(alphal_d/DSQRT(pi))*(2.0D0/3.0D0)*(alphal_d**2)*cor_dip
      

*=======================================================================
*---- Don't forget volume correction if total charge is non zero! ------
*=======================================================================

      VolCor=0.5D0*(TotCharge**2)*pi/alphal_c**2/volume

      self=tot-VolCor

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
