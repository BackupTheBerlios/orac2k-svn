!!$***********************************************************************
!!$   Time-stamp: <2005-02-19 18:13:49 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Mon Feb 14 2005 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This module is part of the program ORAC ----*

MODULE Poisson_Boltzmann
  USE Module_Fourier, ONLY: Fourier_Out, Do_Fourier
  USE Class_Gcharges, ONLY: Gcharges_Init=>Init, Get_Phi,&
       & Gauss_Charges, Get_Density, Make_Copy
CONTAINS
  SUBROUTINE Ions_Boltzmann(Boltz)
!!$======================== DECLARATIONS ================================*

    IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*
  
    TYPE(Gauss_Charges) :: Boltz

!!$----------------------- VARIABLES IN COMMON --------------------------*

!!$------------------------- LOCAL VARIABLES ----------------------------*

    TYPE(Gauss_Charges) :: Boltz_new
    
    REAL(8) :: rms
    INTEGER :: i
!!$----------------------- EXECUTABLE STATEMENTS ------------------------*

    WRITE(*,*) Boltz % label
    WRITE(*,*) Boltz % chg(1:10)

    DO i=1,30
       Boltz_new=Get_Phi(Boltz)
       Boltz_new=Get_Density(Boltz_new)
!!$       WRITE(*,*) Boltz_New % chg(1:2)
!!$       WRITE(*,*) Boltz % chg(1:2)
       rms=Get_Chi()
       Boltz=Boltz_new
       WRITE(*,*) i, rms
    END DO
    WRITE(88,'(f12.6,i8)') (Boltz % chg(i),i,i=1,Boltz % natom)
    STOP
  CONTAINS
    FUNCTION Get_Chi() RESULT(out)
      IMPLICIT none
      REAL(8) :: out

      out=SQRT(SUM((Boltz_New % chg - Boltz % chg)**2))
    END FUNCTION Get_Chi
  END SUBROUTINE Ions_Boltzmann
!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE Poisson_Boltzmann
