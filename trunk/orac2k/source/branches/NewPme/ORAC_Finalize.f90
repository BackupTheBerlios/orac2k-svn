SUBROUTINE ORAC_Finalize

!!$***********************************************************************
!!$   Time-stamp: <2006-02-24 12:33:36 marchi2>                             *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Feb 24 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program  ----*


!!$======================== DECLARATIONS ================================*

 IMPLICIT none

!!$----------------------------- ARGUMENTS ------------------------------*



!!$----------------------- VARIABLES IN COMMON --------------------------*

!!$------------------------- LOCAL VARIABLES ----------------------------*

 INTEGER :: ierr

!!$----------------------- EXECUTABLE STATEMENTS ------------------------*



#if defined PARALLEL
      CALL MPI_FINALIZE(ierr)
#else
      STOP
#endif

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END SUBROUTINE ORAC_Finalize
