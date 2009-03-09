MODULE PBC_Mod

!!$***********************************************************************
!!$   Time-stamp: <2009-03-07 11:24:51 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Fri Feb 10 2006 -                                     *
!!$                                                                      *
!!$***********************************************************************
CONTAINS
  REAL(8) FUNCTION PBC(x)
    REAL(8) :: x
    PBC=DNINT(0.5D0*x)
  END FUNCTION PBC
END MODULE PBC_MOD
