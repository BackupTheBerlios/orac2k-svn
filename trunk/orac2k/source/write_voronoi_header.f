      SUBROUTINE write_voronoi_header(kvoronoi)

************************************************************************
*   Time-stamp: <97/07/08 23:51:07 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Jul  8 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER kvoronoi

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      WRITE(kvoronoi,100)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

100   FORMAT('********************************************************',
     &       '***********************'/
     &       '*                                                       ',
     &       '                      *'/
     &       '*                            VORONOI OUTPUT FILE        ',
     &       '                      *'/
     &       '*                                                       ',
     &       '                      *'/
     &       '*      All volumes and areas are in A^3 and A^2, respect',
     &       'ively.                *'/
     &       '*                                                       ',
     &       '                      *'/
     &       '*      T        :   Simulation time in fs               ',
     &       '                      *'/
     &       '*      VolRes   :   The volume of a reside.             ',
     &       '                      *'/
     &       '*      VolSlt   :   The total volume of the solute      ',
     &       '                      *'/
     &       '*      VolSlv   :   The total volume of the solvent     ',
     &       '                      *'/
     &       '*      VolMolSlt:   The volume of the molecules of the s',
     &       'olute                 *'/
     &       '*      Area     :   The total surface of each residue of',
     &       ' the solute           *'/
     &       '*                   The first number is the fraction of ',
     &       'surface               *'/
     &       '*                   accessible by the solvent           ',
     &       '                      *'/
     &       '*      Contact  :   The contact surface among selected s',
     &       'olute molecules       *'/
     &       '*      Neigh.   :   The Voronoi coordination number rela',
     &       'tive to the whole     *'/
     &       '*                   solute-solute (SLT-SLT), solvent-sol',
     &       'vent (SLV-SLV) and    *'/
     &       '*                   solute-solvent (SLT-SLV) contacts.  ',
     &       '                      *'/
     &       '*                                                       ',
     &       '                      *'/
     &       '********************************************************',
     &       '***********************'/)
      RETURN
      END
