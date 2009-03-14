      SUBROUTINE fix_molecules(kprint,nprot_fix,prot_fix,mass_pfix,mass
     &     ,protl,nprot)

************************************************************************
*   Time-stamp: <05/02/28 12:33:01 gmarchet>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Jan 30 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER kprint,nprot,nprot_fix,protl(*),prot_fix(*)
      REAL*8  mass(*),mass_pfix

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,count,m,atom,mol,get_pointer_protl

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO j=1,prot_fix(1)
         mol=prot_fix(j+1)
         count=get_pointer_protl(mol-1,protl)
         m=protl(count)
         DO i=1,m
            atom=protl(count+i)
            mass(atom)=mass_pfix
         END DO
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
