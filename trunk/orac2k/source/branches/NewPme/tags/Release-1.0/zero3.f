      SUBROUTINE zero3(ax,ay,az,nstart,nend)

************************************************************************
*                                                                      *
*                                                                      *
*     This subroutine set the elements of 3 arrays to zero.            *
*                                                                      *
*     AX      :  Arrays whose element must be set to zero.             *
*     AY         >> real*8 AX(N*M), AY(N*M), AZ(N*M) <<                *
*     AZ                                                               *
*                                                                      *
*     N       :  Physical row dimension of A.                          *
*     M       :  Physical column dimension of A.                       *
*                                                                      *
*---- Last update 04/23/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNALS NONE                                                   *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER nstart,nend
      REAL*8 ax(*),ay(*),az(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i

*==================== EXECUTABLE STATEMENTS ============================

      DO i=nstart,nend
         ax(i)=0.0D0
         ay(i)=0.0D0
         az(i)=0.0D0
      END DO

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
