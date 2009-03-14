      SUBROUTINE P_load_balance(P_dyn_update,node,nprocs,ncube,nbyte
     &     ,rbyte,nstart_h,nend_h,nlocal_h,nstart_ah,nend_ah,nlocal_ah
     &     ,ntap,ngrp,grppt,worka,cpu_h,ncpu_h,nmapnl,mapnl,mapnl_save)

************************************************************************
*   Time-stamp: <04/12/15 03:07:40 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Fri Feb 19 1999 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      LOGICAL P_dyn_update
      INTEGER ntap,ngrp,node,nprocs,ncube,nbyte,rbyte,nstart_h,nend_h
     &     ,nlocal_h,nstart_ah,nend_ah,nlocal_ah,worka(*),grppt(2,*)
     &     ,ncpu_h

      REAL*8  cpu_h(*)
      INTEGER nmapnl,mapnl(*),mapnl_save(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      IF(P_dyn_update .AND.  nprocs .GT. 1) THEN
         CALL P_split_scalar(node,nprocs,ncube,nbyte,nstart_h
     &        ,nend_h,nlocal_h,nstart_ah,nend_ah,nlocal_ah,ngrp,grppt
     &        ,worka,cpu_h,ncpu_h)
         
*---- Divide mapnl -----------------------------------------------------
         
         CALL icopy(nmapnl,mapnl_save,1,mapnl,1)
         CALL mapnl_divide(node,nstart_h,nend_h,grppt,mapnl)
      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
