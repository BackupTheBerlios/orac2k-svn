      SUBROUTINE change_coord_ex(nstart_ex,nlocal_ex,nstart_cme
     &     ,nlocal_cme,protl,co,oc,lx,ly,lz,xp0,yp0,zp0,xpa,ypa,zpa
     &     ,xpcma,ypcma,zpcma)

************************************************************************
*   Time-stamp: <99/02/26 17:19:09 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed Feb 11 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nstart_ex,nlocal_ex,nstart_cme,nlocal_cme,protl(*)
      REAL*8  co(3,3),oc(3,3),lx(*),ly(*),lz(*),xp0(*),yp0(*)
     &     ,zp0(*),xpa(*),ypa(*),zpa(*),xpcma(*),ypcma(*),zpcma(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER ptr,get_pointer_protl

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      ptr=get_pointer_protl(nstart_cme-1,protl)
      
*--- After advancing lx,ly,lz recompute xp0,yp0,zp0

      CALL change_origin(-1,nlocal_cme,protl(ptr),xp0,yp0,zp0,lx,ly
     &     ,lz,xpcma(nstart_cme),ypcma(nstart_cme),zpcma(nstart_cme),co)
      
*---        go from real particle coordinates to scaled ----------------
      
      CALL change_frame(co,oc,-1,nlocal_ex,xp0(nstart_ex),yp0(nstart_ex)
     &     ,zp0(nstart_ex),xpa(nstart_ex),ypa(nstart_ex),zpa(nstart_ex))

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
