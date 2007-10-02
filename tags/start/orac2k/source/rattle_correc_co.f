      SUBROUTINE rattle_correc_co(co,dss,cnst,vco,mass,iret
     &     ,errmsg)

************************************************************************
*   Time-stamp: <05/02/28 11:18:20 gmarchet>                             *
*                                                                      *
*   Constraint the co matrix to support only isotropic changes         *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Feb 13 1997 -                                     *
*                                                                      *
************************************************************************


*======================= DECLARATIONS ==================================

      USE Module_Stress
      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER iret,cnst(2,5)
      REAL*8  co(3,3),vco(3,3),mass(*),dss(5)
      CHARACTER*80 errmsg

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER la,lb,k,degree,iox,iter
      INTEGER start,lc
      REAL*8 tol
      REAL*8  dpp,gg,dcnst(2,5),cc(5)
      DATA tol/1.0D-9/degree/5/

*==================== EXECUTABLE STATEMENTS ============================

      start=1
      IF(FixedAngles_Stress) start=3
*=======================================================================
*---- Compute the vectors associated with each constraint and store ----
*---- them in a temporary array ----------------------------------------
*=======================================================================

      DO k=start,degree
         la=cnst(1,k)
         lb=cnst(2,k)
         IF(FixedAngles_Stress) THEN
            lc=lb
         ELSE
            lc=1
         END IF
         dcnst(1,k)=1.0D0/co(lc,lc)
         dcnst(2,k)=-co(la,lb)/co(lc,lc)**2
      END DO

*=======================================================================
*---- RATTLE loop ------------------------------------------------------
*=======================================================================

      iter=0
      iox=1
      DO WHILE(iox .NE. 0)
         iox=0
         DO k=start,degree
            la=cnst(1,k)
            lb=cnst(2,k)
            IF(DABS(dss(k)) .LT. 1.0D-6) THEN
               vco(la,lb)=0.0D0
            ELSE
               IF(DABS(dpp) .GT. tol) iox=1
               IF(FixedAngles_Stress) THEN
                  lc=lb
               ELSE
                  lc=1
               END IF
               dpp=vco(la,lb)-dss(k)*vco(lc,lc)
               cc(k)=dpp
               gg=dcnst(1,k)/mass(la)-dcnst(2,k)/mass(lc)
               gg=dpp/gg
               vco(la,lb)=vco(la,lb)-gg*dcnst(1,k)/mass(la)
               vco(lc,lc)=vco(lc,lc)-gg*dcnst(2,k)/mass(lc)
            END IF
         END DO
         IF(iox.ne.0) THEN
            iter=iter+1
            IF(iter.GT.5000)THEN
               iret=1
               errmsg=
     & ' While RATTLEing CO matrix: The iteration procedure did not'
     x              //' converge.'
               RETURN
            END IF
         END IF
      END DO

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
