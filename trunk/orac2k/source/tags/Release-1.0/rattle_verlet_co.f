      SUBROUTINE rattle_verlet_co(dt,co1,co,dss,cnst,vco,mass,iret
     &     ,errmsg)

************************************************************************
*   Time-stamp: <05/02/28 11:20:45 gmarchet>                             *
*                                                                      *
*   Constraint the co matrix to support only isotropic changes         *
*                                                                      *
*   The following constraints are satisfied:                           *
*                                                                      *
*        \sigma_{1} = {h_{22} \over h_{11}} - C_1                      *
*        \sigma_{2} = {h_{33} \over h_{11}} - C_2                      *
*        \sigma_{3} = {h_{12} \over h_{11}} - C_3                      *
*        \sigma_{4} = {h_{13} \over h_{11}} - C_4                      *
*        \sigma_{5} = {h_{23} \over h_{11}} - C_5                      *
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
      REAL*8  co1(3,3),co(3,3),vco(3,3),mass(*),dss(5),dt
      CHARACTER*80 errmsg

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER la,lb,k,degree,iox,iter
      INTEGER start,lc
      REAL*8 tol
      REAL*8  dcnt,dti,dp1,dp2,dpp,gg,dps1,dps2,dcnst(2,5)
      DATA tol/1.0D-9/degree/5/

*==================== EXECUTABLE STATEMENTS ============================

      start=1
      IF(FixedAngles_Stress) start=3
      dti=1.0D0/dt

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
*---- SHAKE loop -------------------------------------------------------
*=======================================================================

      iter=0
      iox=1
      DO WHILE(iox .NE. 0)
         iox=0
         DO k=start,degree
            la=cnst(1,k)
            lb=cnst(2,k)
            IF(DABS(dss(k)) .LT. 1.0D-6) THEN
               co1(la,lb)=co(la,lb)
            ELSE
               IF(FixedAngles_Stress) THEN
                  lc=lb
               ELSE
                  lc=1
               END IF
               dp1=1.0D0/co1(lc,lc)
               dp2=-co1(la,lb)/co1(lc,lc)
               dpp=co1(la,lb)/co1(lc,lc)
               dcnt=dpp-dss(k)
               IF(DABS(dcnt) .GT. tol) THEN
                  iox=1
                  dps1=dp1*dcnst(1,k)
                  dps2=dp2*dcnst(2,k)
                  gg=dps1/mass(la)+dps2/mass(lc)
                  gg=-dcnt/gg
                  co1(la,lb)=co1(la,lb)+dcnst(1,k)*gg/mass(la)
                  co1(lc,lc)=co1(lc,lc)+dcnst(2,k)*gg/mass(lc)

                  gg=gg*dti

                  vco(la,lb)=vco(la,lb)+dcnst(1,k)*gg/mass(la)
                  vco(lc,lc)=vco(lc,lc)+dcnst(2,k)*gg/mass(lc)
               END IF
            END IF
         END DO
         IF(iox.ne.0) THEN
            iter=iter+1
            IF(iter.GT.5000)THEN
               iret=1
               errmsg=
     & ' While SHAKEing CO matrix: The iteration procedure did not'
     x              //' converge.'
               RETURN
            END IF
         END IF
      END DO

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
