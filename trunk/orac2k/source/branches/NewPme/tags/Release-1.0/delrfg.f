      SUBROUTINE delrfg(dss,coeff,cnstp,cnstpp,bsitp,nbsitp,ma,
     x                  nrigg,iret,errmsg)

************************************************************************
*                                                                      *
*     Delete secondary sites from the list of constranints.            *
*                                                                      *
*                                                                      *
*---- Last update 11/11/92 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi CECAM, Universite Paris-Sud            *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER ma,iret
      INTEGER cnstpp,nrigg,bsitp(ma,*),nbsitp(*),cnstp(2,*)
      REAL*8  dss(*),coeff(*)
      CHARACTER*80 errmsg

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER n,m,l,a,b,ncount
      LOGICAL ok

*-------------------- Work Array ---------------------------------------

      INCLUDE 'parst.h'
      INTEGER mbs
      PARAMETER (mbs=m1+100)
      INTEGER gather
      COMMON /rag1/ gather(mbs)

*==================== EXECUTABLE STATEMENTS ============================

      ncount=0
      DO l=1,cnstpp
          a=cnstp(1,l)
          b=cnstp(2,l)
          ok=.TRUE.
          DO n=1,nrigg
              DO m=1,nbsitp(n)
                  IF( bsitp(m,n) .EQ. a .OR. bsitp(m,n) .EQ. b) THEN
                      ok=.FALSE.
                  END IF
              END DO
          END DO
          IF(ok) THEN
              ncount=ncount+1
              gather(ncount)=l
              IF(ncount .GT. mbs) THEN
                  iret=1
                  errmsg='Error in DELRFG: Physical dimensions of the'
     x// ' work arrays insufficients. Abort. '
                  RETURN
              END IF
          END IF
      END DO
      DO n=1,ncount
          dss(n)=dss(gather(n))
          coeff(n)=coeff(gather(n))
          cnstp(1,n)=cnstp(1,gather(n))
          cnstp(2,n)=cnstp(2,gather(n))
      END DO
      cnstpp=ncount

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
