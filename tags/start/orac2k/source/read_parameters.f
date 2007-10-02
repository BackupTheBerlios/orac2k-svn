      SUBROUTINE read_parameters(iret,errmsg)

************************************************************************
*   Time-stamp: <1999-11-02 18:28:32 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Nov 19 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER iret
      CHARACTER*80 errmsg

*----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'parst.h'
      INCLUDE 'parameters.h'

*----------------------- EXECUTABLE STATEMENTS ------------------------*

!=======================================================================
!------------ Read in the macromolecular topology and the --------------
!------------ model potential parameters -------------------------------
!=======================================================================

      xnbond=0
      xntor=0
      xnitor=0
      CALL crdrtf(alpha,natype,qchge,natop,jngrp,jgrppt,jbnd,nmbo,jbend
     &     ,nbend,jtor,ntor,jitor,nitor,NOATO,cnat,alphb,BACKP,nback
     &     ,jacc,jdon,nacc,ndon,HYDP,albon,pbon1,pbon2,lpbon,albnd,pbnd1
     &     ,pbnd2,pbnd3,pbnd4,lpbnd,altor,ptor1,ntor2,potj,TORSP,lptor
     &     ,alito,pito1,pito2,pito3,lpito,massu,alhyd,xbond,xnbond,rbond
     &     ,xtor,xntor,rtor,xitor,xnitor,ritor,stpr,stsc,nstpr,nstsc
     &     ,nstrg,NORSRE,NOATRE,BONDP,BENDP,TORSP,ITORP,iret,errmsg)

!----------------- If iret .EQ. 1 STOP !! ------------------------------
!                                 ====

      IF(lpbon .GT. bondp) THEN
          WRITE(errmsg,'(''In JOIN: Dimensions of bond parameters'',
     x   '' exceeded. '',i6,''>'',i6,'' Abort.'')') lpbon,bondp
          iret=1
      END IF

      IF(lpbnd .GT. bendp) THEN
          errmsg='In JOIN: Dimensions of angle parameters exceeded.'//
     x            ' Abort.'
          iret=1
      END IF
      IF(lptor .GT. torsp) THEN
          errmsg='In JOIN: Dimensions of torsional parameters '//
     x            'exceeded. Abort.'
          iret=1
      END IF
      IF(lpito .GT. itorp) THEN
          WRITE(errmsg,'(''In JOIN: Dimensions of i-torsional'',
     x    '' parameters exceeded. LPITO ='',i6)') lpito
          iret=1
      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END

