      SUBROUTINE add_solvent_tpg(mapnl,mapnl_slv,iret,errmsg)

************************************************************************
*   Time-stamp: <99/04/13 15:13:51 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Mar 14 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none
      INCLUDE 'parst.h'
      INCLUDE 'cpropar.h'

*----------------------- ARGUMENTS -------------------------------------

      INTEGER iret,mapnl(*),mapnl_slv(*)
      CHARACTER*80 errmsg

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER idummy,i
      REAL*8  sum

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      CALL add_int2(.TRUE.,ntap,nato_slv,lacc,lacc_slv,llacc,llacc_slv,2
     &     ,m1,nmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In ADD_SOLVENT_TPG: Dimension exceeded for lacc. Abort'
         RETURN
      END IF

      CALL add_int2(.TRUE.,ntap,nato_slv,ldon,ldon_slv,lldon,lldon_slv,2
     &     ,m1,nmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In ADD_SOLVENT_TPG: Dimension exceeded for ldon. Abort'
         RETURN
      END IF

      CALL add_int1(.FALSE.,ntap,nato_slv,ntap,nato_slv,nbtype
     &     ,nbtype_slv,m1,nmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg
     &     ='In ADD_SOLVENT_TPG: Dimension exceeded for nbtype. Abort'

         RETURN
      END IF

      CALL add_int1(.TRUE.,nbun,nbun_slv,ntap,nato_slv,nres(1,1)
     &     ,nres_slv(1,1),m1,nmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In ADD_SOLVENT_TPG: Dimension exceeded for nres. Abort'
         RETURN
      END IF

      CALL add_int1(.FALSE.,nbun,nbun_slv,ntap,nato_slv,nres(1,2)
     &     ,nres_slv(1,2),m1,nmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In ADD_SOLVENT_TPG: Dimension exceeded for nres. Abort'
         RETURN
      END IF

      CALL add_int1(.TRUE.,ntap,nato_slv,ntap,nato_slv,mback,mback_slv
     &     ,m1,nmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg
     &        ='In ADD_SOLVENT_TPG: Dimension exceeded for mback. Abort'
         RETURN
      END IF

      CALL add_int1(.FALSE.,idummy,idummy,nbun,nbun_slv,mend,mend_slv
     &     ,nores,nmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In ADD_SOLVENT_TPG: Dimension exceeded for mend. Abort'
         RETURN
      END IF

      CALL add_int2(.TRUE.,ntap,nato_slv,grppt,grppt_slv,ngrp,ngrp_slv
     &     ,2,m11,nmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg
     &        ='In ADD_SOLVENT_TPG: Dimension exceeded for grptt. Abort'
         RETURN
      END IF

      CALL add_int2(.TRUE.,ntap,nato_slv,lbnd,lbnd_slv,lbond,lbond_slv,2
     &     ,m9,nmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In ADD_SOLVENT_TPG: Dimension exceeded for lbnd. Abort'
         RETURN
      END IF

      CALL add_int2(.TRUE.,ntap,nato_slv,lbndg,lbndg_slv,lbend,lbend_slv
     &     ,3,m2,nmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg
     &    ='In ADD_SOLVENT_TPG: Dimension exceeded for lbndg. Abort'
         RETURN
      END IF

      CALL add_int2(.TRUE.,ntap,nato_slv,ltor,ltor_slv,ltors,ltors_slv,4
     &     ,m3,nmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In ADD_SOLVENT_TPG: Dimension exceeded for ltor. Abort'
         RETURN
      END IF

      CALL add_int2(.TRUE.,ntap,nato_slv,litr,litr_slv,litor,litor_slv,4
     &     ,m4,nmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In ADD_SOLVENT_TPG: Dimension exceeded for litr. Abort'
         RETURN
      END IF
 
      CALL add_int2(.TRUE.,ntap,nato_slv,int14,int14_slv,int14p
     &     ,int14p_slv,2,m3,nmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg
     &    ='In ADD_SOLVENT_TPG: Dimension exceeded for int14. Abort'
         RETURN
      END IF

      CALL add_int1(.FALSE.,ntap,nato_slv,int14p,int14p_slv,type14
     &     ,type14_slv,m3,nmol,iret)

      IF(iret .EQ. 1) THEN
         errmsg
     &    ='In ADD_SOLVENT_TPG: Dimension exceeded for type14. Abort'
         RETURN
      END IF

      CALL add_int2(.TRUE.,ntap,nato_slv,int13,int13_slv,int13p
     &     ,int13p_slv,2,m2,nmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg
     &    ='In ADD_SOLVENT_TPG: Dimension exceeded for int13. Abort'
         RETURN
      END IF

*****Stopped here for the moment ******************************

      CALL add_concta(ntap,nato_slv,concta,concta_slv,ntap,nato_slv,m10
     &     ,m1,slvatm,nmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg
     &    ='In ADD_SOLVENT_TPG: Dimension exceeded for concta. Abort'
         RETURN
      END IF

      CALL add_mapnl(ntap,nato_slv,ntap,nato_slv,mapnl,mapnl_slv,m8,nmol
     &     ,iret)
      IF(iret .EQ. 1) THEN
         errmsg
     &    ='In ADD_SOLVENT_TPG: Dimension exceeded for mapnl. Abort'
         RETURN
      END IF

      CALL add_real1(ntap,nato_slv,chrge,chrge_slv,m1,nmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg
     &    ='In ADD_SOLVENT_TPG: Dimension exceeded for chrge. Abort'
         RETURN
      END IF

      CALL add_real1(ntap,nato_slv,mass,mass_slv,m1,nmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In ADD_SOLVENT_TPG: Dimension exceeded for mass. Abort'
         RETURN
      END IF

****  Need to do arrays type potbe(i,2) e concta(i,j) type of arrays and 
****  nrigg, prsymb, mend and character arrays

      CALL add_real2(2,lbond,lbond_slv,potbo,potbo_slv,m9,slv2,nmol,iret
     &     )
      IF(iret .EQ. 1) THEN
         errmsg
     &   ='In ADD_SOLVENT_TPG: Dimension exceeded for potbo. Abort'
         RETURN
      END IF

      CALL add_real2(4,lbend,lbend_slv,potbe,potbe_slv,m2,slv3,nmol,iret
     &     )
      IF(iret .EQ. 1) THEN
         errmsg
     &   ='In ADD_SOLVENT_TPG: Dimension exceeded for potbe. Abort'
         RETURN
      END IF

      CALL add_real2(2,ltors,ltors_slv,potto,potto_slv,m3,slv4,nmol,iret
     &     )
      IF(iret .EQ. 1) THEN
         errmsg
     &    ='In ADD_SOLVENT_TPG: Dimension exceeded for ptorj. Abort'
         RETURN
      END IF
      
      CALL add_real2(3,litor,litor_slv,potit,potit_slv,m4,slv5,nmol,iret
     &     )
      IF(iret .EQ. 1) THEN
         errmsg
     &    ='In ADD_SOLVENT_TPG: Dimension exceeded for ptorj. Abort'
         RETURN
      END IF

      IF(nrigg .NE. 0) THEN
         iret=1
         errmsg
     &        ='In ADD_SOLVENT_TPG: Rigid atomic groups are '
     &        //'incompatible with symmetry. Abort'
         RETURN
      END IF
      
      CALL add_char1(7,ntap,nato_slv,beta,beta_slv,m1,nmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In ADD_SOLVENT_TPG: Dimension exceeded for beta. Abort'
         RETURN
      END IF

      CALL add_char1(7,ntap,nato_slv,betb,betb_slv,m1,nmol,iret)
      IF(iret .EQ. 1) THEN
         errmsg='In ADD_SOLVENT_TPG: Dimension exceeded for betb. Abort'
         RETURN
      END IF

      nbun=nbun+nbun_slv*nmol
      llacc=llacc+llacc_slv*nmol
      lldon=lldon+lldon_slv*nmol
      nbone=nbone+nbone_slv*nmol
      ntap=ntap+nato_slv*nmol
      lbond=lbond+lbond_slv*nmol
      lbend=lbend+lbend_slv*nmol
      ltors=ltors+ltors_slv*nmol
      litor=litor+litor_slv*nmol
      int14p=int14p+int14p_slv*nmol
      int13p=int13p+int13p_slv*nmol
      ngrp=ngrp+ngrp_slv*nmol
      sum=0.0D0
      DO i=1,nato_slv
         sum=sum+mass_slv(i)
      END DO
      wmtp=wmtp+DFLOAT(nmol)*sum

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
