      SUBROUTINE comp_total_energy(mapnl,mapdn,nmapdn
     &     ,tag_bndg,fudgec,x,xp0,yp0,zp0,utotal)

************************************************************************
*   Time-stamp: <2005-03-05 22:09:44 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sat Feb 14 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER mapnl(*),mapdn(2,*),nmapdn(*),tag_bndg(*)
      REAL*8  fudgec

      REAL*8  xp0(*),yp0(*),zp0(*),x

*----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'parst.h'
      INCLUDE 'cpropar.h'
      INCLUDE 'unit.h'

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8  puconf,pucoul,puhyd,pubnd,ubend,uptors,uitors,uconf,ucoul
     &     ,ureal,urecp,urcs,urcp,urcsp,urcsp_m,eer,ucns,ucos,ucnsp
     &     ,ucosp,ucnp,ucop,conf_bnd_slt,coul_bnd_slt,conf_bnd_slv
     &     ,coul_bnd_slv,self_slt,self_slv,uslvtor,uslvitor,fscnstr_slt
     &     ,fscnstr_slv,conf_bnd_slt_n1,coul_bnd_slt_n1,conf_bnd_slv_n1
     &     ,coul_bnd_slv_n1,virs,virsp,virp,gsin14,fsbend,fsbond,gsbend
     &     ,gsbond,unb14,cnb14,fsin14,ungrp,cngrp,uumb,gr,ubond,uslvbon
     &     ,uslvben,unbond,cnbond,purecp,upconf,upcoul,press
     &     ,stressd(3,3),stressr(3,3),st(3,3),pv,ustot,uptot,upstot
     &     ,utotal,fpx,fpy,fpz
      INTEGER iret,abmd_dir
      CHARACTER*80 errmsg

      CHARACTER*1 rshell,rshk
      DATA   rshell/'h'/rshk/'h'/

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      CALL get_total_energy(.FALSE.,mapnl,mapdn,nmapdn,tag_bndg,abmd_dir
     &     ,fudgec,xp0,yp0,zp0,fpx,fpy,fpz,stressd,stressr,utotal,ucns
     &     ,ucos,urcs,coul_bnd_slv,conf_bnd_slv_n1,coul_bnd_slv_n1
     &     ,self_slv,fsbend,fsbond,fsin14,unb14,cnb14,uslvbon,uslvben
     &     ,uslvtor,uslvitor,uumb,uptors,uitors,ubond,ubend,ucnp,ucop
     &     ,urcp,conf_bnd_slt_n1,coul_bnd_slt,coul_bnd_slt_n1,self_slt
     &     ,ucnsp,ucosp,urcsp,eer)

      pv=0.0D0
      IF(pressure .OR. cpress .OR. lberendsen) THEN
         pv=(pext*volume*efact)/1000.0D0
      END IF

      ustot=0.d0
      ustot=(ucns+ucos+urcs+uslvbon+uslvben+uslvtor+uslvitor
     &     +conf_bnd_slv_n1+coul_bnd_slv+coul_bnd_slv_n1+self_slv)*efact
     &     /1000.0D0
      uptot=0.d0
      uptot=(urcp+ucop+ucnp+ubond+ubend+uptors+uitors
     &     +conf_bnd_slt_n1+coul_bnd_slt+coul_bnd_slt_n1+self_slt)*efact
     &     /1000.d0

      upstot=0.d0
      upstot = (urcsp+ucosp
     &     +ucnsp)*efact/1000.d0

      utotal= ustot+uptot+upstot+eer*efact/1000.d0
      IF(cpress) THEN
         pv=(pext*volume*efact)/1000.0D0
         utotal=utotal+pv
      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
