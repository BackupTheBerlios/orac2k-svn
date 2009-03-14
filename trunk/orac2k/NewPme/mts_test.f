      SUBROUTINE mts_test(node,efact,ntmtss,fpx_h,fpx_l,fpx_m,fpx_n
     &     ,ntmtsp,timesec,urcs_h,urcs_l,urcs_m,ucos_h,ucos_l,ucos_m
     &     ,ucns_h,ucns_l,ucns_m,uslvbon,uslvben,uslvtor,uslvitor
     &     ,conf_bnd_slv_n1,coul_bnd_slv,coul_bnd_slv_n1,self_slv,ucek
     &     ,urcp_h,urcp_l,urcp_m,ucop_h,ucop_l,ucop_m,ucnp_h,ucnp_l
     &     ,ucnp_m,ubond,ubend,uitors,uptors,fsbond,fsbend,fsin14,cnb14
     &     ,unb14,cngrp,ungrp,conf_bnd_slt_n1,coul_bnd_slt
     &     ,coul_bnd_slt_n1,self_slt,pucek,urcsp_h,urcsp_l,urcsp_m
     &     ,ucosp_h,ucosp_l,ucosp_m,ucnsp_h,ucnsp_l,ucnsp_m,eer_m,eer_l
     &     ,eer_h,nmol,nato,ntap,nstep,ktest,lfirst,maxstp,cpress,volume
     &     ,pext,ucepr,thermos,uceh,hpot)

****Beta Version: Procacci at CECAM   **********************************
*                                                                      *
*     MTS_TEST prints out at each step istantaneous values of          *
*     Energies and forces of selected atoms for all shell              *
*     contributions at the end of a full r-RESPA propagation step.     *    
*     Output is designed to be easily "awked"                          *
*                                                                      *
*     INPUT: all energies and forces                                   *
*     OUTPUT: none                                                     *
*                                                                      *
*                                                                      *
************************************************************************

c======================= DECLARATIONS ==================================

      IMPLICIT NONE

c----------------------- ARGUMENTS -------------------------------------

      REAL*8 urcs_h,urcs_l,urcs_m,ucos_h,ucos_l,ucos_m,ucns_h,ucns_l
     &     ,ucns_m,uslvbon,uslvben,ucek,urcp_h,urcp_l,urcp_m
     &     ,ucop_h,ucop_l,ucop_m,ucnp_h,ucnp_l,ucnp_m,ubond,ubend,uitors
     &     ,uptors,fsbond,fsbend,fsin14,cnb14,unb14,cngrp,ungrp
     &     ,pucek,urcsp_h,urcsp_l,urcsp_m,ucosp_h,ucosp_l,ucosp_m
     &     ,ucnsp_h,ucnsp_l,ucnsp_m,eer_m,eer_l,eer_h,conf_bnd_slt_n1
     &     ,coul_bnd_slt,coul_bnd_slt_n1,conf_bnd_slv_n1,coul_bnd_slv
     &     ,coul_bnd_slv_n1,self_slt,self_slv,uslvtor,uslvitor,volume
     &     ,pext,ucepr,uceh,hpot

      REAL*8 efact,timesec
      REAL*8 fpx_h(*),fpx_l(*),fpx_m(*),fpx_n(*)
      INTEGER  node,ntmtss(*),ntmtsp(*),nmol,nato,ntap,nstep,ktest
     &     ,maxstp
      LOGICAL lfirst,cpress,thermos

c----------------------- LOCAL VARIABLES -------------------------------

      REAL*8 urcs_h0,urcs_l0,urcs_m0,ucos_h0,ucos_l0,ucos_m0,ucns_h0
     &     ,ucns_l0,ucns_m0,uslvbon0,uslvben0,urcp_h0,urcp_l0
     &     ,urcp_m0,ucop_h0,ucop_l0,ucop_m0,ucnp_h0,ucnp_l0,ucnp_m0
     &     ,ubond0,ubend0,uptors0,uitors0,fsbond0,fsbend0,cnb140
     &     ,unb140,cngrp0,ungrp0,ucosp_h0
     &     ,ucosp_l0,ucosp_m0,ucnsp_h0,ucnsp_l0,ucnsp_m0,tim,pv
      REAL*8 Nurcs_h,Nurcs_l,Nurcs_m,Nucos_h,Nucos_l,Nucos_m,Nucns_h
     &     ,Nucns_l,Nucns_m,Nuslvbon,Nuslvben,Nurcp_h,Nurcp_l
     &     ,Nurcp_m,Nucop_h,Nucop_l,Nucop_m,Nucnp_h,Nucnp_l,Nucnp_m
     &     ,Nubond,Nubend,Nuptors,Nuitors,Nfsbond,Nfsbend,Nfsin14,Ncnb14
     &     ,Nunb14,Ncngrp,Nungrp,Nurcsp_h,Nurcsp_l,Nurcsp_m,Nucosp_h
     &     ,Nucosp_l,Nucosp_m,Nucnsp_h,Nucnsp_l,Nucnsp_m,ustot,uptot
     &     ,ektot,upstot,utot,fnav,pottot
      REAL*8 fpx_h0(5),fpx_l0(5),fpx_m0(5),fpx_n0(5)
      REAL*8 Nfpx_h(5),Nfpx_l(5),Nfpx_m(5),Nfpx_n(5)
      INTEGER  nav,k,kp,j,i
      logical lskip

      tim=timesec*nstep

      lskip=.false.

      do j=1,5
         if(ntmtsp(j).lt.0) lskip=.true.
      end do
      do j=1,3
         if(ntmtss(j).lt.0) lskip=.true.
      end do
   
      if(lskip) goto 1005
c----------------------------------------------------------------------------
c     find first average values 
c----------------------------------------------------------------------------

      IF(lfirst) then
         nav=nav+1
         urcs_h0=urcs_h0 +  urcs_h
         urcs_l0=urcs_l0 +  urcs_l
         urcs_m0=urcs_m0 +  urcs_m
         urcp_h0=urcp_h0 +  urcp_h
         urcp_l0=urcp_l0 +  urcp_l
         urcp_m0=urcp_m0 +  urcp_m
         ucos_h0=ucos_h0 +  ucos_h
         ucos_l0=ucos_l0 +  ucos_l
         ucos_m0=ucos_m0 +  ucos_m
         ucop_h0=ucop_h0 +  ucop_h
         ucop_l0=ucop_l0 +  ucop_l
         ucop_m0=ucop_m0 +  ucop_m
         ucns_h0=ucns_h0 +  ucns_h
         ucns_l0=ucns_l0 +  ucns_l
         ucns_m0=ucns_m0 +  ucns_m
         ucnp_h0=ucnp_h0 +  ucnp_h
         ucnp_l0=ucnp_l0 +  ucnp_l
         ucnp_m0=ucnp_m0 +  ucnp_m
         ucosp_h0=ucosp_h0 + ucosp_h
         ucosp_l0=ucosp_l0 + ucosp_l
         ucosp_m0=ucosp_m0 + ucosp_m
         ucnsp_h0=ucnsp_h0 + ucnsp_h
         ucnsp_l0=ucnsp_l0 + ucnsp_l
         ucnsp_m0=ucnsp_m0 + ucnsp_m
         ubond0=ubond0 +   ubond
         ubend0=ubend0 +   ubend
         uptors0=uptors0 + uptors
         uitors0=uitors0 + uitors
         uslvbon0=uslvbon0 +uslvbon
         uslvben0=uslvben0 +uslvben
         fsbond0=fsbond0 +  fsbond
         fsbend0=fsbend0 +  fsbend
         cnb140=cnb140 +    cnb14
         unb140=unb140 +    unb14
         cngrp0=cngrp0 +    cngrp
         ungrp0=ungrp0 +    ungrp  
         k=0
         DO j=1,5
            IF(ntmtsp(j).ne.0) then
               i=ntmtsp(j) 
               k=k+1
               fpx_n0(k)=fpx_n0(i)+fpx_n(i)
               fpx_m0(k)=fpx_m0(i)+fpx_m(i)
               fpx_l0(k)=fpx_l0(i)+fpx_l(i)
               fpx_h0(k)=fpx_h0(i)+fpx_h(i)
            ENDIF  
         END DO   

c----    averages are done: LFIRST set to .FALSE.

         if(nav.GT.15) THEN
            lfirst=.false.
            fnav=1.d0/float(nav)
         ENDIF

         RETURN
      END IF

c------------------------------------------------------------------------
c     print out total energies
c------------------------------------------------------------------------

1005  ustot=0.d0
      ustot=(ucns_h+ucos_h+ucns_l+ucos_l+ucns_m+ucos_m+urcs_h+urcs_l
     &     +urcs_m+uslvbon+uslvben+uslvtor+uslvitor+conf_bnd_slv_n1
     &     +coul_bnd_slv+coul_bnd_slv_n1+self_slv)*efact/1000.0D0
      uptot=0.d0
      uptot=(urcp_h+urcp_l+urcp_m+ucop_h+ucop_l+ucop_m+ucnp_h+ucnp_l
     &     +ucnp_m+ubond+ubend+uptors+uitors+conf_bnd_slt_n1
     &     +coul_bnd_slt+coul_bnd_slt_n1+self_slt)*efact/1000.d0

      ektot = (ucek+pucek)/1000.d0  
      
      IF(cpress) ektot=ektot+ucepr/1000.0D0
      IF(thermos) ektot=ektot+uceh/1000.0D0

      upstot=0.d0
      upstot = (urcsp_h+urcsp_l+urcsp_m+ucosp_h+ucosp_l+ucosp_m
     &     +ucnsp_h+ucnsp_l+ucnsp_m)*efact/1000.d0

      
         
      utot= ektot+ustot+uptot+upstot+(eer_m+eer_l+eer_h)*efact/1000.d0
      IF(cpress) THEN
         pv=(pext*volume*efact)/1000.0D0
         utot=utot+pv
      END IF
      IF(thermos) utot=utot+hpot*efact/1000.0D0
      
      pottot = utot-ektot
      WRITE(ktest,300) tim,utot,ustot,uptot,upstot,ektot,pottot
300   FORMAT(' TotalEnergy',f12.3,6f15.3)
      IF(node .EQ. 0) WRITE(*,400) tim,utot,ustot,uptot,upstot,ektot
     &     ,pottot
400   FORMAT(' Total energy test',f12.3,3f15.3/15x,3f15.3)

      if(lskip) return

c------------------------------------------------------------------------
c     shift energies
c------------------------------------------------------------------------

c---  shift energies

      Nurcs_h   = urcs_h   -  urcs_h0 *fnav
      Nurcs_l   = urcs_l   -  urcs_l0 *fnav
      Nurcs_m   = urcs_m   -  urcs_m0 *fnav
      Nurcp_h   = urcp_h   -  urcp_h0 *fnav
      Nurcp_l   = urcp_l   -  urcp_l0 *fnav
      Nurcp_m   = urcp_m   -  urcp_m0 *fnav
      Nucos_h   = ucos_h   -  ucos_h0 *fnav
      Nucos_l   = ucos_l   -  ucos_l0 *fnav
      Nucos_m   = ucos_m   -  ucos_m0 *fnav
      Nucop_h   = ucop_h   -  ucop_h0 *fnav
      Nucop_l   = ucop_l   -  ucop_l0 *fnav
      Nucop_m   = ucop_m   -  ucop_m0 *fnav
      Nucns_h   = ucns_h   -  ucns_h0 *fnav
      Nucns_l   = ucns_l   -  ucns_l0 *fnav
      Nucns_m   = ucns_m   -  ucns_m0 *fnav
      Nucnp_h   = ucnp_h   -  ucnp_h0 *fnav
      Nucnp_l   = ucnp_l   -  ucnp_l0 *fnav
      Nucnp_m   = ucnp_m   -  ucnp_m0 *fnav
      Nucosp_h  = ucosp_h  -  ucosp_h0*fnav
      Nucosp_l  = ucosp_l  -  ucosp_l0*fnav
      Nucosp_m  = ucosp_m  -  ucosp_m0*fnav
      Nucnsp_h  = ucnsp_h -   ucnsp_h0*fnav
      Nucnsp_l  = ucnsp_l -   ucnsp_l0*fnav
      Nucnsp_m  = ucnsp_m -   ucnsp_m0*fnav
      Nubond    = ubond   -   ubond0  *fnav
      Nubend    = ubend   -   ubend0  *fnav
      Nuptors   = uptors -    uptors0 *fnav
      Nuitors   = uitors -    uitors0 *fnav 
      Nuslvbon  = uslvbon -   uslvbon0*fnav
      Nuslvben  = uslvben -   uslvben0*fnav
      Nfsbond   = fsbond  -   fsbond0 *fnav
      Nfsbend   = fsbend  -   fsbend0 *fnav
      Ncnb14    = cnb14   -   cnb140  *fnav
      Nunb14    = unb14   -   unb140  *fnav
      Ncngrp    = cngrp   -   cngrp0  *fnav
      Nungrp    = ungrp   -   ungrp0  *fnav   

c---  shift forces

      kp=0
      DO j=1,5
         IF(ntmtsp(j).ne.0) then
            i=ntmtsp(j) 
            kp=kp+1
            Nfpx_n(kp)=fpx_n(i)-fpx_n0(kp)*fnav
            Nfpx_m(kp)=fpx_m(i)-fpx_m0(kp)*fnav
            Nfpx_l(kp)=fpx_l(i)-fpx_l0(kp)*fnav
            Nfpx_h(kp)=fpx_h(i)-fpx_h0(kp)*fnav
         ENDIF  
      END DO   

c------------------------------------------------------------------------
c     print out "shifted" shells contributions for forces and energies)
c------------------------------------------------------------------------

      DO i=1,kp
         WRITE(ktest,200) i,tim,Nfpx_n(i),Nfpx_m(i),Nfpx_l(i)
     &        ,Nfpx_h(i) 
200      FORMAT('P:Force-',i1,f10.3,4e15.5)
      END DO   
      WRITE(ktest,210) tim,Nurcp_h,Nucop_h,Nucnp_h
      WRITE(ktest,220) tim,Nurcp_l,Nucop_l,Nucnp_l
      WRITE(ktest,230) tim,Nurcp_m,Nucop_m,Nucnp_m,Nfsbond,Nfsbend
     &     ,Nfsin14,cngrp,ungrp,unb14,cnb14
      WRITE(ktest,240) tim,ubond,ubend
210   FORMAT(' P:h-',f10.3,4E15.5)
220   FORMAT(' P:l-',f10.3,4E15.5)
230   FORMAT(' P:m-',f10.3,10E15.5)
240   FORMAT(' P:n-',f10.3,4E15.5)
      WRITE(ktest,110) tim,Nurcs_h,Nucos_h,Nucns_h
      WRITE(ktest,120) tim,Nurcs_l,Nucos_l,Nucns_l
      WRITE(ktest,120) tim,Nurcs_m,Nucos_m,Nucns_m
      WRITE(ktest,140) tim,Nuslvben,Nuslvbon
110   FORMAT(' S:h-',f10.3,4E15.5)
120   FORMAT(' S:l-',f10.3,4E15.5)
130   FORMAT(' S:m-',f10.3,4E15.5)
140   FORMAT(' S:n-',f10.3,4E15.5)
      WRITE(ktest,310) tim,Nurcsp_h,Nucosp_h,Nucnsp_h
      WRITE(ktest,320) tim,Nurcsp_l,Nucosp_l,Nucnsp_l
      WRITE(ktest,330) tim,Nurcsp_m,Nucosp_m,Nucnsp_m
310   FORMAT(' PS:h-',f10.3,4E15.5)
320   FORMAT(' PS:l-',f10.3,4E15.5)
330   FORMAT(' PS:m-',f10.3,4E15.5)
      RETURN
      END
