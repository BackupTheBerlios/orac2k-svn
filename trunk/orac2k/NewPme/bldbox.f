      SUBROUTINE bldbox(xp0,yp0,zp0,xpg,ypg,zpg,gh,mapnl,mapnl_slv,iret
     &     ,errmsg,node,nprocs,ncube)

************************************************************************
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

      INCLUDE 'parst.h'
      INCLUDE 'cpropar.h'

*----------------------- ARGUMENTS -------------------------------------

      INTEGER iret,ncube
      REAL*8  xp0(*),yp0(*),zp0(*),xpg(*),ypg(*),zpg(*),gh(*)
      CHARACTER*80 errmsg
      INTEGER mapnl(*),mapnl_slv(*)


*-------------------- DECLARATION OF A SCRATCH COMMON ------------------

      INTEGER nb,nc,i
      PARAMETER (nc=192)
c$$$      INTEGER list(3,nb),hlist(3,nb),rlist(3,nb),index(nb),node,nprocs
      INTEGER node,nprocs
      REAL(8), POINTER ::  xpa(:),ypa(:),zpa(:),xpb(:),ypb(:),zpb(:)
      INTEGER, POINTER :: list(:,:),hlist(:,:),rlist(:,:),index(:)

      REAL*8 cg(3,3,nc),tg(3,nc)
c$$$      REAL*8  xpa(nb),ypa(nb),zpa(nb),xpb(nb),ypb(nb),zpb(nb),cg(3,3,nc)
c$$$     &     ,tg(3,nc)
      REAL*8  work(3*m1)
      COMMON /rag1/ xpa,ypa,zpa,xpb,ypb,zpb,cg,tg,work,list,hlist,rlist
     &     ,index

      INTEGER mmb,mmc
      PARAMETER(mmb=40)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER m,ib,i1

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      nb=(nmol+1)*nato_slv+ntap*icl*icm*icn
      ALLOCATE(list(3,nb),hlist(3,nb),rlist(3,nb),index(nb))
      ALLOCATE(xpa(nb),ypa(nb),zpa(nb),xpb(nb),ypb(nb),zpb(nb))

      CALL setup_box(xp0,yp0,zp0,xpg,ypg,zpg,gh,mapnl,mapnl_slv,xpa,ypa
     &     ,zpa,xpb,ypb,zpb,index,cg,tg,mmb,mmc,list,hlist,rlist,nb
     &     ,work,iret,errmsg,node,nprocs,ncube)
      IF(iret .EQ. 1) RETURN

      IF(slv_redef) THEN
         CALL redefine_slv(ss_point,m1+1,ss_index,prsymb,mend,nbun
     &        ,slv_redef_type,mres,grppt,iret,errmsg)
      END IF
      DEALLOCATE(list,hlist,rlist,index)
      DEALLOCATE(xpa,ypa,zpa,xpb,ypb,zpb)
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
