MODULE WSC

!!$***********************************************************************
!!$   Time-stamp: <2009-03-07 11:23:27 marchi>                           *
!!$                                                                      *
!!$                                                                      *
!!$                                                                      *
!!$======================================================================*
!!$                                                                      *
!!$              Author:  Massimo Marchi                                 *
!!$              CEA/Centre d'Etudes Saclay, FRANCE                      *
!!$                                                                      *
!!$              - Mon Nov 12 2007 -                                     *
!!$                                                                      *
!!$***********************************************************************

!!$---- This subroutine is part of the program orac2k ----*


!!$======================== DECLARATIONS ================================*

  USE PBC_Mod
  USE INPUT_Mod, ONLY: Read_String, Parser, err_open,err_end,err_unr&
       &,err_fnf,err_args
  IMPLICIT none

  PRIVATE
  PUBLIC Init, Read_it, WSC_, WSC__, simple

  TYPE List
     INTEGER :: n
     INTEGER, ALLOCATABLE :: ind(:)
  END TYPE List
  TYPE(List), ALLOCATABLE, SAVE :: molec_in(:),molec_out(:),molec_tot(:)
  INTEGER, ALLOCATABLE, SAVE :: atoms(:)

  INTEGER, SAVE :: natoms_Tot=0,nats=0,ntot=1,ntap,ntap_in,ntap_out
  LOGICAL, SAVE :: WSC__, simple=.FALSE.
CONTAINS
  SUBROUTINE Init(nprot,protl,atomp)
    IMPLICIT NONE 
    INTEGER nprot,protl(*),atomp(*)
    INTEGER :: count,n,mm,m,counta,countb
    LOGICAL, POINTER :: oks(:)
    TYPE(List), POINTER :: molec_o(:)

    count=0
    ALLOCATE(molec_o(nprot))
    DO n=1,nprot
       mm=protl(count+1)
       molec_o(n) % n=mm
       ALLOCATE(molec_o(n) % ind(mm))
       molec_o(n) % ind(1:mm)=protl(count+1+1:count+1+mm)
       count=count+mm+1
    END DO

    ALLOCATE(oks(nprot))

    oks=.FALSE.
    DO n=1,atoms(1)
       m=atoms(n+1)
       oks(atomp(m))=.TRUE.
    END DO
    m=0
    ntap=0
    DO n=1,nprot
       IF(oks(n)) m=m+1
       ntap=ntap+molec_o(n) % n
    END DO
    ALLOCATE(molec_in(m),molec_out(nprot-m),molec_tot(nprot))
    counta=0
    countb=0
    DO n=1,nprot
       ALLOCATE(molec_tot(n) % ind(molec_o(n) % n))
       molec_tot(n)=molec_o(n)
       IF(.NOT. oks(n)) THEN
          counta=counta+1
          ALLOCATE(molec_out(counta) % ind(molec_o (n) % n))
          molec_out(counta)=molec_o(n)
       ELSE
          countb=countb+1
          ALLOCATE(molec_in(countb) % ind(molec_o(n) % n))
          molec_in(countb) = molec_o(n)
       END IF
    END DO
    ntap_in=SUM(molec_in(:) % n)
    ntap_out=SUM(molec_out(:) % n)
  END SUBROUTINE Init
  SUBROUTINE WSC_(node,xp0,yp0,zp0,coa,aoc)
    IMPLICIT NONE 
    INTEGER :: node
    REAL(8) :: xp0(*),yp0(*),zp0(*),coa(3,3),aoc(3,3)
    INTEGER :: n,n_out,n_in,i,j,mm,l,m,p
    REAL(8) :: xg,yg,zg,rs,tmass,xb,yb,zb,xd,yd,zd,xc,yc&
         &,zc,xcmi,ycmi,zcmi,xcmj,ycmj,zcmj,xcai,ycai,zcai,xcaj,ycaj&
         &,zcaj,xqj,yqj,zqj,rs_min
    TYPE(List), POINTER :: Neigh(:),Neigh_cl(:),Clust_o(:),Clust(:)
    LOGICAL, POINTER :: oks(:),oks_o(:)
    INTEGER, SAVE :: nclust
    REAL(8), POINTER :: xcm(:),ycm(:),zcm(:),xcmc(:),ycmc(:),zcmc(:)
    
    n_in=SIZE(molec_in)
    n_out=SIZE(molec_out)
    ALLOCATE(oks(n_in),oks_o(n_in))

    DO n=1,SIZE(molec_tot)
       mm=molec_tot(n) % n
       xcmi=0.0D0
       ycmi=0.0D0
       zcmi=0.0D0
       DO m=1,mm
          i=molec_tot(n) % ind(m)
          xcmi=xcmi+xp0(i)
          ycmi=ycmi+yp0(i)
          zcmi=zcmi+zp0(i)
       END DO
       xcmi=xcmi/DBLE(molec_tot(n) % n)
       ycmi=ycmi/DBLE(molec_tot(n) % n)
       zcmi=zcmi/DBLE(molec_tot(n) % n)
       xcai=aoc(1,1)*xcmi+aoc(1,2)*ycmi+aoc(1,3)*zcmi
       ycai=aoc(2,1)*xcmi+aoc(2,2)*ycmi+aoc(2,3)*zcmi
       zcai=aoc(3,1)*xcmi+aoc(3,2)*ycmi+aoc(3,3)*zcmi
       xd=-2.0D0*PBC(xcai)
       yd=-2.0D0*PBC(ycai)
       zd=-2.0D0*PBC(zcai)
       xg=coa(1,1)*xd+coa(1,2)*yd+coa(1,3)*zd
       yg=coa(2,1)*xd+coa(2,2)*yd+coa(2,3)*zd
       zg=coa(3,1)*xd+coa(3,2)*yd+coa(3,3)*zd
       DO m=1,mm
          i=molec_tot(n) % ind(m)
          xp0(i)=xp0(i)+xg
          yp0(i)=yp0(i)+yg
          zp0(i)=zp0(i)+zg
       END DO
    END DO

!!$
!!$===== Get C.o.m of first molecule
!!$

    ALLOCATE(xcm(n_in),ycm(n_in),zcm(n_in))
    DO mm=1,n_in
       xcmi=0.0D0
       ycmi=0.0D0
       zcmi=0.0D0
       DO m=1,molec_in(mm) % n
          i=molec_in(mm) % ind(m)
          xcmi=xcmi+xp0(i)
          ycmi=ycmi+yp0(i)
          zcmi=zcmi+zp0(i)
       END DO
       xcmi=xcmi/DBLE(molec_in(mm) % n)
       ycmi=ycmi/DBLE(molec_in(mm) % n)
       zcmi=zcmi/DBLE(molec_in(mm) % n)
       xcai=aoc(1,1)*xcmi+aoc(1,2)*ycmi+aoc(1,3)*zcmi
       ycai=aoc(2,1)*xcmi+aoc(2,2)*ycmi+aoc(2,3)*zcmi
       zcai=aoc(3,1)*xcmi+aoc(3,2)*ycmi+aoc(3,3)*zcmi
       xcm(mm)=xcai
       ycm(mm)=ycai
       zcm(mm)=zcai
    END DO
    CALL Get_Neigh(xcm,ycm,zcm)
    ALLOCATE(Clust_o(n_in))
    nclust=0
    oks_o=.FALSE.
    oks=.FALSE.
    DO 
       m=0
       DO mm=1,n_in
          IF(.NOT. oks(mm)) THEN
             m=mm
             EXIT
          END IF
       END DO
       IF(m == 0) EXIT
       oks(m)=.TRUE.
       CALL Put_Together(oks,m)
       nclust=nclust+1

       p=0
       DO mm=1,n_in
          IF(oks_o(mm) .NEQV. oks(mm)) THEN
             p=p+1
          END IF
       END DO
       Clust_o(nclust) % n=p
       ALLOCATE(Clust_o(nclust) % ind(p))
       p=0
       DO mm=1,n_in
          IF(oks_o(mm) .NEQV. oks(mm)) THEN
             p=p+1
             Clust_o(nclust) % ind(p) = mm
          END IF
       END DO
       oks_o=oks
    END DO
    IF(nclust /= 1 .AND. node == 0) &
         & WRITE(*,'(''==================> '',i6,2x,&
         &'' solute aggregates found '')') nclust

    ALLOCATE(Clust(nclust))
    DO m=1,nclust
       mm=Clust_o(m) % n
       ALLOCATE(Clust(m) % ind(mm))
       Clust(m) = Clust_o(m)
    END DO
    ALLOCATE(xcmc(nclust),ycmc(nclust),zcmc(nclust))
    DO mm=1,nclust
       xcmi=0.0D0
       ycmi=0.0D0
       zcmi=0.0D0
       DO m=1,Clust(mm) % n
          i=Clust(mm) % ind(m)
          xd=coa(1,1)*xcm(i)+coa(1,2)*ycm(i)+coa(1,3)*zcm(i)
          yd=coa(2,1)*xcm(i)+coa(2,2)*ycm(i)+coa(2,3)*zcm(i)
          zd=coa(3,1)*xcm(i)+coa(3,2)*ycm(i)+coa(3,3)*zcm(i)
          xcmi=xcmi+xd
          ycmi=ycmi+yd
          zcmi=zcmi+zd
       END DO
       xcmi=xcmi/DBLE(Clust(mm) % n)
       ycmi=ycmi/DBLE(Clust(mm) % n)
       zcmi=zcmi/DBLE(Clust(mm) % n)
       xcai=aoc(1,1)*xcmi+aoc(1,2)*ycmi+aoc(1,3)*zcmi
       ycai=aoc(2,1)*xcmi+aoc(2,2)*ycmi+aoc(2,3)*zcmi
       zcai=aoc(3,1)*xcmi+aoc(3,2)*ycmi+aoc(3,3)*zcmi
       xcmc(mm)=xcai
       ycmc(mm)=ycai
       zcmc(mm)=zcai
    END DO
    CALL Get_Neigh_Clust(xcmc,ycmc,zcmc)
    DEALLOCATE(oks)
    ALLOCATE(oks(nclust))
    oks=.FALSE.
    p=0
    DO 
       m=0
       DO mm=1,nclust
          IF(.NOT. oks(mm)) THEN
             m=mm
             EXIT
          END IF
       END DO
       IF(m == 0) EXIT
       oks(m)=.TRUE.
       CALL Put_Cluster_Together(oks,m)
       p=p+1
    END DO    

!!$
!!$--- Put the solute in the center of the WSC
!!$

    xcmi=0.0D0
    ycmi=0.0D0
    zcmi=0.0D0
    DO mm=1,n_in
       DO m=1,molec_in(mm) % n
          i=molec_in(mm) % ind(m)
          xcmi=xcmi+xp0(i)
          ycmi=ycmi+yp0(i)
          zcmi=zcmi+zp0(i)
       END DO
    END DO

    xcmi=xcmi/DBLE(ntap_in)
    ycmi=ycmi/DBLE(ntap_in)
    zcmi=zcmi/DBLE(ntap_in)
    xcai=aoc(1,1)*xcmi+aoc(1,2)*ycmi+aoc(1,3)*zcmi
    ycai=aoc(2,1)*xcmi+aoc(2,2)*ycmi+aoc(2,3)*zcmi
    zcai=aoc(3,1)*xcmi+aoc(3,2)*ycmi+aoc(3,3)*zcmi

    xp0(1:ntap)=xp0(1:ntap)-xcmi
    yp0(1:ntap)=yp0(1:ntap)-ycmi
    zp0(1:ntap)=zp0(1:ntap)-zcmi

    xcai=0.0D0
    ycai=0.0D0
    zcai=0.0D0

!!$
!!$--- Displace the remaining molecules
!!$

    DO mm=1,n_out
       xcmj=0.0D0
       ycmj=0.0D0
       zcmj=0.0D0
       DO m=1,molec_out(mm) % n
          j=molec_out(mm) % ind(m)
          xcmj=xcmj+xp0(j)
          ycmj=ycmj+yp0(j)
          zcmj=zcmj+zp0(j)
       END DO
       xcmj=xcmj/DBLE(molec_out(mm) % n)
       ycmj=ycmj/DBLE(molec_out(mm) % n)
       zcmj=zcmj/DBLE(molec_out(mm) % n)

       xcaj=aoc(1,1)*xcmj+aoc(1,2)*ycmj+aoc(1,3)*zcmj
       ycaj=aoc(2,1)*xcmj+aoc(2,2)*ycmj+aoc(2,3)*zcmj
       zcaj=aoc(3,1)*xcmj+aoc(3,2)*ycmj+aoc(3,3)*zcmj
       rs_min=1.0D10
       DO l=-1,1
          DO m=-1,1
             DO n=-1,1
                xg=xcai-xcaj+DBLE(2*l)
                yg=ycai-ycaj+DBLE(2*m)
                zg=zcai-zcaj+DBLE(2*n)
                xc=coa(1,1)*xg+coa(1,2)*yg+coa(1,3)*zg
                yc=            coa(2,2)*yg+coa(2,3)*zg
                zc=                        coa(3,3)*zg
                rs=DSQRT(xc**2+yc**2+zc**2)
                IF(rs .LT. rs_min) THEN
                   rs_min=rs
                   xd=-DBLE(2*l)
                   yd=-DBLE(2*m)
                   zd=-DBLE(2*n)
                END IF
             END DO
          END DO
       END DO
       xqj=xd
       yqj=yd
       zqj=zd
       xg=coa(1,1)*xqj+coa(1,2)*yqj+coa(1,3)*zqj
       yg=coa(2,1)*xqj+coa(2,2)*yqj+coa(2,3)*zqj
       zg=coa(3,1)*xqj+coa(3,2)*yqj+coa(3,3)*zqj
       DO m=1,molec_out(mm) % n
          j=molec_out(mm) % ind(m)
          xp0(j)=xp0(j)+xg
          yp0(j)=yp0(j)+yg
          zp0(j)=zp0(j)+zg
       END DO
    END DO
  CONTAINS
    SUBROUTINE Get_Neigh(xcm,ycm,zcm)
      REAL(8) :: xcm(:),ycm(:),zcm(:)
      INTEGER :: n,m,nn,o,j
      REAL(8), POINTER :: rs(:)
      INTEGER, POINTER :: inda(:),indb(:)
      REAL(8) :: xcmi,ycmi,zcmi,xcmj,ycmj,zcmj,xd,yd,zd,xc,yc,zc,rscut2,rss
      rscut2=10.0D0**2

      nn=SIZE(xcm)
      ALLOCATE(Neigh(nn),rs(nn-1),inda(nn-1),indb(nn-1))

      DO n=1,nn
         xcmi=xcm(n)
         ycmi=ycm(n)
         zcmi=zcm(n)
         o=0
         DO m=1,nn
            IF(m == n) CYCLE
            xcmj=xcm(m)
            ycmj=ycm(m)
            zcmj=zcm(m)
            xd=xcmi-xcmj
            yd=ycmi-ycmj
            zd=zcmi-zcmj
            xd=xd-2.0D0*PBC(xd)
            yd=yd-2.0D0*PBC(yd)
            zd=zd-2.0D0*PBC(zd)
            xc=coa(1,1)*xd+coa(1,2)*yd+coa(1,3)*zd
            yc=coa(2,1)*xd+coa(2,2)*yd+coa(2,3)*zd
            zc=coa(3,1)*xd+coa(3,2)*yd+coa(3,3)*zd
            rss=xc*xc+yc*yc+zc*zc
            IF(rss <= rscut2) THEN
               o=o+1
               rs(o)=rss
               inda(o)=m
            END IF
         END DO
         IF(o > 0) THEN
            Neigh(n) % n = o
            ALLOCATE(Neigh(n) % ind(o))
            IF( o /= 1) THEN
               CALL indexx(o,rs,indb)
               Neigh(n) % ind(:)=inda(indb(1:o))               
            ELSE
               Neigh(n) % ind(1)=inda(1)
            END IF
         END IF
      END DO
    END SUBROUTINE Get_Neigh
    RECURSIVE SUBROUTINE Put_Together(oks,m)
      LOGICAL :: oks(:)
      INTEGER :: m
      INTEGER :: n,j,jj,l,o,p
      REAL(8) :: xcmi,ycmi,zcmi,xcmj,ycmj,zcmj,xcai,ycai,zcai,xcaj,ycaj&
           &,zcaj,xqj,yqj,zqj,rs_min

      xcai=xcm(m)
      ycai=ycm(m)
      zcai=zcm(m)

      DO n=1,Neigh(m) % n
         j=Neigh(m) % ind(n)

         IF(.NOT. oks(j)) THEN
            oks(j)=.TRUE.
            xcaj=xcm(j)
            ycaj=ycm(j)
            zcaj=zcm(j)
            rs_min=1.0D10

            DO l=-1,1
               DO o=-1,1
                  DO p=-1,1
                     xg=xcai-xcaj+DBLE(2*l)
                     yg=ycai-ycaj+DBLE(2*o)
                     zg=zcai-zcaj+DBLE(2*p)
                     xc=coa(1,1)*xg+coa(1,2)*yg+coa(1,3)*zg
                     yc=            coa(2,2)*yg+coa(2,3)*zg
                     zc=                        coa(3,3)*zg
                     rs=DSQRT(xc**2+yc**2+zc**2)
                     IF(rs .LT. rs_min) THEN
                        rs_min=rs
                        xd=-DBLE(2*l)
                        yd=-DBLE(2*o)
                        zd=-DBLE(2*p)
                     END IF
                  END DO
               END DO
            END DO
            xqj=xd
            yqj=yd
            zqj=zd
            xg=coa(1,1)*xqj+coa(1,2)*yqj+coa(1,3)*zqj
            yg=coa(2,1)*xqj+coa(2,2)*yqj+coa(2,3)*zqj
            zg=coa(3,1)*xqj+coa(3,2)*yqj+coa(3,3)*zqj

            DO p=1,molec_in(j) % n
               jj=molec_in(j) % ind(p)
               xp0(jj)=xp0(jj)+xg
               yp0(jj)=yp0(jj)+yg
               zp0(jj)=zp0(jj)+zg
            END DO

            xcm(j)=xcm(j)+xqj
            ycm(j)=ycm(j)+yqj
            zcm(j)=zcm(j)+zqj
            CALL Put_Together(oks,j)
         END IF
      END DO
      RETURN
    END SUBROUTINE Put_Together
    SUBROUTINE Get_Neigh_Clust(xcmc,ycmc,zcmc)
      REAL(8) :: xcmc(:),ycmc(:),zcmc(:)
      INTEGER :: n,m,nn,o,j
      REAL(8), POINTER :: rs(:)
      INTEGER, POINTER :: inda(:),indb(:)
      REAL(8) :: xcmi,ycmi,zcmi,xcmj,ycmj,zcmj,xd,yd,zd,xc,yc,zc,rscut2,rss
      rscut2=40.0D0**2

      nn=SIZE(xcmc)
      ALLOCATE(Neigh_cl(nn),rs(nn-1),inda(nn-1),indb(nn-1))
      DO n=1,nn
         xcmi=xcmc(n)
         ycmi=ycmc(n)
         zcmi=zcmc(n)
         o=0
         DO m=1,nn
            IF(m == n) CYCLE
            xcmj=xcmc(m)
            ycmj=ycmc(m)
            zcmj=zcmc(m)
            xd=xcmi-xcmj
            yd=ycmi-ycmj
            zd=zcmi-zcmj
            xd=xd-2.0D0*PBC(xd)
            yd=yd-2.0D0*PBC(yd)
            zd=zd-2.0D0*PBC(zd)
            xc=coa(1,1)*xd+coa(1,2)*yd+coa(1,3)*zd
            yc=coa(2,1)*xd+coa(2,2)*yd+coa(2,3)*zd
            zc=coa(3,1)*xd+coa(3,2)*yd+coa(3,3)*zd
            rss=xc*xc+yc*yc+zc*zc
            IF(rss <= rscut2) THEN
               o=o+1
               rs(o)=rss
               inda(o)=m
            END IF
         END DO
         IF(o > 0) THEN
            Neigh_cl(n) % n = o
            ALLOCATE(Neigh_cl(n) % ind(o))
            IF( o /= 1) THEN
               CALL indexx(o,rs,indb)
               Neigh_cl(n) % ind(:)=inda(indb(1:o))               
            ELSE
               Neigh_cl(n) % ind(1)=inda(1)
            END IF
         END IF
      END DO
    END SUBROUTINE Get_Neigh_Clust
    RECURSIVE SUBROUTINE Put_Cluster_Together(oks,m)
      LOGICAL :: oks(:)
      INTEGER :: m
      INTEGER :: n,j,jj,l,o,p,pp
      REAL(8) :: xcmi,ycmi,zcmi,xcmj,ycmj,zcmj,xcai,ycai,zcai,xcaj,ycaj&
           &,zcaj,xqj,yqj,zqj,rs_min

      xcai=xcmc(m)
      ycai=ycmc(m)
      zcai=zcmc(m)

      DO n=1,Neigh_Cl(m) % n
         j=Neigh_Cl(m) % ind(n)

         IF(.NOT. oks(j)) THEN
            oks(j)=.TRUE.
            xcaj=xcmc(j)
            ycaj=ycmc(j)
            zcaj=zcmc(j)
            rs_min=1.0D10

            DO l=-1,1
               DO o=-1,1
                  DO p=-1,1
                     xg=xcai-xcaj+DBLE(2*l)
                     yg=ycai-ycaj+DBLE(2*o)
                     zg=zcai-zcaj+DBLE(2*p)
                     xc=coa(1,1)*xg+coa(1,2)*yg+coa(1,3)*zg
                     yc=            coa(2,2)*yg+coa(2,3)*zg
                     zc=                        coa(3,3)*zg
                     rs=DSQRT(xc**2+yc**2+zc**2)
                     IF(rs .LT. rs_min) THEN
                        rs_min=rs
                        xd=-DBLE(2*l)
                        yd=-DBLE(2*o)
                        zd=-DBLE(2*p)
                     END IF
                  END DO
               END DO
            END DO
            xqj=xd
            yqj=yd
            zqj=zd
            xg=coa(1,1)*xqj+coa(1,2)*yqj+coa(1,3)*zqj
            yg=coa(2,1)*xqj+coa(2,2)*yqj+coa(2,3)*zqj
            zg=coa(3,1)*xqj+coa(3,2)*yqj+coa(3,3)*zqj

            DO p=1,Clust(j) % n
               pp=Clust(j) % ind(p)
               xcm(pp)=xcm(pp)+xqj
               ycm(pp)=ycm(pp)+yqj
               zcm(pp)=zcm(pp)+zqj
               DO o=1,Molec_in(pp) % n
                  jj=Molec_in(pp) % ind(o)
                  xp0(jj)=xp0(jj)+xg
                  yp0(jj)=yp0(jj)+yg
                  zp0(jj)=zp0(jj)+zg
               END DO
            END DO
            xcmc(j)=xcmc(j)+xqj
            ycmc(j)=ycmc(j)+yqj
            zcmc(j)=zcmc(j)+zqj
            CALL Put_Cluster_Together(oks,j)
         END IF
      END DO
      RETURN
    END SUBROUTINE Put_Cluster_Together
  END SUBROUTINE WSC_
    
  INCLUDE 'WSC__Read.f90'

!!$----------------- END OF EXECUTABLE STATEMENTS -----------------------*

END MODULE WSC
