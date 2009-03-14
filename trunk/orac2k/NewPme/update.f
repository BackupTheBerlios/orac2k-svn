      SUBROUTINE update(co,xp0,yp0,zp0,llacc,lacc,lldon,ldon,hrcut
     &     ,nmapdn,mapdn,o2)

************************************************************************
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
*              - Wed Feb 12 1997 -                                     *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER o2
      INTEGER lldon,ldon(2,*),llacc,lacc(2,*)
      INTEGER mapdn(*),nmapdn(*)
      REAL*8  co(3,3),xp0(*),yp0(*),zp0(*),hrcut

*-------------------- DECLARATION OF A SCRATCH COMMON ------------------

      INTEGER nb
      INCLUDE 'parst.h'
      PARAMETER (nb=m11)
      REAL*8  xmap(nb),ymap(nb),zmap(nb)
      INTEGER mmap(nb)
      LOGICAL maplg(2*nb)
      COMMON / rag1 / xmap,ymap,zmap,mmap,maplg
      INCLUDE 'unit.h'

*-------------------- LOCAL VARIABLES ----------------------------------

      CHARACTER*80 errmsg
      INTEGER i,j,m,n,ia,ira,id,ird,map,ncount,l1,l2,
     x        l3,l4
      REAL*8  radius,x32,y32,z32,xc32,yc32,zc32,rsq32,rsp32
      REAL*8  xpi,ypi,zpi,xsa,ysa,zsa,xca,yca,zca,hrcut2
      INTEGER  npp
      INCLUDE 'pbc.h'

*==================== EXECUTABLE STATEMENTS ============================

      hrcut2=hrcut**2

*=======================================================================
*----- Compute maps for the hydrogen bonds
*=======================================================================

*=======================================================================
*----- Start to loop over all the acceptors of the protein
*=======================================================================

      ncount=0
      DO i=1,llacc
         l2=lacc(2,i)
         map=0
         DO m=1,lldon
            
*=======================================================================
*----- Check if this is a valid hydrogen bond --------------------------
*=======================================================================
            
            l3=ldon(1,m)
            x32=xp0(l3)-xp0(l2)
            y32=yp0(l3)-yp0(l2)
            z32=zp0(l3)-zp0(l2)
            x32=x32-2.0D0*PBC(x32)
            y32=y32-2.0D0*PBC(y32)
            z32=z32-2.0D0*PBC(z32)
            
            xc32=co(1,1)*x32+co(1,2)*y32+co(1,3)*z32
            yc32=            co(2,2)*y32+co(2,3)*z32
            zc32=                        co(3,3)*z32
            
            rsq32=xc32**2+yc32**2+zc32**2
            
*=======================================================================
*----- Now check if interaction is to be included or not ---------------
*=======================================================================
            
            IF(rsq32 .LE. hrcut2) THEN
               map=map+1
               mapdn(map+ncount)=m
            END IF
         END DO
         
*=======================================================================
*----- Nmapdn(i) contains the number of hydrogen bonds -----------------
*----- to be included in the calculation for acceptor i ----------------
*=======================================================================
         
         nmapdn(i)=map
         ncount=ncount+map
         IF(o2.LT.ncount) THEN
            errmsg=' IN UPDATE : Dimensions of the MAPDN are'/ /
     &           ' insufficient.'
            WRITE(kprint,*) o2,ncount
            CALL xerror(errmsg,80,1,2)
         END IF
      END DO
      
      WRITE(kprint,11000) ncount
      
*================= END OF EXECUTABLE STATEMENTS ========================

11000 FORMAT('Hydrogen Bond Neighbor Dimensions *neighbor(',i7
     &     ,')* '/)

      RETURN
      END
