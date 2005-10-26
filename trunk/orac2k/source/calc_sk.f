      SUBROUTINE calc_sk(ss_index,co,oc,xp0,yp0,zp0,whe,ntap,rkcut,delsk
     &     ,sofk,nsofk,psofk)

************************************************************************
*                                                                      *
*                                                                      *
*     The contributions from the sum over the reciprocal lattice       *
*     to the solute-solute electrostatic energy and forces             *
*     are computed by FURIPP. The energy is evaluated by               *
*     taking the square of the geometric structure factor              *
*     multiplied by the sites charges. Addition and subtraction        *
*     formulas for cosines and sines are used.                         *
*                                                                      *
*                                                                      *
*     OC      :  Transform the crordinates to simulation box      (I)  *
*                frame.                                                *
*                >> real*8 OC(3,3) <<                                  *
*     XP0     :  Coordinates for the solvent molecules, packed    (I)  *
*     YP0        by molecule.                                          *
*     ZP0        >> real*8  XP0(NTAP), ... <<                          *
*     CHARGE  :  List of atomic charges on the solute molecules   (I)  *
*                as seen by the solvent, packed by molecule.           *
*                >> real*8  CHARGE(NTAP) <<                            *
*     NTAP    :  Total number of sites on the molecules of the    (I)  *
*                solute.                                          (I)  *
*     ALPHAL  :  Ewald sum exponential parameter.                 (I)  *
*     RKCUT   :  Ewald sum cutoff over the k-vectors.             (I)  *
*     V       :  Volume of the simulation box.                    (I)  *
*     URECP   :  Contribution from the reciprocal lattice sum     (O)  *
*                to the electrostatic energy.                          *
*     VIRMM   :  Contribution from the reciprocal lattice sum     (O)  *
*                to the virial.                                        *
*     STRESSTC(3,3) : Coulomb contribution to the stress tensor        *
*     EMVIR         :  Virial                                          *
*                                                                      *
*                                                                      *
*=================== Definition of the sum boundaries =================*
*                                                                      *
*     The sum over the reciprocal space is carried out in the 3        *
*     dimensions from -LMAX, 0, LMAX  -MMAX, 0, MMAX  -NMAX, 0, NMAX   *
*                                                                      *
*     However only half of the reciprocal lattice vectors are used     *
*     The constant lofset, mofset, nofset, lmaxt, mmaxt and nmaxt      *
*     are defined from LMAX, MMAX and NMAX.                            *
*                                                                      *
*     lofset=1                                                         *
*     mofset = MMAX + 1                                                *
*     nofset = NMAX + 1                                                *
*     lmaxt  = LMAX + 1                                                *
*     mmaxt  = 2*MMAX + 1                                              *
*     nmaxt  = 2*NMAX + 1                                              *
*                                                                      *
*========================== WORK ARRAYS ===============================*
*                                                                      *
*     CLM     :  >> real*8  CLM(NA) <<                                 *
*     SLM     :  >> real*8  SLM(NA) <<                                 *
*     EXPIKR  :  >> real*8  EXPIKR(NA) <<                              *
*     EXPIKI  :  >> real*8  EXPIKI(NA) <<                              *
*     ELR     :  >> real*8  ELR(NA,LMAX+1) <<                          *
*     ELI     :  >> real*8  ELI(NA,LMAX+1) <<                          *
*     EMR     :  >> real*8  EMR(NA,MMAX+1) <<                          *
*     EMI     :  >> real*8  EMI(NA,MMAX+1) <<                          *
*     ENR     :  >> real*8  ENR(NA,NMAX+1) <<                          *
*     ENI     :  >> real*8  ENI(NA,NMAX+1) <<                          *
*     NA .GE. NTAP                                                     *
*     If this is not satisfied the routine aborts the run.             *
*                                                                      *
*---- Last update 01/07/93 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi CECAM, Orsay France 1993               *
*                                                                      *
*     EXTERNALS XERROR.                                                *
*                                                                      *
************************************************************************

C======================= DECLARATIONS ==================================

      IMPLICIT none

C----------------------- ARGUMENTS -------------------------------------

      INTEGER ntap,psofk
      INTEGER ss_index(*)
      REAL*8  oc(3,3),co(3,3),rkcut,delsk
      REAL*8  xp0(*),yp0(*),zp0(*),sofk(psofk,2),whe(*),nsofk(psofk,2)

C----------------------- PARAMETERS ------------------------------------

!============================================================================
!---- Parameters of the sum over the reciprocal lattice ---------------------
!============================================================================
      INTEGER lmax,mmax,nmax
      INTEGER lofset,mofset,nofset,lmaxt,mmaxt,nmaxt,mmint,nmint
      PARAMETER (LMAX = 2, MMAX = 2, NMAX = 2, 
     x           lofset = 1, mofset = MMAX + 1, nofset = NMAX + 1,
     x           lmaxt = LMAX + 1, mmaxt = 2*MMAX + 1, 
     x           nmaxt = 2*NMAX + 1)
!============================================================================

C-------------------- DECLARATION OF A SCRATCH COMMON ------------------

      INTEGER nd
      INCLUDE 'parst.h'
      REAL*8 zero,one,two
      PARAMETER(zero=0.0D0,one=1.0D0,two=2.0D0)
      PARAMETER (nd=m1)
      REAL*8  clm(nd),slm(nd),expikr(nd),expiki(nd),elr(nd,LMAX+1)
     &     ,eli(nd,LMAX+1),emr(nd,MMAX+1),emi(nd,MMAX+1),enr(nd,NMAX+1)
     &     ,eni(nd,NMAX+1)
      COMMON / rag1 / clm,slm,expikr,expiki,elr,eli,emr,emi,enr,eni

C-------------------- VARIABLES IN COMMONS -----------------------------

      INCLUDE 'unit.h'

C-------------------- LOCAL VARIABLES ----------------------------------

      CHARACTER*80 errmsg
      INTEGER i,j,l,m,n,mp,np,loff,moff,noff,lmn,lmn1,typei,kbox
      REAL*8  ro(3,3),qpefac,qfrfac,qa,rkscut,tworcl,tworcm,tworcn,qx
     &     ,qy,qz,rl,rm,rn,rolm1,rolm2,rolm3,gl,gm,gn,dsnm,emii,enii
     &     ,rksq,rk,akm,dsnn,sumrp(2),sumip(2)

      REAL*8  sx,sy,sz,aux,delski
      INTEGER k,count,ia,nweight
      REAL*8 weight

C================ EXECUTABLE STATEMENTS ================================

C=======================================================================
C---- Check the dimensions of the work arrays --------------------------
C=======================================================================

      IF(nd.LT.ntap) THEN
         errmsg=' IN CALC_SK : Dimensions of the work arrays are'
     &   / /'insufficient.'
          WRITE(6,*) errmsg
          STOP
      END IF

      IF(DABS(rkcut) .LT. 1.0D-2) RETURN

      nweight=0
      DO i=1,ntap
         nweight=nweight+IDINT(whe(i))
      END DO
      weight=1.0D0/DBLE(nweight)

      delski=1.0D0/delsk
      qpefac=two
      qfrfac=-two
      rkscut=rkcut**2
      tworcl=two*pi/boxl
      tworcm=two*pi/boxl
      tworcn=two*pi/boxl
      mmint=mofset
      nmint=nofset+1
      DO i=1,3
         ro(i,1)=oc(1,i)
         ro(i,2)=oc(2,i)
         ro(i,3)=oc(3,i)
      END DO
      DO i=1,ntap
         qx=tworcl*xp0(i)
         qy=tworcm*yp0(i)
         qz=tworcn*zp0(i)
         elr(i,1)=one
         emr(i,1)=one
         enr(i,1)=one
         eli(i,1)=zero
         emi(i,1)=zero
         eni(i,1)=zero
         elr(i,2)=DCOS(qx)
         emr(i,2)=DCOS(qy)
         enr(i,2)=DCOS(qz)
         eli(i,2)=DSIN(qx)
         emi(i,2)=DSIN(qy)
         eni(i,2)=DSIN(qz)
      END DO
      DO lmn=2,lmax
         lmn1=lmn+1
         DO i=1,ntap
            elr(i,lmn1)=elr(i,lmn)*elr(i,2)-eli(i,lmn)*eli(i,2)
            eli(i,lmn1)=elr(i,lmn)*eli(i,2)+eli(i,lmn)*elr(i,2)
         END DO
      END DO
      DO lmn=2,mmax
         lmn1=lmn+1
         DO i=1,ntap
            emr(i,lmn1)=emr(i,lmn)*emr(i,2)-emi(i,lmn)*emi(i,2)
            emi(i,lmn1)=emr(i,lmn)*emi(i,2)+emi(i,lmn)*emr(i,2)
         END DO
      END DO
      DO lmn=2,nmax
         lmn1=lmn+1
         DO i=1,ntap
            enr(i,lmn1)=enr(i,lmn)*enr(i,2)-eni(i,lmn)*eni(i,2)
            eni(i,lmn1)=enr(i,lmn)*eni(i,2)+eni(i,lmn)*enr(i,2)
         END DO
      END DO
      DO loff=1,lmaxt
         l=loff-lofset
         rl=tworcl*DBLE(l)
         IF(loff.NE.1) THEN
            mmint=1
            nmint=1
         END IF
         DO moff=mmint,mmaxt
            m=moff-mofset
            mp=IABS(m)+1
            rm=tworcm*DBLE(m)
            rolm1=ro(1,1)*rl+ro(1,2)*rm
            rolm2=ro(2,1)*rl+ro(2,2)*rm
            rolm3=ro(3,1)*rl+ro(3,2)*rm
            dsnm=DBLE(SIGN(1,m))
            DO i=1,ntap
               emii=dsnm*emi(i,mp)
               clm(i)=elr(i,loff)*emr(i,mp)-eli(i,loff)*emii
               slm(i)=elr(i,loff)*emii     +eli(i,loff)*emr(i,mp)
            END DO
            DO noff=nmint,nmaxt
               n=noff-nofset
               np=IABS(n)+1
               rn=tworcn*DBLE(n)
               gl=rolm1+ro(1,3)*rn
               gm=rolm2+ro(2,3)*rn
               gn=rolm3+ro(3,3)*rn
               rksq = gl*gl + gm*gm + gn*gn        
               IF(rksq.LE.rkscut) THEN
                  rk=DSQRT(rksq)
                  DO i=1,2
                     sumrp(i)=zero
                     sumip(i)=zero
                  END DO
                  dsnn=dble(sign(1,n))
                  DO i=1,ntap
                     typei=ss_index(i)
                     enii=dsnn*eni(i,np)
                     expikr(i)=enr(i,np)*clm(i)-enii*slm(i)
                     expiki(i)=enii*clm(i)    +enr(i,np)*slm(i)
                     sumrp(typei)=sumrp(typei)+expikr(i)*whe(i)
                     sumip(typei)=sumip(typei)+expiki(i)*whe(i)
                  END DO
                  kbox=MIN0(IDINT(delski*rk+0.5D0),psofk)
                  DO typei=1,2
                     sofk(kbox,typei)=sofk(kbox,typei)+(sumrp(typei)**2
     &                    +sumip(typei)**2)*weight
                     nsofk(kbox,typei)=nsofk(kbox,typei)+1.0D0
                  END DO
               END IF
            END DO
            nmint=1
         END DO
      END DO


C================ END EXECUTABLE STATEMENTS ============================

      RETURN
      END
