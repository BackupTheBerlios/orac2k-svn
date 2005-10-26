      SUBROUTINE restim(anprot,annpro,anpoint,res,dres,mres,ntap,nmol,
     x           mend,nbun,grppt,ncorr,nmax,tstar,time,beta,prsymb)

************************************************************************
*                                                                      *
*     Calculates residence time per group at the end of the analisys   *
*     run.                                                             *
*                                                                      *
*                                                                      *
*---- Last update 11/13/92 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi CEA Saclay 1993                        *
*                                                                      *
*     EXTERNAL NONE                                                    *
*                                                                      *
************************************************************************

      IMPLICIT none

      INTEGER nmax,ncorr,dres
      INTEGER ntap,nmol,nbun,annpro,mend(*),res(dres,*),mres(2,*),
     x        anpoint(2,*),grppt(2,*),tstar
      LOGICAL anprot
      CHARACTER*8 prsymb(*)
      CHARACTER*7 beta(*)
      REAL*8  time

      INCLUDE 'resid.h'
      INCLUDE 'unit.h'

      INTEGER dcorr,k1,k2
      PARAMETER (dcorr=150,k1=20,k2=150)
      REAL*8  corr(0:dcorr,dlst),tcorr(0:dcorr),der(0:dcorr)
      INTEGER ncc(0:dcorr,dlst),tlst(k1,dtm),p(dtm,k2),boundr(dlst),
     x        boundw(dlst),boundrr(dlst),boundww(dlst)
      LOGICAL label(dlst)
      COMMON /rag1/ corr,tcorr,boundr,boundw,boundww,boundrr,der,
     x              ncc,tlst,p,label

      CHARACTER*80 errmsg
      INTEGER i,j,n,l,k,mstart,mfin,nstart,nend,nc,pk,nmc,tn,idermin,
     x        nbd,nbd1,m
      REAL*8  dermin
      LOGICAL ok, near0

      IF(ncorr .GT. dcorr) THEN
          errmsg=
     x'In RESTIM: Correlation array underdimensioned. Abort'
          CALL xerror(errmsg,80,1,2)
      END IF

      nbd=0      
      nbd1=0      
      DO n=1,annpro
          nstart=anpoint(1,n)
          nend=anpoint(2,n)
          mstart=mres(1,res(nstart,1))
          mfin=mres(2,res(nend,1))
          DO l=mstart,mfin

!=======================================================================
!----- Zero correlation function ---------------------------------------
!=======================================================================
 
              DO k=0,ncorr
                  corr(k,l)=0.0D0
                  ncc(k,l)=0
              END DO
              label(l)=.FALSE.

!=======================================================================
!----- For each group calculate the interacting water molecules --------
!----- and put it in tlist(i,k) ----------------------------------------
!=======================================================================

              DO k=1,nmax
                  IF(wnlst(l,k) .GT. k1) THEN
                      errmsg=
     x'In RESTIM: Temporary array tlst dimensions exceeded. Abort.'
                      CALL xerror(errmsg,80,1,2)
                  END IF
                  nc=0
                  DO i=1,l-1
                      nc=nc+wnlst(i,k)
                  END DO
                  DO i=1,wnlst(l,k)
                      nc=nc+1
                      tlst(i,k)=wlst(nc,k)
                  END DO
              END DO

!=======================================================================
!----- Calculate how many and which water molecules have been in  ------
!----- the first solv. shell of group l. Put results in p(k,nc)   ------
!----- number of waters is nmc -----------------------------------------
!=======================================================================

              nc=1
              DO k=1,nmax
                  p(k,nc)=0
              END DO
              DO i=1,nmol
                  ok=.FALSE.
                  DO k=1,nmax
                      DO j=1,wnlst(l,k)
                          IF(i .EQ. tlst(j,k)) THEN
                              ok=.TRUE.
                              p(k,nc)=1
                          END IF
                      END DO
                  END DO
                  IF(ok) THEN

!=======================================================================
!----- Calculate water molecules which always reside near a group ------
!=======================================================================

                      pk=0 
                      DO k=1,nmax 
                          pk=pk+p(k,nc)
                      END DO
                      IF(pk .EQ. nmax) THEN
                          nbd=nbd+1
                          IF(nbd .GT. dlst) THEN
                              errmsg=
     x'In RESTIM: Temporary array boundr dimensions exceeded. Abort.'
                              CALL xerror(errmsg,80,1,2)
                          END IF
                          boundr(nbd)=l
                          boundw(nbd)=i
                      END IF
                      IF(pk .NE. nmax .AND. pk .GT. 0.9*nmax) THEN
                          nbd1=nbd1+1
                          IF(nbd1 .GT. dlst) THEN
                              errmsg=
     x'In RESTIM: Temporary array boundr dimensions exceeded. Abort.'
                              CALL xerror(errmsg,80,1,2)
                          END IF
                          boundrr(nbd1)=l
                          boundww(nbd1)=i
                      END IF

!=======================================================================
!----- Prepare for the next water molecule -----------------------------
!=======================================================================

                      nc=nc+1
                      IF(nc .GT. k2) THEN
                          errmsg=
     x'In RESTIM: Temporary array p(i,j) dimensions exceeded. Abort.'
                          CALL xerror(errmsg,80,1,2)
                      END IF
       
                      DO k=1,nmax
                          p(k,nc)=0
                      END DO
                  END IF
              END DO
              nmc=nc-1
              IF(nmc .EQ. 0) GOTO 100

!=======================================================================
!----- Calculate correlation function for group l ----------------------
!=======================================================================
            
              DO i=1,nmax
                  tn=ncorr+i
                  IF(tn .GT. nmax) tn=nmax
                  DO j=i,tn
                      ncc(j-i,l)=ncc(j-i,l)+1
                  END DO
                  DO k=1,nmc
                      pk=0
                      DO j=i,tn
                          pk=pk+p(j,k)
                          IF(pk .NE. 0) THEN
                              IF(pk+tstar .GE. j-i) THEN
                                  corr(j-i,l)=corr(j-i,l)+1.0D0
                              END IF
                          END IF
                      END DO
                  END DO
              END DO
              label(l)=.TRUE.

100           CONTINUE
          END DO
      END DO

!=======================================================================
!----- Write residence times and correlation ---------------------------
!=======================================================================
    
      tn=ncorr
      IF(tn .GT. nmax) tn=nmax
      WRITE(ksol1,500) time
      nc=0
      DO n=1,annpro
          nstart=anpoint(1,n)
          nend=anpoint(2,n)
          mstart=mres(1,res(nstart,1))
          mfin=mres(2,res(nend,1))
          DO l=mstart,mfin
              i=grppt(1,l)
              IF(label(l)) THEN
                  pk=0
                  DO k=1,nmax
                      pk=pk+wnlst(l,k)
                  END DO

                  DO k=0,ncorr
                      corr(k,l)=corr(k,l)/DBLE(ncc(k,l))
                  END DO


                  DO k=1,ncorr
                      tcorr(k)=corr(k,l)/corr(0,l)
                      tcorr(k)=-time*DBLE(k)/DLOG(tcorr(k))
                  END DO
                  IF(ksol2 .NE. -1) THEN
                      WRITE(ksol2,2500) res(grppt(1,l),1),l
                      DO k=0,ncorr
                          IF(.NOT. near0(corr(k,l)))
     x                        WRITE(ksol2,3000) time*DBLE(k),corr(k,l)
                      END DO
                  END IF

                  dermin=1.0D+10
                  idermin=-1
                  DO k=1,ncorr-1
                      IF(k .NE. 0 .AND. (.NOT. near0(corr(k,l))) ) THEN
                          der(k-1)=(tcorr(k)-tcorr(k-1))/time
                          IF(dermin .GT. DABS(der(k-1)) ) THEN
                              dermin=DABS(der(k-1))
                              idermin=k-1
                          END IF
                      END IF
                  END DO
                  k=grppt(1,l)
                  IF(idermin .NE. -1) THEN
                      WRITE(ksol1,1000) prsymb(res(k,2)),
     x                  res(k,1),res(k,1)-nc,corr(0,l),l,
     x                  tcorr(idermin),der(idermin),
     x                  (beta(i),i=grppt(1,l),grppt(2,l))
                  ELSE
                      WRITE(ksol1,2000) prsymb(res(k,2)),
     x                  res(k,1),res(k,1)-nc,corr(0,l),l,
     x                  (beta(i),i=grppt(1,l),grppt(2,l))
                  END IF
              END IF
          END DO
          nc=res(grppt(2,mfin),1)
      END DO

!=======================================================================
!----- Write waters only bound more than 90 % of the time --------------
!=======================================================================

      WRITE(ksol1,5500) 
      nc=0
      DO i=1,nbd1
          n=boundww(i)
          m=boundrr(i)
          ok=.TRUE.
          DO j=i-1,1,-1
              IF(n .EQ. boundww(j)) ok=.FALSE.
          END DO
          IF(ok) THEN
              nc=nc+1
              WRITE(ksol1,6000) nc,prsymb(res(grppt(1,m),2)),
     x              res(grppt(1,m),1),n,m
          END IF
      END DO

!=======================================================================
!----- Write bound waters ----------------------------------------------
!=======================================================================

      WRITE(ksol1,5000) 
      nc=0
      DO i=1,nbd
          n=boundw(i)
          m=boundr(i)
          ok=.TRUE.
          DO j=i-1,1,-1
              IF(n .EQ. boundw(j)) ok=.FALSE.
          END DO
          IF(ok) THEN
              nc=nc+1
              WRITE(ksol1,6000) nc,prsymb(res(grppt(1,m),2)),
     x              res(grppt(1,m),1),n,m
          END IF
      END DO
      
500   FORMAT(' TIME = ',f12.5)
1000  FORMAT('RESIDUE ',a5,2i5,' AVG WATER =',f12.5,'  GROUP ',i6,
     x ' TAU = ',f13.5,'  DER = ',e12.5,' ATOMS  ',15a5)
2000  FORMAT('RESIDUE ',a5,2i5,' AVG WATER =',f12.5,'  GROUP ',i6,
     x ' TAU =  NA ',' ATOMS  ',15a5)
2500  FORMAT('CORRELATION ',2i5)
3000  FORMAT(f12.5,f12.7)
5000  FORMAT('BOUND WATER MOLECULES: ',
     x'RESIDE ON THE RESIDUE ALL THE TIME')
6000  FORMAT(i6,2x,'RESIDUE ',a5,i6,'  WATER NUM. = ',i6,
     x' GROUP NUM. = ',i6)
5500  FORMAT('WATER MOLECULES BOUND TO A RESIDUE ',
     x'MORE THAN 90 % OF THE TIME')

      RETURN
      END
