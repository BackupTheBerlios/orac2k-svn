      SUBROUTINE comp_abmd_fnative(unbias,flag,iflag,nato,xp0,yp0,zp0
     &     ,vpx,vpy,vpz,co,atomp,xpcm,ypcm,zpcm,uconf,enout,omega,alpha
     &     ,atres_map1,atres_map2,listp,list,native_dist,gr,gra,dfree,fx
     &     ,fy,fz,stress,nstart,nend,nlocal,node,nprocs,ncube)

************************************************************************
*                                                                      *
*     Umbrella subroutine.                                             *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*---- Last update 03/10/91 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi UCB, Berkeley 1991                     *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      
      REAL*8  uconf,enout,omega,alpha,dfree,co(3,3),xp0(*),yp0(*),zp0(*)
     &     ,xpcm(*),ypcm(*),zpcm(*),stress(3,3)
      REAL*8  vpx(*),vpy(*),vpz(*),fx(*),fy(*),fz(*)
      REAL*8  gr,gra,native_dist
      INTEGER nstart,nend,nlocal,node,nprocs,ncube,nato
     &     ,iflag,atres_map1(*),atres_map2(*),listp,list(2,*),atomp(*)
      LOGICAL flag,unbias

*-------------------- COMMON VARIABLES ---------------------------------

      include 'parst.h'

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,n,m,nbtot,initial,count1,count2,ma,mb,ia,ii,jj
      REAL*8  upper,lower,coef,grp,fnlistp,rsq,dexpa,dexpb
      REAL*8  fpx,fpy,fpz,sumix,sumiy,sumiz,sumjx,sumjy,sumjz,work
      LOGICAL near0
      INTEGER mma,nstart_1,nend_1,nlocal_1,listp_1,idlt
      PARAMETER (mma=300)
      COMMON /rag1/ fpx(m1),fpy(m1),fpz(m1),sumix(mma),sumiy(mma)
     &     ,sumiz(mma),sumjx(mma),sumjy(mma),sumjz(mma),work(m1)
      DATA initial/0/
      SAVE initial
      REAL*8  cut2,xpi,ypi,zpi,xc,yc,zc,xd,yd,zd,rs,r,al,sum,nat,auxa
     &     ,auxb,sgx,sgy,sgz,sx,sy,sz
      REAL*8 stc1,stc2,stc3,stc4,stc5,stc6,stc7,stc8,stc9,dlt
      LOGICAL ok

*-------------------- LOCAL VARIABLES ----------------------------------

      DO i=1,nato
         fpx(i)=0.0D0
         fpy(i)=0.0D0
         fpz(i)=0.0D0
      END DO

      uconf=0.0D0
      gr=0.0D0
      gra=0.0D0

      count1=0
      count2=0
      cut2=native_dist**2
      fnlistp=1.0D0/DFLOAT(listp)

      listp_1=listp/nprocs
      nstart_1=listp_1*node+1
      nend_1=listp_1*(node+1)
      IF(node .EQ. nprocs-1) nend_1=listp
      nlocal_1=nend_1-nstart_1+1

      DO ia=1,nstart_1-1
         ma=atres_map1(count1+1)
         mb=atres_map2(count2+1)
         count1=count1+ma+1
         count2=count2+mb+1
      END DO

      DO ia=nstart_1,nend_1
         ii=list(1,ia)
         jj=list(2,ia)
         ok=.FALSE.
         sum=1.0D0
         ma=atres_map1(count1+1)
         mb=atres_map2(count2+1)
         DO ii=1,ma
            sumix(ii)=0.0D0
            sumiy(ii)=0.0D0
            sumiz(ii)=0.0D0
         END DO
         DO jj=1,mb
            sumjx(jj)=0.0D0
            sumjy(jj)=0.0D0
            sumjz(jj)=0.0D0
         END DO
         DO ii=1,ma
            i=atres_map1(count1+1+ii)
            xpi=xp0(i)
            ypi=yp0(i)
            zpi=zp0(i)
            DO jj=1,mb
               j=atres_map2(count2+1+jj)
               xc=xpi-xp0(j)
               yc=ypi-yp0(j)
               zc=zpi-zp0(j)
               rs=xc*xc+yc*yc+zc*zc
               rsq=DSQRT(rs)
               r=(rsq-native_dist)*alpha
               rsq=1.0D0/rsq
               IF(r .LT. 30.0D0) THEN
                  dexpa=DEXP(r)
                  dexpb=DEXP(-r)
               ELSE
                  dexpa=1.0D20
                  dexpb=1.0D-20
               END IF
               auxa=0.5D0+0.5D0*(dexpa-dexpb)/(dexpa+dexpb)
               sum=sum*auxa
               auxb=2.0D0*alpha/(dexpa+dexpb)**2
               auxb=rsq*auxb/auxa
               sumix(ii)=sumix(ii)+xc*auxb
               sumiy(ii)=sumiy(ii)+yc*auxb
               sumiz(ii)=sumiz(ii)+zc*auxb
               sumjx(jj)=sumjx(jj)+xc*auxb
               sumjy(jj)=sumjy(jj)+yc*auxb
               sumjz(jj)=sumjz(jj)+zc*auxb
            END DO
         END DO
         DO ii=1,ma
            i=atres_map1(count1+1+ii)
            auxa=sum*fnlistp
            fpx(i)=fpx(i)-sumix(ii)*auxa
            fpy(i)=fpy(i)-sumiy(ii)*auxa
            fpz(i)=fpz(i)-sumiz(ii)*auxa
         END DO
         DO jj=1,mb
            j=atres_map2(count2+1+jj)
            auxa=sum*fnlistp
            fpx(j)=fpx(j)+sumjx(jj)*auxa
            fpy(j)=fpy(j)+sumjy(jj)*auxa
            fpz(j)=fpz(j)+sumjz(jj)*auxa
         END DO
         gr=gr+1.0D0-sum
         gra=gra+DNINT(1.0D0-sum)
         count1=count1+ma+1
         count2=count2+mb+1
      END DO

#ifdef PARALLEL
      CALL P_merge_r8(gr)
      CALL P_merge_r8(gra)
#endif

      gr=gr*fnlistp
      gra=gra*fnlistp
      grp=gr
      IF(near0(enout)) THEN
         enout=grp
      END IF
      IF(initial .EQ. 0) THEN
         IF(node .EQ. 0) WRITE(*,1000) enout
      END IF

#ifdef PARALLEL
      CALL P_fold_r8(nato,fpx,nstart,nend,nlocal,node,nprocs)
      CALL P_fold_r8(nato,fpy,nstart,nend,nlocal,node,nprocs)
      CALL P_fold_r8(nato,fpz,nstart,nend,nlocal,node,nprocs)
#endif

      IF(flag .OR. unbias) THEN
         IF(.NOT. unbias) THEN
            IF(iflag .LT. 0) THEN
            
*=======================================================================
*--- If the solute must be folded --------------------------------------
*=======================================================================

               IF(enout .GE. grp) THEN
                  enout=grp
                  coef=0.0D0
                  uconf=0.0D0
               ELSE
                  coef=2.0D0*omega*(enout-grp)
                  uconf=omega*(enout-grp)**2
               END IF
               DO i=nstart,nend
                  fpx(i)=coef*fpx(i)
                  fpy(i)=coef*fpy(i)
                  fpz(i)=coef*fpz(i)
               END DO

*=======================================================================
*--- If the solute must be unfolded ------------------------------------
*=======================================================================

            ELSE
               IF(enout .LE. grp) THEN
                  enout=grp
                  coef=0.0D0
                  uconf=0.0D0
               ELSE
                  coef=2.0D0*omega*(enout-grp)
                  uconf=omega*(enout-grp)**2
               END IF
               DO i=nstart,nend
                  fpx(i)=coef*fpx(i)
                  fpy(i)=coef*fpy(i)
                  fpz(i)=coef*fpz(i)
               END DO
            END IF
         ELSE
            coef=2.0D0*omega*(enout-grp)
            uconf=omega*(enout-grp)**2
            DO i=nstart,nend
               fpx(i)=coef*fpx(i)
               fpy(i)=coef*fpy(i)
               fpz(i)=coef*fpz(i)
            END DO
         END IF
      END IF
      gr=grp
      initial=1
      
      stc1=0.0D0
      stc2=0.0D0
      stc3=0.0D0
      stc4=0.0D0
      stc5=0.0D0
      stc6=0.0D0
      stc7=0.0D0
      stc8=0.0D0
      stc9=0.0D0
      DO i=nstart,nend
         j=atomp(i)
         sx=xpcm(j)
         sy=ypcm(j)
         sz=zpcm(j)
         sgx=sx*co(1,1)+sy*co(1,2)+sz*co(1,3)
         sgy=sx*co(2,1)+sy*co(2,2)+sz*co(2,3)
         sgz=sx*co(3,1)+sy*co(3,2)+sz*co(3,3)
         stc1=stc1+sgx*fpx(i)
         stc2=stc2+sgy*fpx(i)
         stc3=stc3+sgz*fpx(i)
         stc4=stc4+sgx*fpy(i)
         stc5=stc5+sgy*fpy(i)
         stc6=stc6+sgz*fpy(i)
         stc7=stc7+sgx*fpz(i)
         stc8=stc8+sgy*fpz(i)
         stc9=stc9+sgz*fpz(i)
      END DO
#ifdef PARALLEL
      CALL P_merge_r8(stc1)
      CALL P_merge_r8(stc2)
      CALL P_merge_r8(stc3)
      CALL P_merge_r8(stc4)
      CALL P_merge_r8(stc5)
      CALL P_merge_r8(stc6)
      CALL P_merge_r8(stc7)
      CALL P_merge_r8(stc8)
      CALL P_merge_r8(stc9)
#endif
      stress(1,1)=stress(1,1)+stc1
      stress(1,2)=stress(1,2)+stc2
      stress(1,3)=stress(1,3)+stc3
      stress(2,1)=stress(2,1)+stc4
      stress(2,2)=stress(2,2)+stc5
      stress(2,3)=stress(2,3)+stc6
      stress(3,1)=stress(3,1)+stc7
      stress(3,2)=stress(3,2)+stc8
      stress(3,3)=stress(3,3)+stc9

      DO i=nstart,nend
         fx(i)=fx(i)+fpx(i)
         fy(i)=fy(i)+fpy(i)
         fz(i)=fz(i)+fpz(i)
      END DO

      CALL comp_noneq_free(vpx,vpy,vpz,fpx,fpy,fpz,nstart,nend,dfree)
#ifdef PARALLEL
      CALL P_merge_r8(dfree)
#endif

*================= END OF EXECUTABLE STATEMENTS ========================

1000  FORMAT(/ /'----------- N_c = ', f12.4,' -------------'/ /)
      RETURN
      END
