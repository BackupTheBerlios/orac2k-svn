      SUBROUTINE write_hbonds(update_anl,hbonds_tot,hbonds_res,fstep
     &     ,khbonds,ss_index,nbun,nres,o1,atomp,co,xp0,yp0,zp0,lacc,ldon
     &     ,llacc,lldon,hrcut,hacut,ha2cut,nmapdn,mapdn,hhisto
     &     ,hhisto_list,hhisto_bin,hhisto_dim,hhisto_count,beta,prsymb
     &     ,read_pdb)

************************************************************************
*   Time-stamp: <99/09/27 20:18:56 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed Jul  9 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER o1,khbonds,nbun,update_anl,hhisto_dim,hhisto_count
      INTEGER hhisto_list(3,*),ss_index(*),nres(o1,*),llacc,lldon,lacc(2
     &     ,*),ldon(2,*),atomp(*),mapdn(*),nmapdn(*)
      REAL*8 hrcut,hacut,ha2cut,xp0(*),yp0(*),zp0(*),co(3,3),fstep
     &     ,hhisto_bin
      LOGICAL hbonds_tot,hbonds_res,hhisto,read_pdb
      CHARACTER*7 beta(*)
      CHARACTER*8 prsymb(*)

*----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'parst.h'
      INTEGER count_res(nores),index(nores),mapp(m1)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,ia,ib,j,ja,jb,typei,typeij,type_resi,type_resj,count(3)
     &     ,map,n,m,counta,jj,ibin,GET_IND
      REAL*8  xr1,xr2,xr3,xr4,yr1,yr2,yr3,yr4,zr1,zr2,zr3,zr4,x21,x32
     &     ,x43,x42,y21,y32,y43,y42,z21,z32,z43,z42,xc32,xc43,xc42
     &     ,yc32,yc43,yc42,zc32,zc43,zc42,rsp32,rsp43,rsq32,rsq43
     &     ,rsq42,cb234,hrcut2,xc21,yc21,zc21,rsp21,rsq21,cb123,y
     &     ,x1,pi
      LOGICAL ok
      COMMON /rag1/ xr1,xr2,xr3,xr4,yr1,yr2,yr3,yr4,zr1,zr2,zr3,zr4,x21
     &     ,x32,x43,x42,y21,y32,y43,y42,z21,z32,z43,z42,xc32,xc43,xc42
     &     ,yc32,yc43,yc42,zc32,zc43,zc42,rsp32,rsp43,rsq32,rsq43,rsq42
     &     ,cb234,hrcut2,xc21,yc21,zc21,rsp21,rsq21,cb123,y,x1,pi
     &     ,count_res,index,mapp
      INCLUDE 'pbc.h'
      GET_IND(y)=INT(y)+(SIGN(1.0D0,y-INT(y))-1.0D0)/2+1

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      IF(read_pdb) WRITE(khbonds,100)
      pi=4.0D0*DATAN(1.0D0)
      hrcut2=hrcut**2
      DO i=1,nbun
         count_res(i)=0
      END DO      
      DO i=1,3
         count(i)=0
      END DO
      counta=0
      DO i=1,llacc
         ia=lacc(1,i)
         ib=lacc(2,i)
         typei=ss_index(ib)
         type_resi=nres(ib,1)
         xr1=xp0(ia)
         yr1=yp0(ia)
         zr1=zp0(ia)
         xr2=xp0(ib)   
         yr2=yp0(ib)
         zr2=zp0(ib)
         IF(update_anl .NE. 0) THEN
            m=nmapdn(i)
            DO j=1,m
               mapp(j)=mapdn(counta+j)
            END DO
         ELSE
            m=lldon
            DO j=1,m
               mapp(j)=j
            END DO
         END IF
         
         DO jj=1,m
            j=mapp(jj)
            ja=ldon(1,j)
            jb=ldon(2,j)
            typeij=typei+ss_index(ja)-1
            type_resj=nres(ja,1)
            IF(IABS(type_resj-type_resi) .GT. 1 .OR. (atomp(ia) .NE.
     &           atomp(ja))) THEN
               xr3=xp0(ja)
               yr3=yp0(ja)
               zr3=zp0(ja)   
               xr4=xp0(jb)   
               yr4=yp0(jb)
               zr4=zp0(jb)

               x32=xr3-xr2
               y32=yr3-yr2
               z32=zr3-zr2
               x32=x32-2.0D0*PBC(x32)
               y32=y32-2.0D0*PBC(y32)
               z32=z32-2.0D0*PBC(z32)
               xc32=co(1,1)*x32+co(1,2)*y32+co(1,3)*z32
               yc32=            co(2,2)*y32+co(2,3)*z32
               zc32=                        co(3,3)*z32
               rsq32=xc32**2+yc32**2+zc32**2
               IF(rsq32 .LE. hrcut2) THEN
                  ok=.TRUE.
                  IF(DABS(hacut-1.0D0) .GT. 1.0D-6) THEN
                     x43=xr4-xr3
                     y43=yr4-yr3
                     z43=zr4-zr3
                     x43=x43-2.0D0*PBC(x43)
                     y43=y43-2.0D0*PBC(y43)
                     z43=z43-2.0D0*PBC(z43)
                     xc43=co(1,1)*x43+co(1,2)*y43+co(1,3)*z43
                     yc43=            co(2,2)*y43+co(2,3)*z43
                     zc43=                        co(3,3)*z43
                     rsq32=xc32**2+yc32**2+zc32**2
                     rsq43=xc43**2+yc43**2+zc43**2
                     rsp32=DSQRT(rsq32)
                     rsp43=DSQRT(rsq43)
                     cb234=-(xc43*xc32+yc43*yc32+zc43*zc32)/(rsp43*rsp32
     &                    )
                     
*=======================================================================
*---- 123 = AA-A-H        234 = A-H-D ----------------------------------
*=======================================================================

                     IF(cb234 .GT. hacut) ok=.FALSE.
                  END IF
                  IF(DABS(ha2cut-1.0D0) .GT. 1.0D-6 .AND. ia .NE. ib)
     &                 THEN
                     x21=xr2-xr1
                     y21=yr2-yr1
                     z21=zr2-zr1
                     x21=x21-2.0D0*PBC(x21)
                     y21=y21-2.0D0*PBC(y21)
                     z21=z21-2.0D0*PBC(z21)
                     xc21=co(1,1)*x21+co(1,2)*y21+co(1,3)*z21
                     yc21=            co(2,2)*y21+co(2,3)*z21
                     zc21=                        co(3,3)*z21
                     rsq21=xc21**2+yc21**2+zc21**2
                     rsp21=DSQRT(rsq21)
                     IF(ia .NE. ib) THEN
                        cb123=-(xc21*xc32+yc21*yc32+zc21*zc32)/(rsp21
     &                       *rsp32)
                        IF(cb123 .GT. ha2cut) ok=.FALSE.
                     ELSE
                        ok=.FALSE.
                     END IF
                  END IF
                  IF(read_pdb) THEN
                     WRITE(khbonds,'(2(i3,1x,a8,1x,a7,1x),3f12.4)')
     &                    nres(ia,1),prsymb(nres(ia,2)),beta(ib)
     &                    ,nres(ja,1),prsymb(nres(ja,2)),beta(ja)
     &                    ,DSQRT(rsq32),DACOS(cb234)*180.0D0/pi
     &                    ,DACOS(cb123)*180.0D0/pi
                  END IF
                  IF(ok .AND. (.NOT. read_pdb)) THEN
                     count(typeij)=count(typeij)+1
                     count_res(type_resi)=count_res(type_resi)+1
                     count_res(type_resj)=count_res(type_resj)+1
                     IF(hhisto) THEN
                        x1=DSQRT(rsq32)/hhisto_bin
                        ibin=GET_IND(x1)
                        hhisto_list(typeij,ibin)=hhisto_list(typeij,ibin
     &                       )+1
                     END IF
                  END IF
               END IF
            END IF
         END DO
         counta=counta+m
      END DO
      
      IF(.NOT. read_pdb) THEN
         IF(hbonds_tot)
     &        WRITE(khbonds,'(a1,f11.2,a16,i7,a9,i7,a9,i7)') 'T',fstep,
     &        ' Hbonds SLT-SLT ',count(1),' SLV-SLV ',count(3)
     &        ,' SLT-SLV ',count(2)
         
         IF(hbonds_res) THEN
            map=0
            DO i=1,nbun
               IF(count_res(i) .GT. 0) THEN
                  map=map+1
                  index(map)=i
               END IF
            END DO
            DO j=1,map,4
               n=3
               IF(map-j .LT. 3) n=map-j
               WRITE(khbonds,'(a1,f11.2,1x,a8,4(1x,i6,1x,i6))') 'T'
     &              ,fstep,' Hbonds ',(count_res(index(i)),index(i),i=j
     &              ,j+n)
            END DO
         END IF
      END IF
      IF(hhisto) hhisto_count=hhisto_count+1
         
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

100   FORMAT(' NA RESA     ACC      ND RESD     DON   '
     &     ,'        A...H     A...H--D   AA--A...H '/
     &     '******************************************'
     &     ,'************************************')
      RETURN
      END
