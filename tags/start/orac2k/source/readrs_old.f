      SUBROUTINE readrs_old(kdump,file,nstep,temp,ntap,ngrp,nprot,xp0
     &     ,yp0,zp0,vpx,vpy,vpz,xpg,ypg,zpg,xcm,ycm,zcm,vpx1,vpy1,vpz1
     &     ,vcx,vcy,vcz,gh,vh,vh1,nh,hoover,avg1,avg2,na1,na2,co,oc
     &     ,cpress,zz,zz1,lzz,lzz1,igr,krdf,igrn,igrnmax)

************************************************************************
*                                                                      *
*                                                                      *
*     READRS will read all data necessary to reastrt an MD run.        *
*                                                                      *
*     N1      :  Number of elements of AVG1 read by READRS.       (O)  *
*     N2      :  Number of elements of AVG2 read by READRS.       (O)  *
*     M1      :  Physical row dimension og IGR1 and IGR2.         (I)  *
*     IG1     :  Number of columns  of IGR1 read by READRS.       (O)  *
*     IG2     :  Number of columns  of IGR2 read by READRS.       (O)  *
*     MAXINT  :  Number of rows of IGR2 read by READRS.           (O)  *
*                                                                      *
*----- Last update 06/5/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNALS none.                                                  *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER nstep,ntap,kdump,ngrp,m1,nprot
      INTEGER nh,na1,na2,igrn,igrnmax
      LOGICAL hoover,cpress,igr,krdf(*)
      REAL*8  xp0(*),yp0(*),zp0(*),vpx(*),vpy(*),vpz(*),xpg(*),ypg(*)
     &     ,zpg(*),vpx1(*),vpy1(*),vpz1(*),gh(*),vh(*),vh1(*),avg1(*)
     &     ,avg2(*),co(3,3),oc(3,3),zz(3,3),zz1(3,3),lzz,lzz1,temp,xcm(
     &     *),ycm(*),zcm(*),vcx(*),vcy(*),vcz(*)
      CHARACTER*80 file

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,k,n1,n2
      LOGICAL igr1

*==================== EXECUTABLE STATEMENTS ============================

      OPEN(unit=kdump,file=file,form='UNFORMATTED',status='OLD')
      REWIND kdump
      READ(kdump,END=100) ((co(i,j),j=1,3),i=1,3),((oc(i,j),j=1,3),i=1,3
     &     )
      READ(kdump,END=100,ERR=100) nstep,temp

      READ(kdump,END=100,ERR=100) ntap,ngrp,nprot
      READ(kdump,END=100,ERR=100) (xp0(i),i=1,ntap),(yp0(i),i=1,ntap)
     &     ,(zp0(i),i=1,ntap),(vpx(i),i=1,ntap),(vpy(i),i=1,ntap)
     &     ,(vpz(i),i=1,ntap),(xpg(i),i=1,ngrp),(ypg(i),i=1,ngrp)
     &     ,(zpg(i),i=1,ngrp),(xcm(i),i=1,nprot),(ycm(i),i=1,nprot)
     &     ,(zcm(i),i=1,nprot),(vcx(i),i=1,nprot),(vcy(i),i=1,nprot)
     &     ,(vcz(i),i=1,nprot)
      READ(kdump,END=100,ERR=100) (vpx1(i),i=1,ntap),(vpy1(i),i=1
     &        ,ntap),(vpz1(i),i=1,ntap)

      READ(kdump,END=100,ERR=100) n1,n2

      IF(na1 .LT. n1 .OR. na2 .LT. n2) THEN
         WRITE(6,'(a)') 'In READRS: Dimension of the averages array'
     &        //' insufficient. Abort.'
         STOP
      END IF
      READ(kdump,END=100,ERR=100) (avg1(i),i=1,n1),(avg2(i),i=1,n2)

      READ(kdump,ERR=200,END=200) hoover,cpress
      IF(hoover) THEN
         READ(kdump,END=100,ERR=100) nh
         READ(kdump,END=100,ERR=100) (gh(i),i=1,nh),(vh(i),i=1,nh)
     &        ,(vh1(i),i=1,nh)
      END IF
      IF(cpress) THEN
         READ(kdump,END=100,ERR=100) ((zz(i,j),j=1,3),i=1,3),((zz1(i,j)
     &        ,j=1,3),i=1,3),lzz,lzz1
      END IF

      READ(kdump,END=100,ERR=100) igr1
      IF(igr .AND. igr1) THEN
         READ(kdump,END=100,ERR=100) igrn
         IF(igrn .LE. igrnmax) THEN
            READ(kdump,END=100,ERR=100) (krdf(i),i=1,igrn)
         ELSE
            GOTO 300
         END IF
      END IF

      GOTO 200
100   CONTINUE
      WRITE(6,*)
      WRITE(6,'(a)') ' End of file encountered. Wrong restart'//
     &     ' file or file has been corrupted. Abort.'
      WRITE(6,*)
      STOP
200   CONTINUE

      CLOSE(kdump)

      RETURN

300   CONTINUE
      WRITE(6,*)
      WRITE(6,'(a)') ' The G(r) is read in incorrectly. I hope it is'//
     &               ' ok. Run continues'
      WRITE(6,*)

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
