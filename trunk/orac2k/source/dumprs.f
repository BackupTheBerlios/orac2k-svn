      SUBROUTINE dumprs(kdump,restart_file,nstep,temp,ntap,ngrp,nprot
     &     ,xp0,yp0,zp0,vpx,vpy,vpz,gh,vh,nh,hoover,avg1,avg2,na1,na2,co
     &     ,oc,cpress,zz,igr,krdf,igrn)

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

      INTEGER nstep,ntap,kdump,ngrp,nprot
      INTEGER nh,na1,na2,igrn,krdf(*)
      LOGICAL hoover,cpress,igr
      REAL*8  xp0(*),yp0(*),zp0(*),vpx(*),vpy(*),vpz(*),gh(*),vh(*)
     &     ,avg1(*),avg2(*),co(3,3),oc(3,3),zz(3,3),temp
      CHARACTER*80 restart_file

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,k
      LOGICAL  exist

*==================== EXECUTABLE STATEMENTS ============================

      INQUIRE(FILE=restart_file,EXIST=exist)
      IF(exist) THEN
         CALL openf(kdump,restart_file,'UNFORMATTED','OLD',0)
      ELSE
         CALL openf(kdump,restart_file,'UNFORMATTED','NEW',0)
      END IF
      REWIND kdump
      WRITE(kdump) ((co(i,j),j=1,3),i=1,3),((oc(i,j),j=1,3),i=1,3)
      WRITE(kdump) nstep,temp
      WRITE(kdump) ntap,ngrp,nprot
      WRITE(kdump) (xp0(i),i=1,ntap),(yp0(i),i=1,ntap),
     &     (zp0(i),i=1,ntap),(vpx(i),i=1,ntap),(vpy(i),i=1,ntap)
     &     ,(vpz(i),i=1,ntap)
      WRITE(kdump) na1,na2
      WRITE(kdump) (avg1(i),i=1,na1),(avg2(i),i=1,na2)
      WRITE(kdump) hoover,cpress

      IF(hoover) THEN
         WRITE(kdump) nh
         WRITE(kdump) (gh(i),i=1,nh),(vh(i),i=1,nh)
      END IF
      IF(cpress) THEN
         WRITE(kdump) ((zz(i,j),j=1,3),i=1,3)
      END IF

      WRITE(kdump) igr
      IF(igr) THEN
         WRITE(kdump) igrn
         WRITE(kdump) (krdf(i),i=1,igrn)
      END IF
      CLOSE(kdump)

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
