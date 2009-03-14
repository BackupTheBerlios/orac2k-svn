      SUBROUTINE read_confc_rows(co,xp0,yp0,zp0,xau,yau,zau,nato,fstep
     &     ,record,end,err,divide_records,atom_record)

************************************************************************
*                                                                      *
*                                                                      *
*     CONF will dump the hystory of the MD run. The writing on         *
*     disk is performed by the IOLIB routines.                         *
*                                                                      *
*     XP0     :  Solute molecules site coordinates, packed        (I)  *
*     YP0        by molecule.                                          *
*     ZP0        >> real*8  XP0(*), ... <<                             *
*                                                                      *
*     NSTEP   :  Number of the current configuration being        (I)  *
*                dumped.                                               *
*                                                                      *
*---- Last update 05/03/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNALS IOLIB                                                  *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

#if defined T3E | defined _CRAY_
#define _REAL_   REAL*8
#else
#define _REAL_   REAL*4
#endif
      INTEGER nato,record,divide_records,atom_record
      REAL*8  xp0(*),yp0(*),zp0(*),fstep
      _REAL_  xau(*),yau(*),zau(*),faux
      REAL*8  co(3,3)
      LOGICAL end,err
      
*------------------- VARIABLES IN COMMON -------------------------------

      INCLUDE 'unit.h'

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,rec,count,unit,reca,nstart,nend

*==================== EXECUTABLE STATEMENTS ============================

      end=.FALSE.
      err=.FALSE.
*=======================================================================
*----- Read box and timestep -------------------------------------------
*=======================================================================

      rec=record
      unit=kwrite_dump_i(1)
      READ(unit=unit,rec=rec,ERR=100) faux,((co(i,j),j=1,3),i=1,3)

*=======================================================================
*----- Determines which file to read from ------------------------------
*=======================================================================

      rec=divide_records*(record-1)+1
      count=0
      i=1
      DO WHILE(count .LT. rec)
         i=i+1
         count=count+no_records_i(i)
      END DO

      unit=kwrite_dump_i(i)

      DO j=2,i-1
         rec=rec-no_records_i(j)
      END DO

*=======================================================================
*----- Dump the coordinates of the solute molecules when required ------
*=======================================================================

      DO i=1,divide_records
         nstart=(i-1)*atom_record+1
         nend=nstart-1+atom_record
         IF(nend .GT. nato) nend=nato
         reca=rec+i-1
         READ(unit=unit,rec=reca,ERR=200) (xau(j),yau(j),zau(j),j
     &        =nstart,nend)
      END DO

      DO i=1,nato
         xp0(i)=xau(i)
         yp0(i)=yau(i)
         zp0(i)=zau(i)
      END DO
      fstep=faux
      RETURN

100   end=.TRUE.
      RETURN
200   err=.TRUE.

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
