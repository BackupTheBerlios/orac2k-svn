      SUBROUTINE write_confc(co,xp0,yp0,zp0,nato,fstep,ninner,nconf
     &     ,divide_records,atom_record)

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

      INTEGER ninner,nconf,nato,divide_records,atom_record
      REAL*8  xp0(*),yp0(*),zp0(*),fstep
      REAL*8  co(3,3)

*------------------- VARIABLES IN COMMON -------------------------------

      INCLUDE 'unit.h'

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,rec,count,unit,nstart,nend,reca,recl
      LOGICAL exist
      CHARACTER*80 form,file

*==================== EXECUTABLE STATEMENTS ============================

*=======================================================================
*----- Write box and timestep ------------------------------------------
*=======================================================================

#if defined OPTERON
      DO i=1,nwrite_dump_o
         INQUIRE(FILE=file_names_o(i),EXIST=exist)
         unit=kwrite_dump_o(i)
         file=file_names_o(i)
         form='UNFORMATTED'
         CLOSE(unit)
         IF(i .NE. 1) THEN
            recl=recl2_dump_o
            IF(exist) THEN
               OPEN(unit=unit,file=file,access='DIRECT',form=form,status
     &              ='OLD',recl=recl)
            ELSE
               OPEN(unit=unit,file=file,access='DIRECT',form=form,status
     &              ='NEW',recl=recl)
            END IF
         ELSE 
            recl=recl1_dump_o
            IF(exist) THEN
               OPEN(unit=unit,file=file,access='DIRECT',form=form,status
     &              ='OLD',recl=recl)
            ELSE
               OPEN(unit=unit,file=file,access='DIRECT',form=form,status
     &              ='NEW',recl=recl)
            END IF
         END IF
      END DO
#endif
      rec=ninner/nconf
      unit=kwrite_dump_o(1)
#if defined T3E | defined _CRAY_
      WRITE(unit=unit,rec=rec) fstep,((co(i,j),j=1,3),i=1,3)
#else
      WRITE(unit=unit,rec=rec) SNGL(fstep),((co(i,j),j=1,3),i=1,3)
#endif

*=======================================================================
*----- Determines which file to write to -------------------------------
*=======================================================================

      rec=divide_records*(ninner/nconf-1)+1
      count=0
      i=1
      DO WHILE(count .LT. rec)
         i=i+1
         count=count+no_records_o(i)
      END DO

      unit=kwrite_dump_o(i)
      DO j=2,i-1
         rec=rec-no_records_o(j)
      END DO

*=======================================================================
*----- Dump the coordinates of the solute molecules when required ------
*=======================================================================

      DO i=1,divide_records
         nstart=(i-1)*atom_record+1
         nend=nstart-1+atom_record
         IF(nend .GT. nato) nend=nato
         reca=rec+i-1
#if defined T3E | defined _CRAY_
         WRITE(unit=unit,rec=reca) (xp0(j),yp0(j),zp0(j),j=nstart
     &        ,nend)
#else
         WRITE(unit=unit,rec=reca) (SNGL(xp0(j)),SNGL(yp0(j))
     &        ,SNGL(zp0(j)),j=nstart,nend)
#endif
      END DO

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
