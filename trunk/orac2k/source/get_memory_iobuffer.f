      SUBROUTINE get_memory_iobuffer(start_anl,stop_anl,start_time
     &     ,end_time,length_run,atom_record,length_tot,length_fft)

************************************************************************
*   Time-stamp: <2009-03-09 12:05:30 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed Dec  3 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none
#include  "iobuffer.h"
*----------------------------- ARGUMENTS ------------------------------*

      

      INTEGER start_anl,stop_anl,length_run,length_tot
     &     ,length_fft,start_time,end_time,iret
     &     ,atom_record
      CHARACTER*80 errmsg

*------------------------- LOCAL VARIABLES ----------------------------*

      CHARACTER*8 string
      INTEGER len,M_get_length

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      iret=0
      start_time=start_anl
      end_time=stop_anl
      length_run=stop_anl
      length_tot=end_time-start_time+1
      buffer_time=length_tot
      nbuffer=length_tot*atom_record
      IF(MOD(length_tot,2) .NE. 0) THEN
         end_time=end_time-1
         length_tot=length_tot-1
      END IF
      length_fft=length_tot*2
      buffer_fft=length_fft

      ALLOCATE(xpb(nbuffer))
      ALLOCATE(ypb(nbuffer))
      ALLOCATE(zpb(nbuffer))
      ALLOCATE(xpc(buffer_time))
      ALLOCATE(ypc(buffer_time))
      ALLOCATE(zpc(buffer_time))
      ALLOCATE(spline_x(buffer_time))
      ALLOCATE(spline_y(4,buffer_time))
      ALLOCATE(xpcc(buffer_time))
      ALLOCATE(ypcc(buffer_time))
      ALLOCATE(zpcc(buffer_time))
      ALLOCATE(RotMat(3,3,buffer_time))

      ALLOCATE(vxpc(buffer_fft))
      ALLOCATE(vypc(buffer_fft))
      ALLOCATE(vzpc(buffer_fft))

      ALLOCATE(vacf_data(0:buffer_fft*3))
      ALLOCATE(vacf_spectra(0:buffer_fft*3))
      ALLOCATE(rms_disp(0:buffer_fft*2))
      ALLOCATE(tot_rms_disp(0:buffer_fft))

      ALLOCATE(wsave1(buffer_fft*4+15))
      ALLOCATE(w1(buffer_fft*2))
      ALLOCATE(w2(buffer_fft*2))

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
