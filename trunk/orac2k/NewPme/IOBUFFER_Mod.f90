MODULE IOBUFFER_Mod


!!$*------------------ PARAMETER  STATEMENT -------------------------------

  IMPLICIT none
  INTEGER(4), SAVE ::  nbuffer,buffer_time,buffer_fft
  REAL(4), DIMENSION (:), ALLOCATABLE, SAVE :: xpb,ypb,zpb
  REAL(8), DIMENSION (:), ALLOCATABLE, SAVE :: xpc,ypc,zpc,vxpc,vypc&
       &,vzpc,vxpd,vypd,vzpd,spline_x,xpcc,ypcc,zpcc
  REAL(8), ALLOCATABLE, SAVE :: spline_y(:,:),RotMat(:,:,:)
  REAL(8), ALLOCATABLE, SAVE ::  wsave1(:),vacf_data(:)&
       &,vacf_spectra(:),rms_disp(:),tot_rms_disp(:)
  DOUBLE COMPLEX, ALLOCATABLE, SAVE ::  w1(:),w2(:)
CONTAINS
  SUBROUTINE get_memory_iobuffer(start_anl,stop_anl,start_time&
       &,end_time,length_run,atom_record,length_tot,length_fft)
    IMPLICIT none
    INTEGER :: start_anl,stop_anl,length_run,length_tot,length_fft&
         &,start_time,end_time,iret,atom_record
    CHARACTER(80) ::  errmsg
    CHARACTER(8) ::  string
    INTEGER  :: len,M_get_length
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
    
  END SUBROUTINE get_memory_iobuffer
    
END MODULE IOBUFFER_Mod
