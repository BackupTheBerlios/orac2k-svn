#define _REAL_ REAL(4)
*------------------ PARAMETER  STATEMENT -------------------------------

      INTEGER(4) ::  nbuffer,buffer_time,buffer_fft
      COMMON /iobuff/ nbuffer,buffer_time,buffer_fft
      _REAL_, DIMENSION (:), POINTER :: xpb,ypb,zpb
      REAL(8), DIMENSION (:), POINTER :: xpc,ypc,zpc,vxpc,vypc,vzpc,vxpd
     &     ,vypd,vzpd,spline_x,xpcc,ypcc,zpcc
      REAL(8), POINTER :: spline_y(:,:),RotMat(:,:,:)
      REAL(8), POINTER ::  wsave1(:),vacf_data(:),vacf_spectra(:)
     &     ,rms_disp(:),tot_rms_disp(:)
      DOUBLE COMPLEX, POINTER ::  w1(:),w2(:)
