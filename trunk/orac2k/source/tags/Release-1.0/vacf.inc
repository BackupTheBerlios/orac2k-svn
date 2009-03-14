
#include "config.h"


      INCLUDE 'parst.h'
      INTEGER nb,nc
      REAL*4 vpl
      REAL*4 vafp,weip,vacf_tot,vacf_slt,vacf_slv,psp_tot
     &     ,psp_slv,psp_slt
      PARAMETER (nb=_NAT_SOLU_,nc=_NUM_CORR_)
      COMMON /velox/ weip(0:nc,2),vafp(0:nc,2),vpl(3,nb,nc)
     &     ,vacf_slt(0:nc),vacf_slv(0:nc),psp_slv(0:nc),psp_slt(0:nc)
     &     ,vacf_tot(0:nc),psp_tot(0:nc)

