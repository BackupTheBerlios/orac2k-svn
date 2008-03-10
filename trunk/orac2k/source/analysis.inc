#include "config.h"

      INTEGER  maxbin,hhisto_maxbin,maxcav_atom,maxcav_bin,maxcav_nres,
     &     maxcav

      PARAMETER(maxbin=_MAX_BIN_,
     &     hhisto_maxbin=_MAX_HHISTO_BIN_,
     &     maxcav_atom=_MAX_CAVITIES_ATOM_,
     &     maxcav_bin=_MAX_CAVITIES_BIN_,
     &     maxcav_nres=_MAX_CAVITIES_NRES_,
     &     maxcav=_MAX_CAVITIES_)

      INTEGER(4), POINTER :: cavity_n(:,:)
      REAL(8), POINTER :: cavity_h(:,:),cavity_r(:)
      COMMON /cavity/ cavity_n,cavity_h,cavity_r


