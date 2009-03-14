
#include "config.h"


!==== Used only to calculate the residence time ===============

      INTEGER vic,dlst,dtm
      PARAMETER (dlst=_DLST_, dtm=_DTM_, vic=_VIC_)
      INTEGER*2 wnlst(dlst,dtm),wlst(vic,dtm)
 
      COMMON /resid/  wnlst,wlst

