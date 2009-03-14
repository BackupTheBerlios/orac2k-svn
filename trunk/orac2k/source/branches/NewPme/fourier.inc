
#include "config.h"


!============================================================================
!---- Parameters of the sum over the reciprocal lattice ---------------------
!============================================================================

      INTEGER lmax,mmax,nmax
      INTEGER lofset,mofset,nofset,lmaxt,mmaxt,nmaxt,mmint,nmint

      PARAMETER (LMAX = _LMAX_, MMAX = _MMAX_, NMAX = _NMAX_, 
     x           lofset = 1, mofset = MMAX + 1, nofset = NMAX + 1,
     x           lmaxt = LMAX + 1, mmaxt = 2*MMAX + 1, 
     x           nmaxt = 2*NMAX + 1)

!============================================================================
