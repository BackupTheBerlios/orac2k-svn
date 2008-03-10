#include "config.h"
      INTEGER max_neigh,pig_nnl,pnnlpp_vor,maxpla,maxver

      PARAMETER (max_neigh=_MAX_NEIGH_,
     &          pig_nnl=_MAX_IG_NNL_,
     &          pnnlpp_vor=MAXPLA+1,
     &          maxpla=MAXPLA,
     &          maxver=MAXVER)
      
      INTEGER, DIMENSION (:,:), POINTER :: nnlpp_vor,ig_nnl
      REAL(8), POINTER :: area_vor(:,:), volume_vor(:)
      COMMON /voronoi/ ig_nnl,nnlpp_vor,area_vor,volume_vor
